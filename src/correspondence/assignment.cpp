#include "vallescope2/correspondence/assignment.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace vallescope2 {
namespace {

struct GroupedAnchor {
    std::string sequence_id;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::string anchor_group_id;
    bool groupable = false;
};

struct Anchor {
    std::string anchor_id;
    std::string sample;
    std::string sequence_id;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::string anchor_group_id;
    std::string context_key;
    std::vector<std::uint32_t> forward_tokens;
    std::vector<std::uint64_t> alpha_prime;
};

struct Pair {
    std::uint32_t target_sample = 0;
    std::uint32_t query_sample = 0;
};

struct Candidate {
    const Anchor* target = nullptr;
    std::uint32_t edit_distance = 0;
    char strand = '.';
    std::int64_t score = 0;
    std::uint64_t q_alpha_prime = 0;
    std::uint64_t r_alpha_prime = 0;
};

struct DirectedPrimary {
    std::uint64_t pair_id = 0;
    std::uint32_t target_sample = 0;
    std::uint32_t query_sample = 0;
    std::string query_anchor_id;
    std::string target_anchor_id;
    std::string method;
    char strand = '.';
    std::uint32_t edit_distance = 0;
    std::int64_t score = 0;
};

struct CanonicalPrimary {
    const DirectedPrimary* record = nullptr;
    std::string anchor_a;
    std::string anchor_b;
    bool forward = true;
};

struct DataSet {
    std::vector<std::string> samples;
    std::unordered_map<std::string, std::uint32_t> sample_ids;
    std::unordered_map<std::string, std::string> sample_roles;
    std::vector<Anchor> anchors;
    std::vector<std::vector<const Anchor*>> anchors_by_sample;
    std::unordered_map<std::string, std::uint32_t> token_ids;
};

std::string pair_merge_mode_name(const PairMergeMode mode) {
    return mode == PairMergeMode::reciprocal ? "reciprocal" : "union";
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            fields.emplace_back(line.substr(begin));
            break;
        }
        fields.emplace_back(line.substr(begin, tab - begin));
        begin = tab + 1;
    }
    return fields;
}

std::uint64_t parse_u64(const std::string& text, const std::string& field) {
    try {
        std::size_t parsed = 0;
        const auto value = std::stoull(text, &parsed);
        if (text.empty() || parsed != text.size()) throw std::out_of_range("bad");
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer in " + field + ": " + text);
    }
}

std::unordered_map<std::string, std::size_t>
header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < header.size(); ++i) index.emplace(header[i], i);
    return index;
}

std::size_t require_column(
    const std::unordered_map<std::string, std::size_t>& columns,
    const std::string& name) {
    const auto found = columns.find(name);
    if (found == columns.end()) throw std::runtime_error("missing column: " + name);
    return found->second;
}

std::unordered_map<std::string, GroupedAnchor>
load_grouped_anchors(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open grouped anchors: " + path.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty grouped anchors file");
    const auto columns = header_index(split_tab(line));
    const auto anchor_id_column = require_column(columns, "anchor_id");
    const auto sequence_column = require_column(columns, "sequence_id");
    const auto start_column = require_column(columns, "start");
    const auto end_column = require_column(columns, "end");
    const auto group_column = require_column(columns, "anchor_group_id");
    const auto groupable_column = require_column(columns, "groupable");

    std::unordered_map<std::string, GroupedAnchor> anchors;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= groupable_column)
            throw std::runtime_error("invalid grouped anchor row");
        GroupedAnchor anchor;
        anchor.sequence_id = fields[sequence_column];
        anchor.start = parse_u64(fields[start_column], "start");
        anchor.end = parse_u64(fields[end_column], "end");
        anchor.anchor_group_id = fields[group_column];
        anchor.groupable = fields[groupable_column] == "true";
        anchors.emplace(fields[anchor_id_column], std::move(anchor));
    }
    return anchors;
}

std::uint32_t intern_token(DataSet& data, const std::string& token) {
    const auto found = data.token_ids.find(token);
    if (found != data.token_ids.end()) return found->second;
    const auto id = static_cast<std::uint32_t>(data.token_ids.size());
    data.token_ids.emplace(token, id);
    return id;
}

std::vector<std::uint32_t> parse_tokens(DataSet& data, const std::string& text) {
    std::vector<std::uint32_t> tokens;
    std::size_t begin = 0;
    while (begin <= text.size()) {
        const std::size_t separator = text.find('|', begin);
        const std::string token =
            separator == std::string::npos
                ? text.substr(begin)
                : text.substr(begin, separator - begin);
        if (!token.empty()) tokens.push_back(intern_token(data, token));
        if (separator == std::string::npos) break;
        begin = separator + 1;
    }
    return tokens;
}

void load_sequence_table(const std::filesystem::path& path, DataSet& data) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open sequence table: " + path.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty sequence table");
    const auto columns = header_index(split_tab(line));
    const auto sample_column = require_column(columns, "sample");
    const auto role_column = require_column(columns, "role");
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= role_column) throw std::runtime_error("invalid sequence table row");
        data.sample_roles.emplace(fields[sample_column], fields[role_column]);
    }
}

std::string strip_alpha_suffix(const std::string& header) {
    constexpr const char* suffix = "_alpha_prime";
    constexpr std::size_t suffix_size = 12;
    if (header.size() < suffix_size ||
        header.compare(header.size() - suffix_size, suffix_size, suffix) != 0) {
        throw std::runtime_error("invalid alpha prime column: " + header);
    }
    return header.substr(0, header.size() - suffix_size);
}

DataSet load_assignment_data(const std::filesystem::path& grouped_anchors,
                             const std::filesystem::path& anchor_contexts,
                             const std::filesystem::path& sequence_table) {
    DataSet data;
    load_sequence_table(sequence_table, data);
    const auto grouped = load_grouped_anchors(grouped_anchors);

    std::ifstream input(anchor_contexts);
    if (!input) throw std::runtime_error("cannot open anchor contexts: " + anchor_contexts.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty anchor contexts file");
    const auto header = split_tab(line);
    const auto columns = header_index(header);
    const auto anchor_id_column = require_column(columns, "anchor_id");
    const auto sequence_column = require_column(columns, "sequence_id");
    const auto sample_column = require_column(columns, "sample");
    const auto context_key_column = require_column(columns, "canonical_context_key");
    const auto forward_column = require_column(columns, "forward_context_tokens");

    std::vector<std::size_t> alpha_columns;
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i].size() >= 12 &&
            header[i].compare(header[i].size() - 12, 12, "_alpha_prime") == 0) {
            const auto sample = strip_alpha_suffix(header[i]);
            data.sample_ids.emplace(sample, static_cast<std::uint32_t>(data.samples.size()));
            data.samples.push_back(sample);
            alpha_columns.push_back(i);
        }
    }
    if (data.samples.empty()) throw std::runtime_error("anchor contexts have no alpha prime columns");
    data.anchors_by_sample.resize(data.samples.size());

    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= forward_column) throw std::runtime_error("invalid anchor context row");
        const auto grouped_found = grouped.find(fields[anchor_id_column]);
        if (grouped_found == grouped.end())
            throw std::runtime_error("anchor context has no grouped anchor: " + fields[anchor_id_column]);
        if (!grouped_found->second.groupable ||
            grouped_found->second.anchor_group_id == ".") {
            continue;
        }
        const auto sample_found = data.sample_ids.find(fields[sample_column]);
        if (sample_found == data.sample_ids.end())
            throw std::runtime_error("anchor context sample has no alpha column: " + fields[sample_column]);

        Anchor anchor;
        anchor.anchor_id = fields[anchor_id_column];
        anchor.sample = fields[sample_column];
        anchor.sequence_id = fields[sequence_column];
        anchor.start = grouped_found->second.start;
        anchor.end = grouped_found->second.end;
        anchor.anchor_group_id = grouped_found->second.anchor_group_id;
        anchor.context_key = fields[context_key_column];
        anchor.forward_tokens = parse_tokens(data, fields[forward_column]);
        anchor.alpha_prime.reserve(alpha_columns.size());
        for (const auto column : alpha_columns) {
            if (fields.size() <= column) throw std::runtime_error("missing alpha prime value");
            anchor.alpha_prime.push_back(parse_u64(fields[column], "alpha_prime"));
        }
        data.anchors.push_back(std::move(anchor));
        data.anchors_by_sample[sample_found->second].push_back(&data.anchors.back());
    }

    // Pointers into data.anchors may have been invalidated while the vector grew.
    for (auto& sample_anchors : data.anchors_by_sample) sample_anchors.clear();
    for (const auto& anchor : data.anchors) {
        const auto sample = data.sample_ids.at(anchor.sample);
        data.anchors_by_sample[sample].push_back(&anchor);
    }
    return data;
}

std::vector<Pair> make_ordered_pairs(const DataSet& data) {
    std::vector<std::uint32_t> targets;
    std::vector<std::uint32_t> queries;
    for (std::uint32_t i = 0; i < data.samples.size(); ++i) {
        const auto role_found = data.sample_roles.find(data.samples[i]);
        const auto role = role_found == data.sample_roles.end() ? "self" : role_found->second;
        if (role == "target") targets.push_back(i);
        else if (role == "query") queries.push_back(i);
    }
    std::vector<Pair> pairs;
    if (!targets.empty() || !queries.empty()) {
        for (const auto target : targets)
            for (const auto query : queries)
                pairs.push_back({target, query});
        return pairs;
    }
    if (data.samples.size() == 1) {
        pairs.push_back({0, 0});
        return pairs;
    }
    for (std::uint32_t target = 0; target < data.samples.size(); ++target) {
        for (std::uint32_t query = 0; query < data.samples.size(); ++query) {
            if (target != query) pairs.push_back({target, query});
        }
    }
    return pairs;
}

std::uint32_t bounded_edit_distance(const std::vector<std::uint32_t>& left,
                                    const std::vector<std::uint32_t>& right,
                                    const std::uint32_t limit) {
    const std::uint32_t fail = limit + 1;
    const auto n = static_cast<std::int64_t>(left.size());
    const auto m = static_cast<std::int64_t>(right.size());
    if (std::llabs(n - m) > static_cast<long long>(limit)) return fail;

    const std::uint32_t inf = fail;
    std::vector<std::uint32_t> previous(right.size() + 1, inf);
    std::vector<std::uint32_t> current(right.size() + 1, inf);
    const std::uint32_t first_end =
        std::min<std::uint32_t>(static_cast<std::uint32_t>(right.size()), limit);
    for (std::uint32_t j = 0; j <= first_end; ++j) previous[j] = j;

    for (std::uint32_t i = 1; i <= left.size(); ++i) {
        std::fill(current.begin(), current.end(), inf);
        const std::uint32_t begin = i > limit ? i - limit : 1;
        const std::uint32_t end =
            std::min<std::uint32_t>(static_cast<std::uint32_t>(right.size()), i + limit);
        if (begin == 1) current[0] = i <= limit ? i : inf;
        std::uint32_t row_min = inf;
        for (std::uint32_t j = begin; j <= end; ++j) {
            const std::uint32_t substitution =
                previous[j - 1] + (left[i - 1] == right[j - 1] ? 0 : 1);
            const std::uint32_t deletion = previous[j] + 1;
            const std::uint32_t insertion = current[j - 1] + 1;
            current[j] = std::min({substitution, deletion, insertion, inf});
            row_min = std::min(row_min, current[j]);
        }
        if (row_min > limit) return fail;
        previous.swap(current);
    }
    return previous[right.size()] <= limit ? previous[right.size()] : fail;
}

std::vector<std::uint32_t> reversed_tokens(const std::vector<std::uint32_t>& tokens) {
    return {tokens.rbegin(), tokens.rend()};
}

Candidate make_candidate(const Anchor& query,
                         const std::vector<std::uint32_t>& reverse_query,
                         const Anchor& target,
                         const std::uint64_t q_alpha_prime,
                         const std::uint64_t r_alpha_prime,
                         const std::uint32_t edit_distance_limit) {
    const auto forward_ed = bounded_edit_distance(
        query.forward_tokens, target.forward_tokens, edit_distance_limit);
    const auto reverse_ed = bounded_edit_distance(
        reverse_query, target.forward_tokens, edit_distance_limit);
    const bool use_reverse = reverse_ed < forward_ed;
    const auto edit_distance = use_reverse ? reverse_ed : forward_ed;
    Candidate candidate;
    candidate.target = &target;
    candidate.edit_distance = edit_distance;
    candidate.strand = use_reverse ? '-' : '+';
    candidate.q_alpha_prime = q_alpha_prime;
    candidate.r_alpha_prime = r_alpha_prime;
    candidate.score = static_cast<std::int64_t>(candidate.q_alpha_prime) +
                      static_cast<std::int64_t>(candidate.r_alpha_prime) -
                      static_cast<std::int64_t>(edit_distance);
    return candidate;
}

bool candidate_less(const Candidate& left, const Candidate& right) {
    if (left.score != right.score) return left.score > right.score;
    if (left.edit_distance != right.edit_distance)
        return left.edit_distance < right.edit_distance;
    return left.target->anchor_id < right.target->anchor_id;
}

void write_unmatched(std::ofstream& output,
                     const std::uint64_t pair_id,
                     const std::string& target_sample,
                     const std::string& query_sample,
                     const Anchor& query) {
    output << pair_id << '\t' << target_sample << '\t' << query_sample << '\t'
           << query.anchor_id << "\t.\tunmatched\tunmatched\t.\t.\t.\t.\t0\t.\t.\t"
           << query.anchor_group_id << "\t.\t" << query.context_key << "\t.\n";
}

void write_assignment(std::ofstream& output,
                      const std::uint64_t pair_id,
                      const std::string& target_sample,
                      const std::string& query_sample,
                      const Anchor& query,
                      const Candidate& candidate,
                      const std::string& status,
                      const std::string& method,
                      const std::int64_t score_margin,
                      const std::uint64_t candidate_count) {
    output << pair_id << '\t' << target_sample << '\t' << query_sample << '\t'
           << query.anchor_id << '\t' << candidate.target->anchor_id << '\t'
           << status << '\t' << method << '\t' << candidate.strand << '\t'
           << candidate.edit_distance << '\t' << candidate.score << '\t'
           << score_margin << '\t' << candidate_count << '\t'
           << candidate.q_alpha_prime << '\t' << candidate.r_alpha_prime << '\t'
           << query.anchor_group_id << '\t' << candidate.target->anchor_group_id
           << '\t' << query.context_key << '\t' << candidate.target->context_key
           << '\n';
}

std::string merge_key(const std::uint32_t sample_a,
                      const std::uint32_t sample_b,
                      const std::string& anchor_a,
                      const std::string& anchor_b) {
    constexpr char separator = '\x1f';
    return std::to_string(sample_a) + separator + std::to_string(sample_b) +
           separator + anchor_a + separator + anchor_b;
}

std::string anchor_key(const std::uint32_t sample_a,
                       const std::uint32_t sample_b,
                       const std::string& anchor) {
    constexpr char separator = '\x1f';
    return std::to_string(sample_a) + separator + std::to_string(sample_b) +
           separator + anchor;
}

CanonicalPrimary canonical_primary(const DirectedPrimary& primary) {
    CanonicalPrimary item;
    item.record = &primary;
    if (primary.target_sample < primary.query_sample) {
        item.anchor_a = primary.target_anchor_id;
        item.anchor_b = primary.query_anchor_id;
        item.forward = true;
    } else {
        item.anchor_a = primary.query_anchor_id;
        item.anchor_b = primary.target_anchor_id;
        item.forward = false;
    }
    return item;
}

void write_correspondence_metadata(const std::filesystem::path& path,
                                   const AssignmentParameters& parameters,
                                   const AssignmentResult& result) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create correspondence metadata");
    output << "{\n"
           << "  \"pair_merge_mode\": \""
           << pair_merge_mode_name(parameters.pair_merge_mode) << "\",\n"
           << "  \"correspondence_count\": " << result.correspondence_count << ",\n"
           << "  \"reciprocal_correspondence_count\": "
           << result.reciprocal_correspondence_count << ",\n"
           << "  \"forward_only_correspondence_count\": "
           << result.forward_only_correspondence_count << ",\n"
           << "  \"reverse_only_correspondence_count\": "
           << result.reverse_only_correspondence_count << ",\n"
           << "  \"discordant_correspondence_count\": "
           << result.discordant_correspondence_count << "\n"
           << "}\n";
}

void write_correspondences(const std::filesystem::path& output_path,
                           const std::filesystem::path& metadata_path,
                           const DataSet& data,
                           const std::vector<DirectedPrimary>& primaries,
                           const AssignmentParameters& parameters,
                           AssignmentResult& result) {
    std::unordered_map<std::string, const DirectedPrimary*> forward;
    std::unordered_map<std::string, const DirectedPrimary*> reverse;
    std::unordered_map<std::string, std::vector<std::string>> forward_by_anchor_a;
    std::unordered_map<std::string, std::vector<std::string>> reverse_by_anchor_b;

    for (const auto& primary : primaries) {
        if (primary.target_sample == primary.query_sample) continue;
        const auto sample_a = std::min(primary.target_sample, primary.query_sample);
        const auto sample_b = std::max(primary.target_sample, primary.query_sample);
        const auto canonical = canonical_primary(primary);
        const auto key = merge_key(sample_a, sample_b, canonical.anchor_a,
                                   canonical.anchor_b);
        if (canonical.forward) {
            forward.emplace(key, &primary);
            forward_by_anchor_a[anchor_key(sample_a, sample_b, canonical.anchor_a)]
                .push_back(canonical.anchor_b);
        } else {
            reverse.emplace(key, &primary);
            reverse_by_anchor_b[anchor_key(sample_a, sample_b, canonical.anchor_b)]
                .push_back(canonical.anchor_a);
        }
    }

    std::vector<std::string> keys;
    keys.reserve(forward.size() + reverse.size());
    for (const auto& item : forward) keys.push_back(item.first);
    for (const auto& item : reverse) {
        if (forward.find(item.first) == forward.end()) keys.push_back(item.first);
    }
    std::sort(keys.begin(), keys.end());

    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("cannot create correspondence output");
    output << "sample_a\tsample_b\tanchor_a\tanchor_b\tsupport_direction"
              "\tforward_pair_id\treverse_pair_id\tforward_score\treverse_score"
              "\tforward_edit_distance\treverse_edit_distance"
              "\tforward_strand\treverse_strand"
              "\tforward_method\treverse_method\tdiscordant_anchor\n";

    for (const auto& key : keys) {
        const auto forward_found = forward.find(key);
        const auto reverse_found = reverse.find(key);
        const auto* f = forward_found == forward.end() ? nullptr : forward_found->second;
        const auto* r = reverse_found == reverse.end() ? nullptr : reverse_found->second;
        const auto& primary = f ? *f : *r;
        const auto sample_a = std::min(primary.target_sample, primary.query_sample);
        const auto sample_b = std::max(primary.target_sample, primary.query_sample);
        const auto canonical = canonical_primary(primary);

        std::string support_direction;
        std::string discordant_anchor = ".";
        if (f && r) {
            support_direction = "both";
            ++result.reciprocal_correspondence_count;
        } else if (f) {
            const auto opposite = reverse_by_anchor_b.find(
                anchor_key(sample_a, sample_b, canonical.anchor_b));
            if (opposite != reverse_by_anchor_b.end() &&
                std::find(opposite->second.begin(), opposite->second.end(),
                          canonical.anchor_a) == opposite->second.end()) {
                support_direction = "discordant";
                discordant_anchor = opposite->second.front();
                ++result.discordant_correspondence_count;
            } else {
                support_direction = "forward_only";
                ++result.forward_only_correspondence_count;
            }
        } else {
            const auto opposite = forward_by_anchor_a.find(
                anchor_key(sample_a, sample_b, canonical.anchor_a));
            if (opposite != forward_by_anchor_a.end() &&
                std::find(opposite->second.begin(), opposite->second.end(),
                          canonical.anchor_b) == opposite->second.end()) {
                support_direction = "discordant";
                discordant_anchor = opposite->second.front();
                ++result.discordant_correspondence_count;
            } else {
                support_direction = "reverse_only";
                ++result.reverse_only_correspondence_count;
            }
        }

        if (parameters.pair_merge_mode == PairMergeMode::reciprocal &&
            support_direction != "both") {
            if (support_direction == "forward_only")
                --result.forward_only_correspondence_count;
            else if (support_direction == "reverse_only")
                --result.reverse_only_correspondence_count;
            else if (support_direction == "discordant")
                --result.discordant_correspondence_count;
            continue;
        }

        output << data.samples[sample_a] << '\t' << data.samples[sample_b] << '\t'
               << canonical.anchor_a << '\t' << canonical.anchor_b << '\t'
               << support_direction << '\t'
               << (f ? std::to_string(f->pair_id) : ".") << '\t'
               << (r ? std::to_string(r->pair_id) : ".") << '\t'
               << (f ? std::to_string(f->score) : ".") << '\t'
               << (r ? std::to_string(r->score) : ".") << '\t'
               << (f ? std::to_string(f->edit_distance) : ".") << '\t'
               << (r ? std::to_string(r->edit_distance) : ".") << '\t'
               << (f ? std::string(1, f->strand) : ".") << '\t'
               << (r ? std::string(1, r->strand) : ".") << '\t'
               << (f ? f->method : ".") << '\t'
               << (r ? r->method : ".") << '\t'
               << discordant_anchor << '\n';
        ++result.correspondence_count;
    }

    write_correspondence_metadata(metadata_path, parameters, result);
}

void write_metadata(const std::filesystem::path& path,
                    const AssignmentParameters& parameters,
                    const AssignmentResult& result) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create assignment metadata");
    output << "{\n"
           << "  \"legacy_beta_tolerance\": " << parameters.beta_tolerance << ",\n"
           << "  \"candidate_filter\": \"candidate_score > min_candidate_score\",\n"
           << "  \"min_candidate_score\": " << parameters.min_candidate_score << ",\n"
           << "  \"pair_merge_mode\": \""
           << pair_merge_mode_name(parameters.pair_merge_mode) << "\",\n"
           << "  \"primary_margin\": " << parameters.primary_margin << ",\n"
           << "  \"ordered_pair_count\": " << result.ordered_pair_count << ",\n"
           << "  \"query_anchor_count\": " << result.query_anchor_count << ",\n"
           << "  \"primary_count\": " << result.primary_count << ",\n"
           << "  \"ambiguous_query_count\": " << result.ambiguous_query_count << ",\n"
           << "  \"ambiguous_row_count\": " << result.ambiguous_row_count << ",\n"
           << "  \"unmatched_count\": " << result.unmatched_count << ",\n"
           << "  \"exact_context_primary_count\": "
           << result.exact_context_primary_count << ",\n"
           << "  \"edit_distance_primary_count\": "
           << result.edit_distance_primary_count << "\n"
           << "}\n";
}

}  // namespace

AssignmentResult assign_context_correspondences(
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& anchor_contexts,
    const std::filesystem::path& sequence_table,
    const std::filesystem::path& assignment_output,
    const std::filesystem::path& correspondence_output,
    const std::filesystem::path& correspondence_metadata_output,
    const std::filesystem::path& metadata_output,
    const AssignmentParameters& parameters) {
    const auto data = load_assignment_data(grouped_anchors, anchor_contexts, sequence_table);
    const auto pairs = make_ordered_pairs(data);

    std::ofstream output(assignment_output);
    if (!output) throw std::runtime_error("cannot create assignment output");
    output << "pair_id\ttarget_sample\tquery_sample\tquery_anchor_id"
              "\ttarget_anchor_id\tstatus\tmethod\tassign_strand"
              "\tedit_distance\tcandidate_score\tscore_margin\tn_candidates"
              "\tq_alpha_prime\tr_alpha_prime\tquery_anchor_group_id"
              "\ttarget_anchor_group_id\tquery_context_key\ttarget_context_key\n";

    AssignmentResult result;
    result.ordered_pair_count = pairs.size();
    std::vector<DirectedPrimary> primaries;

    for (std::uint64_t pair_index = 0; pair_index < pairs.size(); ++pair_index) {
        const auto pair = pairs[pair_index];
        const auto& target_sample_name = data.samples[pair.target_sample];
        const auto& query_sample_name = data.samples[pair.query_sample];

        std::unordered_map<std::string, std::vector<const Anchor*>> by_context;
        std::unordered_map<std::string, std::vector<const Anchor*>> by_group;
        for (const auto* target : data.anchors_by_sample[pair.target_sample]) {
            by_context[target->context_key].push_back(target);
            by_group[target->anchor_group_id].push_back(target);
        }

        for (const auto* query_ptr : data.anchors_by_sample[pair.query_sample]) {
            const auto& query = *query_ptr;
            ++result.query_anchor_count;

            const auto exact_found = by_context.find(query.context_key);
            if (exact_found != by_context.end() && exact_found->second.size() == 1) {
                const auto* target = exact_found->second.front();
                if (target->anchor_id == query.anchor_id) {
                    write_unmatched(output, pair_index, target_sample_name,
                                    query_sample_name, query);
                    ++result.unmatched_count;
                    continue;
                }
                Candidate candidate;
                candidate.target = target;
                candidate.edit_distance = 0;
                candidate.strand = '.';
                candidate.q_alpha_prime = query.alpha_prime[pair.target_sample];
                candidate.r_alpha_prime = target->alpha_prime[pair.query_sample];
                candidate.score = static_cast<std::int64_t>(candidate.q_alpha_prime) +
                                  static_cast<std::int64_t>(candidate.r_alpha_prime);
                if (candidate.q_alpha_prime > 0 || candidate.r_alpha_prime > 0) {
                    write_assignment(output, pair_index, target_sample_name,
                                     query_sample_name, query, candidate,
                                     "primary", "exact_context", 0, 1);
                    primaries.push_back(
                        {pair_index, pair.target_sample, pair.query_sample,
                         query.anchor_id, target->anchor_id, "exact_context",
                         candidate.strand, candidate.edit_distance,
                         candidate.score});
                    ++result.primary_count;
                    ++result.exact_context_primary_count;
                    continue;
                }
            }

            std::vector<Candidate> candidates;
            const auto q_alpha_prime = query.alpha_prime[pair.target_sample];
            const auto group_found = by_group.find(query.anchor_group_id);
            if (q_alpha_prime > 0 && group_found != by_group.end()) {
                const auto reverse_query = reversed_tokens(query.forward_tokens);
                for (const auto* target : group_found->second) {
                    if (target->anchor_id == query.anchor_id) continue;
                    const auto r_alpha_prime = target->alpha_prime[pair.query_sample];
                    if (r_alpha_prime == 0) continue;
                    const auto maximum_passing_edit_distance =
                        static_cast<std::int64_t>(q_alpha_prime) +
                        static_cast<std::int64_t>(r_alpha_prime) -
                        parameters.min_candidate_score - 1;
                    if (maximum_passing_edit_distance < 0) continue;
                    const auto maximum_possible_edit_distance =
                        std::max(query.forward_tokens.size(),
                                 target->forward_tokens.size());
                    const auto edit_distance_limit =
                        static_cast<std::uint32_t>(std::min<std::uint64_t>(
                            static_cast<std::uint64_t>(
                                maximum_passing_edit_distance),
                            maximum_possible_edit_distance));
                    auto candidate = make_candidate(
                        query, reverse_query, *target, q_alpha_prime,
                        r_alpha_prime, edit_distance_limit);
                    if (candidate.score > parameters.min_candidate_score) {
                        candidates.push_back(candidate);
                    }
                }
            }

            if (candidates.empty()) {
                write_unmatched(output, pair_index, target_sample_name,
                                query_sample_name, query);
                ++result.unmatched_count;
                continue;
            }
            std::sort(candidates.begin(), candidates.end(), candidate_less);
            if (candidates.size() == 1) {
                write_assignment(output, pair_index, target_sample_name,
                                 query_sample_name, query, candidates.front(),
                                 "primary", "edit_distance", 0, 1);
                primaries.push_back(
                    {pair_index, pair.target_sample, pair.query_sample,
                     query.anchor_id, candidates.front().target->anchor_id,
                     "edit_distance", candidates.front().strand,
                     candidates.front().edit_distance, candidates.front().score});
                ++result.primary_count;
                ++result.edit_distance_primary_count;
                continue;
            }

            const auto margin = candidates[0].score - candidates[1].score;
            if (margin >= parameters.primary_margin) {
                write_assignment(output, pair_index, target_sample_name,
                                 query_sample_name, query, candidates.front(),
                                 "primary", "edit_distance", margin,
                                 candidates.size());
                primaries.push_back(
                    {pair_index, pair.target_sample, pair.query_sample,
                     query.anchor_id, candidates.front().target->anchor_id,
                     "edit_distance", candidates.front().strand,
                     candidates.front().edit_distance, candidates.front().score});
                ++result.primary_count;
                ++result.edit_distance_primary_count;
            } else {
                ++result.ambiguous_query_count;
                for (const auto& candidate : candidates) {
                    if (candidates.front().score - candidate.score >
                        parameters.primary_margin) {
                        break;
                    }
                    write_assignment(output, pair_index, target_sample_name,
                                     query_sample_name, query, candidate,
                                     "ambiguous", "edit_distance", margin,
                                     candidates.size());
                    ++result.ambiguous_row_count;
                }
            }
        }
    }

    write_correspondences(correspondence_output, correspondence_metadata_output,
                          data, primaries, parameters, result);
    write_metadata(metadata_output, parameters, result);
    return result;
}

}  // namespace vallescope2
