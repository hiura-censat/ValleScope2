#include "vallescope2/bundle/chaining.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
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

struct AnchorPosition {
    std::string sequence_id;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::uint64_t center = 0;
};

struct Candidate {
    std::string sample_a;
    std::string sample_b;
    std::string anchor_a;
    std::string anchor_b;
    std::string support_direction;
    std::string sequence_a;
    std::string sequence_b;
    std::uint64_t ref_center = 0;
    std::uint64_t query_center = 0;
    double score = 0.0;
    char strand = '+';
};

struct ChainState {
    double dp = 0.0;
    int predecessor = -1;
};

struct EmittedChain {
    std::uint64_t id = 0;
    std::string sample_a;
    std::string sample_b;
    std::string sequence_a;
    std::string sequence_b;
    char strand = '+';
    std::uint64_t ref_start = 0;
    std::uint64_t ref_end = 0;
    std::uint64_t query_start = 0;
    std::uint64_t query_end = 0;
    double score = 0.0;
    std::vector<Candidate> anchors;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(begin));
            break;
        }
        fields.push_back(line.substr(begin, tab - begin));
        begin = tab + 1;
    }
    return fields;
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

std::uint64_t parse_u64(const std::string& value, const std::string& column) {
    try {
        std::size_t parsed = 0;
        const auto number = std::stoull(value, &parsed);
        if (value.empty() || parsed != value.size()) throw std::out_of_range("bad");
        return number;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer in " + column + ": " + value);
    }
}

double parse_optional_score(const std::string& value,
                            const double missing_value) {
    if (value == ".") return missing_value;
    try {
        std::size_t parsed = 0;
        const double number = std::stod(value, &parsed);
        if (value.empty() || parsed != value.size()) throw std::out_of_range("bad");
        return number;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid score: " + value);
    }
}

double parse_required_score(const std::string& value,
                            const std::string& column) {
    if (value == ".") throw std::runtime_error("missing score in " + column);
    try {
        std::size_t parsed = 0;
        const double number = std::stod(value, &parsed);
        if (value.empty() || parsed != value.size()) throw std::out_of_range("bad");
        return number;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid score in " + column + ": " + value);
    }
}

char parse_optional_strand(const std::string& value) {
    if (value == "+" || value == "-") return value.front();
    return '+';
}

std::unordered_map<std::string, AnchorPosition>
load_anchor_positions(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open grouped anchors: " + path.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty grouped anchors");
    const auto columns = header_index(split_tab(line));
    const auto anchor_column = require_column(columns, "anchor_id");
    const auto sequence_column = require_column(columns, "sequence_id");
    const auto start_column = require_column(columns, "start");
    const auto end_column = require_column(columns, "end");

    std::unordered_map<std::string, AnchorPosition> positions;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= end_column) throw std::runtime_error("invalid grouped anchor row");
        AnchorPosition position;
        position.sequence_id = fields[sequence_column];
        position.start = parse_u64(fields[start_column], "start");
        position.end = parse_u64(fields[end_column], "end");
        position.center = (position.start + position.end) / 2;
        positions.emplace(fields[anchor_column], std::move(position));
    }
    return positions;
}

std::pair<double, char> choose_score_and_strand(
    const std::string& support_direction,
    const std::string& forward_score_text,
    const std::string& reverse_score_text,
    const std::string& forward_strand_text,
    const std::string& reverse_strand_text) {
    const double missing = -std::numeric_limits<double>::infinity();
    const double forward_score = parse_optional_score(forward_score_text, missing);
    const double reverse_score = parse_optional_score(reverse_score_text, missing);
    const char forward_strand = parse_optional_strand(forward_strand_text);
    const char reverse_strand = parse_optional_strand(reverse_strand_text);

    if (support_direction == "forward_only") return {forward_score, forward_strand};
    if (support_direction == "reverse_only") return {reverse_score, reverse_strand};
    if (reverse_score > forward_score) return {reverse_score, reverse_strand};
    return {forward_score, forward_strand};
}

std::vector<Candidate> load_candidates(
    const std::filesystem::path& path,
    const std::unordered_map<std::string, AnchorPosition>& positions) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open correspondences: " + path.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty correspondences file");
    const auto columns = header_index(split_tab(line));
    const auto sample_a_column = require_column(columns, "sample_a");
    const auto sample_b_column = require_column(columns, "sample_b");
    const auto anchor_a_column = require_column(columns, "anchor_a");
    const auto anchor_b_column = require_column(columns, "anchor_b");
    const auto support_column = require_column(columns, "support_direction");
    const auto forward_score_column = require_column(columns, "forward_score");
    const auto reverse_score_column = require_column(columns, "reverse_score");
    const auto forward_strand_column = require_column(columns, "forward_strand");
    const auto reverse_strand_column = require_column(columns, "reverse_strand");

    std::vector<Candidate> candidates;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= reverse_strand_column)
            throw std::runtime_error("invalid correspondence row");
        const auto anchor_a = fields[anchor_a_column];
        const auto anchor_b = fields[anchor_b_column];
        const auto position_a = positions.find(anchor_a);
        const auto position_b = positions.find(anchor_b);
        if (position_a == positions.end() || position_b == positions.end()) {
            throw std::runtime_error("correspondence refers to unknown anchor");
        }
        const auto score_and_strand = choose_score_and_strand(
            fields[support_column], fields[forward_score_column],
            fields[reverse_score_column], fields[forward_strand_column],
            fields[reverse_strand_column]);
        if (!std::isfinite(score_and_strand.first)) continue;
        Candidate candidate;
        candidate.sample_a = fields[sample_a_column];
        candidate.sample_b = fields[sample_b_column];
        candidate.anchor_a = anchor_a;
        candidate.anchor_b = anchor_b;
        candidate.support_direction = fields[support_column];
        candidate.sequence_a = position_a->second.sequence_id;
        candidate.sequence_b = position_b->second.sequence_id;
        candidate.ref_center = position_a->second.center;
        candidate.query_center = position_b->second.center;
        candidate.score = score_and_strand.first;
        candidate.strand = score_and_strand.second == '-' ? '-' : '+';
        candidates.push_back(std::move(candidate));
    }
    return candidates;
}

std::vector<Candidate> load_assignment_candidates(
    const std::filesystem::path& path,
    const std::unordered_map<std::string, AnchorPosition>& positions) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open assignments: " + path.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty assignments file");
    const auto columns = header_index(split_tab(line));
    const auto target_sample_column = require_column(columns, "target_sample");
    const auto query_sample_column = require_column(columns, "query_sample");
    const auto query_anchor_column = require_column(columns, "query_anchor_id");
    const auto target_anchor_column = require_column(columns, "target_anchor_id");
    const auto status_column = require_column(columns, "status");
    const auto strand_column = require_column(columns, "assign_strand");
    const auto score_column = require_column(columns, "candidate_score");

    std::vector<Candidate> candidates;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= score_column)
            throw std::runtime_error("invalid assignment row");
        const auto& status = fields[status_column];
        if (status != "primary" && status != "ambiguous") continue;
        const auto& target_anchor = fields[target_anchor_column];
        const auto& query_anchor = fields[query_anchor_column];
        if (target_anchor == "." || query_anchor == ".") continue;
        const auto position_a = positions.find(target_anchor);
        const auto position_b = positions.find(query_anchor);
        if (position_a == positions.end() || position_b == positions.end()) {
            throw std::runtime_error("assignment refers to unknown anchor");
        }

        Candidate candidate;
        candidate.sample_a = fields[target_sample_column];
        candidate.sample_b = fields[query_sample_column];
        candidate.anchor_a = target_anchor;
        candidate.anchor_b = query_anchor;
        candidate.support_direction = status;
        candidate.sequence_a = position_a->second.sequence_id;
        candidate.sequence_b = position_b->second.sequence_id;
        candidate.ref_center = position_a->second.center;
        candidate.query_center = position_b->second.center;
        candidate.score = parse_required_score(fields[score_column],
                                               "candidate_score");
        candidate.strand = fields[strand_column] == "-" ? '-' : '+';
        candidates.push_back(std::move(candidate));
    }
    return candidates;
}

std::string group_key(const Candidate& candidate) {
    return candidate.sample_a + '\x1f' + candidate.sample_b + '\x1f' +
           candidate.sequence_a + '\x1f' + candidate.sequence_b + '\x1f' +
           candidate.strand;
}

std::string candidate_key(const Candidate& candidate) {
    return candidate.anchor_a + '\x1f' + candidate.anchor_b;
}

bool predecessor_allowed(const Candidate& predecessor,
                         const Candidate& current) {
    if (predecessor.ref_center >= current.ref_center) return false;
    if (current.strand == '+') {
        return predecessor.query_center < current.query_center;
    }
    return predecessor.query_center > current.query_center;
}

double gap_cost(const Candidate& predecessor,
                const Candidate& current,
                const ChainingParameters& parameters) {
    const auto dr = current.ref_center - predecessor.ref_center;
    const auto dq = current.query_center > predecessor.query_center
                        ? current.query_center - predecessor.query_center
                        : predecessor.query_center - current.query_center;
    const double normalized =
        static_cast<double>(dr > dq ? dr - dq : dq - dr) /
        static_cast<double>(parameters.gap_unit);
    return parameters.gap_weight * std::log1p(normalized);
}

std::vector<std::size_t> chain_group(std::vector<Candidate>& candidates,
                                     const ChainingParameters& parameters,
                                     double& best_score) {
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& left, const Candidate& right) {
                  if (left.ref_center != right.ref_center)
                      return left.ref_center < right.ref_center;
                  if (left.strand == '-')
                      return left.query_center > right.query_center;
                  return left.query_center < right.query_center;
              });
    std::vector<ChainState> states(candidates.size());
    best_score = -std::numeric_limits<double>::infinity();
    int best_index = -1;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        states[i].dp = candidates[i].score;
        const std::size_t begin =
            i > parameters.predecessor_count ? i - parameters.predecessor_count : 0;
        for (std::size_t j = i; j-- > begin;) {
            if (!predecessor_allowed(candidates[j], candidates[i])) continue;
            const double score =
                candidates[i].score + states[j].dp -
                gap_cost(candidates[j], candidates[i], parameters);
            if (score > states[i].dp) {
                states[i].dp = score;
                states[i].predecessor = static_cast<int>(j);
            }
        }
        if (states[i].dp > best_score) {
            best_score = states[i].dp;
            best_index = static_cast<int>(i);
        }
    }

    std::vector<std::size_t> chain;
    for (int index = best_index; index >= 0; index = states[index].predecessor) {
        chain.push_back(static_cast<std::size_t>(index));
    }
    std::reverse(chain.begin(), chain.end());
    return chain;
}

EmittedChain make_emitted_chain(const std::uint64_t chain_id,
                                const std::vector<Candidate>& group,
                                const std::vector<std::size_t>& chain,
                                const double chain_score) {
    const auto& first = group[chain.front()];
    const auto& last = group[chain.back()];
    EmittedChain emitted;
    emitted.id = chain_id;
    emitted.sample_a = first.sample_a;
    emitted.sample_b = first.sample_b;
    emitted.sequence_a = first.sequence_a;
    emitted.sequence_b = first.sequence_b;
    emitted.strand = first.strand;
    emitted.ref_start = first.ref_center;
    emitted.ref_end = last.ref_center;
    emitted.query_start = first.query_center;
    emitted.query_end = first.query_center;
    emitted.score = chain_score;
    emitted.anchors.reserve(chain.size());
    for (const auto index : chain) {
        emitted.query_start = std::min(emitted.query_start, group[index].query_center);
        emitted.query_end = std::max(emitted.query_end, group[index].query_center);
        emitted.anchors.push_back(group[index]);
    }
    return emitted;
}

void write_chain(std::ostream& chains,
                 std::ostream& chain_anchors,
                 const EmittedChain& chain) {
    chains << chain.id << '\t' << chain.sample_a << '\t'
           << chain.sample_b << '\t' << chain.sequence_a << '\t'
           << chain.sequence_b << '\t' << chain.strand << '\t'
           << chain.anchors.size() << '\t' << chain.score << '\t'
           << chain.ref_start << '\t' << chain.ref_end << '\t'
           << chain.query_start << '\t' << chain.query_end << '\n';
    std::uint64_t rank = 0;
    for (const auto& candidate : chain.anchors) {
        chain_anchors << chain.id << '\t' << rank++ << '\t'
                      << candidate.sample_a << '\t' << candidate.sample_b << '\t'
                      << candidate.anchor_a << '\t' << candidate.anchor_b << '\t'
                      << candidate.ref_center << '\t' << candidate.query_center << '\t'
                      << candidate.strand << '\t' << candidate.score << '\t'
                      << candidate.support_direction << '\n';
    }
}

void write_refined_chain(std::ostream& chains,
                         std::ostream& chain_anchors,
                         const EmittedChain& chain,
                         const std::string& refinement_type,
                         const std::string& parent_chain_id,
                         const std::string& left_chain_id,
                         const std::string& right_chain_id) {
    chains << chain.id << '\t' << refinement_type << '\t'
           << parent_chain_id << '\t' << left_chain_id << '\t'
           << right_chain_id << '\t' << chain.sample_a << '\t'
           << chain.sample_b << '\t' << chain.sequence_a << '\t'
           << chain.sequence_b << '\t' << chain.strand << '\t'
           << chain.anchors.size() << '\t' << chain.score << '\t'
           << chain.ref_start << '\t' << chain.ref_end << '\t'
           << chain.query_start << '\t' << chain.query_end << '\n';
    std::uint64_t rank = 0;
    for (const auto& candidate : chain.anchors) {
        chain_anchors << chain.id << '\t' << rank++ << '\t'
                      << refinement_type << '\t' << candidate.sample_a << '\t'
                      << candidate.sample_b << '\t' << candidate.anchor_a << '\t'
                      << candidate.anchor_b << '\t' << candidate.ref_center << '\t'
                      << candidate.query_center << '\t' << candidate.strand << '\t'
                      << candidate.score << '\t' << candidate.support_direction << '\n';
    }
}

bool same_chain_group(const EmittedChain& chain,
                      const std::string& key) {
    Candidate candidate;
    candidate.sample_a = chain.sample_a;
    candidate.sample_b = chain.sample_b;
    candidate.sequence_a = chain.sequence_a;
    candidate.sequence_b = chain.sequence_b;
    candidate.strand = chain.strand;
    return group_key(candidate) == key;
}

std::uint64_t subtract_window(const std::uint64_t value,
                              const std::uint32_t window) {
    return value > window ? value - window : 0;
}

bool candidate_in_box(const Candidate& candidate,
                      const std::uint64_t ref_start,
                      const std::uint64_t ref_end,
                      const std::uint64_t query_start,
                      const std::uint64_t query_end,
                      const std::unordered_set<std::string>& used) {
    if (used.find(candidate_key(candidate)) != used.end()) return false;
    return ref_start <= candidate.ref_center && candidate.ref_center <= ref_end &&
           query_start <= candidate.query_center && candidate.query_center <= query_end;
}

bool emit_refined_subset(
    std::vector<Candidate> subset,
    const ChainingParameters& parameters,
    std::uint64_t& next_chain_id,
    ChainingResult& result,
    std::unordered_set<std::string>& used,
    std::ostream& chains,
    std::ostream& chain_anchors,
    const std::string& refinement_type,
    const std::string& parent_chain_id,
    const std::string& left_chain_id,
    const std::string& right_chain_id) {
    if (subset.size() < parameters.refinement_min_chain_anchors) return false;
    double chain_score = 0.0;
    const auto chain = chain_group(subset, parameters, chain_score);
    if (chain.size() < parameters.refinement_min_chain_anchors ||
        chain_score < parameters.min_chain_score) {
        return false;
    }
    auto emitted = make_emitted_chain(next_chain_id++, subset, chain, chain_score);
    write_refined_chain(chains, chain_anchors, emitted, refinement_type,
                        parent_chain_id, left_chain_id, right_chain_id);
    for (const auto& candidate : emitted.anchors) used.insert(candidate_key(candidate));
    ++result.refined_chain_count;
    result.refined_chain_anchor_count += emitted.anchors.size();
    return true;
}

void refine_chains(
    const std::vector<std::pair<std::string, std::vector<Candidate>>>& leftovers,
    std::vector<EmittedChain> primary_chains,
    const std::filesystem::path& refined_chain_output,
    const std::filesystem::path& refined_chain_anchor_output,
    const ChainingParameters& parameters,
    ChainingResult& result) {
    std::ofstream chains(refined_chain_output);
    if (!chains) throw std::runtime_error("cannot create refined chain output");
    std::ofstream chain_anchors(refined_chain_anchor_output);
    if (!chain_anchors) throw std::runtime_error("cannot create refined chain anchor output");
    chains << "chain_id\trefinement_type\tparent_chain_id\tleft_chain_id"
              "\tright_chain_id\tsample_a\tsample_b\tsequence_a\tsequence_b"
              "\tassign_strand\tn_candidates\tchain_score"
              "\tref_start\tref_end\tquery_start\tquery_end\n";
    chain_anchors << "chain_id\trank\trefinement_type\tsample_a\tsample_b"
                     "\tanchor_a\tanchor_b\tref_center\tquery_center"
                     "\tassign_strand\tcandidate_score\tsupport_direction\n";

    std::uint64_t next_chain_id = 0;
    std::unordered_set<std::string> refined_used;
    for (const auto& primary : primary_chains) {
        for (const auto& candidate : primary.anchors) {
            refined_used.insert(candidate_key(candidate));
        }
    }

    for (const auto& item : leftovers) {
        const auto& key = item.first;
        const auto& candidates = item.second;
        for (const auto& primary : primary_chains) {
            if (!same_chain_group(primary, key)) continue;
            std::vector<Candidate> subset;
            const auto ref_start = subtract_window(primary.ref_start,
                                                   parameters.refinement_window);
            const auto ref_end = primary.ref_end + parameters.refinement_window;
            const auto query_start = subtract_window(primary.query_start,
                                                     parameters.refinement_window);
            const auto query_end = primary.query_end + parameters.refinement_window;
            for (const auto& candidate : candidates) {
                if (candidate_in_box(candidate, ref_start, ref_end, query_start,
                                     query_end, refined_used)) {
                    subset.push_back(candidate);
                }
            }
            emit_refined_subset(std::move(subset), parameters, next_chain_id, result,
                                refined_used, chains, chain_anchors, "around",
                                std::to_string(primary.id), ".", ".");
        }
    }

    std::sort(primary_chains.begin(), primary_chains.end(),
              [](const EmittedChain& left, const EmittedChain& right) {
                  const auto left_key = left.sample_a + '\x1f' + left.sample_b + '\x1f' +
                                        left.sequence_a + '\x1f' + left.sequence_b + '\x1f' +
                                        left.strand;
                  const auto right_key = right.sample_a + '\x1f' + right.sample_b + '\x1f' +
                                         right.sequence_a + '\x1f' + right.sequence_b + '\x1f' +
                                         right.strand;
                  if (left_key != right_key) return left_key < right_key;
                  return left.ref_start < right.ref_start;
              });
    for (std::size_t index = 1; index < primary_chains.size(); ++index) {
        const auto& left = primary_chains[index - 1];
        const auto& right = primary_chains[index];
        Candidate left_group;
        left_group.sample_a = left.sample_a;
        left_group.sample_b = left.sample_b;
        left_group.sequence_a = left.sequence_a;
        left_group.sequence_b = left.sequence_b;
        left_group.strand = left.strand;
        const auto key = group_key(left_group);
        if (!same_chain_group(right, key)) continue;
        if (left.ref_end > right.ref_start) continue;

        std::uint64_t query_start = 0;
        std::uint64_t query_end = 0;
        if (left.strand == '+') {
            if (left.query_end > right.query_start) continue;
            query_start = left.query_end;
            query_end = right.query_start;
        } else {
            if (right.query_end > left.query_start) continue;
            query_start = right.query_end;
            query_end = left.query_start;
        }

        const auto found = std::find_if(
            leftovers.begin(), leftovers.end(),
            [&](const auto& item) { return item.first == key; });
        if (found == leftovers.end()) continue;
        std::vector<Candidate> subset;
        for (const auto& candidate : found->second) {
            if (candidate_in_box(candidate, left.ref_end, right.ref_start,
                                 query_start, query_end, refined_used)) {
                subset.push_back(candidate);
            }
        }
        emit_refined_subset(std::move(subset), parameters, next_chain_id, result,
                            refined_used, chains, chain_anchors, "between", ".",
                            std::to_string(left.id), std::to_string(right.id));
    }
}

void write_metadata(const std::filesystem::path& path,
                    const ChainingParameters& parameters,
                    const ChainingResult& result) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create chain metadata");
    output << "{\n"
           << "  \"chain_predecessors\": " << parameters.predecessor_count << ",\n"
           << "  \"gap_weight\": " << parameters.gap_weight << ",\n"
           << "  \"gap_unit\": " << parameters.gap_unit << ",\n"
           << "  \"min_chain_anchors\": " << parameters.min_chain_anchors << ",\n"
           << "  \"min_chain_score\": " << parameters.min_chain_score << ",\n"
           << "  \"refinement_window\": " << parameters.refinement_window << ",\n"
           << "  \"refinement_min_chain_anchors\": "
           << parameters.refinement_min_chain_anchors << ",\n"
           << "  \"candidate_count\": " << result.candidate_count << ",\n"
           << "  \"chain_count\": " << result.chain_count << ",\n"
           << "  \"chain_anchor_count\": " << result.chain_anchor_count << ",\n"
           << "  \"refined_chain_count\": " << result.refined_chain_count << ",\n"
           << "  \"refined_chain_anchor_count\": "
           << result.refined_chain_anchor_count << "\n"
           << "}\n";
}

}  // namespace

ChainingResult build_anchor_chains(
    const std::filesystem::path& correspondences,
    const std::filesystem::path& assignments,
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& chain_output,
    const std::filesystem::path& chain_anchor_output,
    const std::filesystem::path& metadata_output,
    const std::filesystem::path& refined_chain_output,
    const std::filesystem::path& refined_chain_anchor_output,
    const std::filesystem::path& refined_metadata_output,
    const ChainingParameters& parameters) {
    if (parameters.predecessor_count == 0) {
        throw std::runtime_error("chain predecessor count must be greater than zero");
    }
    if (parameters.gap_unit == 0) {
        throw std::runtime_error("gap unit must be greater than zero");
    }
    if (parameters.min_chain_anchors == 0) {
        throw std::runtime_error("minimum chain anchor count must be greater than zero");
    }
    if (parameters.refinement_min_chain_anchors == 0) {
        throw std::runtime_error("minimum refinement chain anchor count must be greater than zero");
    }
    const auto positions = load_anchor_positions(grouped_anchors);
    auto candidates = load_candidates(correspondences, positions);
    auto assignment_candidates = load_assignment_candidates(assignments, positions);
    ChainingResult result;
    result.candidate_count = candidates.size();

    std::unordered_map<std::string, std::vector<Candidate>> groups;
    for (auto& candidate : candidates) {
        groups[group_key(candidate)].push_back(std::move(candidate));
    }

    std::ofstream chains(chain_output);
    if (!chains) throw std::runtime_error("cannot create chain output");
    std::ofstream chain_anchors(chain_anchor_output);
    if (!chain_anchors) throw std::runtime_error("cannot create chain anchor output");
    chains << "chain_id\tsample_a\tsample_b\tsequence_a\tsequence_b"
              "\tassign_strand\tn_candidates\tchain_score"
              "\tref_start\tref_end\tquery_start\tquery_end\n";
    chain_anchors << "chain_id\trank\tsample_a\tsample_b\tanchor_a\tanchor_b"
                     "\tref_center\tquery_center\tassign_strand"
                     "\tcandidate_score\tsupport_direction\n";

    std::vector<std::pair<std::string, std::vector<Candidate>>> ordered_groups(
        std::make_move_iterator(groups.begin()),
        std::make_move_iterator(groups.end()));
    std::sort(ordered_groups.begin(), ordered_groups.end(),
              [](const auto& left, const auto& right) {
                  return left.first < right.first;
              });

    std::uint64_t next_chain_id = 0;
    std::vector<EmittedChain> emitted_chains;
    for (auto& item : ordered_groups) {
        auto& group = item.second;
        std::vector<Candidate> refinement_candidates;
        while (!group.empty()) {
            double chain_score = 0.0;
            const auto chain = chain_group(group, parameters, chain_score);
            if (chain.empty()) break;
            if (chain_score < parameters.min_chain_score) {
                break;
            }
            std::vector<bool> used(group.size(), false);
            for (const auto index : chain) used[index] = true;
            if (chain.size() < parameters.min_chain_anchors) {
                std::vector<Candidate> remaining;
                remaining.reserve(group.size() - chain.size());
                for (std::size_t index = 0; index < group.size(); ++index) {
                    if (used[index])
                        refinement_candidates.push_back(std::move(group[index]));
                    else
                        remaining.push_back(std::move(group[index]));
                }
                group = std::move(remaining);
                continue;
            }

            auto emitted = make_emitted_chain(next_chain_id, group, chain, chain_score);
            write_chain(chains, chain_anchors, emitted);
            emitted_chains.push_back(std::move(emitted));
            ++next_chain_id;
            ++result.chain_count;
            result.chain_anchor_count += chain.size();

            std::vector<Candidate> remaining;
            remaining.reserve(group.size() - chain.size());
            for (std::size_t index = 0; index < group.size(); ++index) {
                if (!used[index]) remaining.push_back(std::move(group[index]));
            }
            group = std::move(remaining);
        }
        group.insert(group.end(),
                     std::make_move_iterator(refinement_candidates.begin()),
                     std::make_move_iterator(refinement_candidates.end()));
    }

    std::unordered_map<std::string, std::size_t> leftover_index;
    for (std::size_t index = 0; index < ordered_groups.size(); ++index) {
        leftover_index.emplace(ordered_groups[index].first, index);
    }
    for (auto& candidate : assignment_candidates) {
        const auto key = group_key(candidate);
        const auto found = leftover_index.find(key);
        if (found == leftover_index.end()) {
            leftover_index.emplace(key, ordered_groups.size());
            ordered_groups.emplace_back(key, std::vector<Candidate>{});
            ordered_groups.back().second.push_back(std::move(candidate));
        } else {
            ordered_groups[found->second].second.push_back(std::move(candidate));
        }
    }
    std::sort(ordered_groups.begin(), ordered_groups.end(),
              [](const auto& left, const auto& right) {
                  return left.first < right.first;
              });

    refine_chains(ordered_groups, std::move(emitted_chains),
                  refined_chain_output, refined_chain_anchor_output,
                  parameters, result);
    write_metadata(metadata_output, parameters, result);
    write_metadata(refined_metadata_output, parameters, result);
    return result;
}

}  // namespace vallescope2
