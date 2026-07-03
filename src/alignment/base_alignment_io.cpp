#include "base_alignment_internal.hpp"

#include "vallescope2/context/anchor_grouping.hpp"

#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <unordered_set>

namespace vallescope2 {
namespace alignment_detail {

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

double parse_f64(const std::string& value, const std::string& column) {
    try {
        std::size_t parsed = 0;
        const auto number = std::stod(value, &parsed);
        if (value.empty() || parsed != value.size()) throw std::out_of_range("bad");
        return number;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid number in " + column + ": " + value);
    }
}

struct AnchorPosition {
    std::string sequence_id;
    std::string sequence;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::uint64_t center = 0;
};

std::unordered_map<std::string, AnchorPosition> load_anchor_positions(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open grouped anchors: " + path.string());
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("empty grouped anchors file");
    }
    const auto columns = header_index(split_tab(line));
    const auto anchor_column = require_column(columns, "anchor_id");
    const auto sequence_column = require_column(columns, "sequence_id");
    const auto start_column = require_column(columns, "start");
    const auto end_column = require_column(columns, "end");
    const auto sequence_text_column = require_column(columns, "anchor_seq");

    std::unordered_map<std::string, AnchorPosition> positions;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= end_column) {
            throw std::runtime_error("invalid grouped anchor row");
        }
        AnchorPosition position;
        position.sequence_id = fields[sequence_column];
        position.sequence = fields[sequence_text_column];
        position.start = parse_u64(fields[start_column], "start");
        position.end = parse_u64(fields[end_column], "end");
        position.center = (position.start + position.end) / 2;
        positions.emplace(fields[anchor_column], std::move(position));
    }
    return positions;
}

std::vector<ChainBundle> load_chains(const std::filesystem::path& path,
                                     const std::string& source) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open chains: " + path.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty chains file");
    const auto columns = header_index(split_tab(line));
    const auto chain_column = require_column(columns, "chain_id");
    const auto sample_a_column = require_column(columns, "sample_a");
    const auto sample_b_column = require_column(columns, "sample_b");
    const auto sequence_a_column = require_column(columns, "sequence_a");
    const auto sequence_b_column = require_column(columns, "sequence_b");
    const auto strand_column = require_column(columns, "assign_strand");
    const auto score_column = require_column(columns, "chain_score");
    const auto ref_start_column = require_column(columns, "ref_start");
    const auto ref_end_column = require_column(columns, "ref_end");
    const auto query_start_column = require_column(columns, "query_start");
    const auto query_end_column = require_column(columns, "query_end");

    std::vector<ChainBundle> chains;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= query_end_column)
            throw std::runtime_error("invalid chain row");
        ChainBundle chain;
        chain.source = source;
        chain.chain_id = parse_u64(fields[chain_column], "chain_id");
        chain.sample_a = fields[sample_a_column];
        chain.sample_b = fields[sample_b_column];
        chain.sequence_a = fields[sequence_a_column];
        chain.sequence_b = fields[sequence_b_column];
        chain.strand = fields[strand_column] == "-" ? '-' : '+';
        chain.score = parse_f64(fields[score_column], "chain_score");
        chain.ref_start = parse_u64(fields[ref_start_column], "ref_start");
        chain.ref_end = parse_u64(fields[ref_end_column], "ref_end");
        chain.query_start = parse_u64(fields[query_start_column], "query_start");
        chain.query_end = parse_u64(fields[query_end_column], "query_end");
        chains.push_back(std::move(chain));
    }
    return chains;
}

void attach_chain_anchors(
    std::vector<ChainBundle>& chains,
    const std::filesystem::path& path,
    const std::unordered_map<std::string, AnchorPosition>& positions,
    AnchorStore& store) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open chain anchors: " + path.string());
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("empty chain anchors file");
    }
    const auto columns = header_index(split_tab(line));
    const auto chain_column = require_column(columns, "chain_id");
    const auto anchor_a_column = require_column(columns, "anchor_a");
    const auto anchor_b_column = require_column(columns, "anchor_b");

    std::unordered_map<std::uint64_t, ChainBundle*> chain_by_id;
    for (auto& chain : chains) chain_by_id.emplace(chain.chain_id, &chain);

    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= anchor_b_column) {
            throw std::runtime_error("invalid chain anchor row");
        }
        const auto chain_id = parse_u64(fields[chain_column], "chain_id");
        const auto chain_found = chain_by_id.find(chain_id);
        if (chain_found == chain_by_id.end()) continue;

        const auto anchor_a = fields[anchor_a_column];
        const auto anchor_b = fields[anchor_b_column];
        const auto position_a = positions.find(anchor_a);
        const auto position_b = positions.find(anchor_b);
        if (position_a == positions.end() || position_b == positions.end()) {
            throw std::runtime_error("chain anchor is absent from grouped anchors");
        }

        AnchorPair anchor;
        anchor.anchor_a = anchor_a;
        anchor.anchor_b = anchor_b;
        anchor.ref_start = position_a->second.start;
        anchor.ref_end = position_a->second.end;
        anchor.query_start = position_b->second.start;
        anchor.query_end = position_b->second.end;
        const auto index = store.anchors.size();
        store.anchors.push_back(std::move(anchor));
        chain_found->second->anchor_indices.push_back(index);
    }
}

std::vector<ChainBundle> load_chain_files(
    const std::vector<std::filesystem::path>& paths,
    const std::vector<std::filesystem::path>& anchor_paths,
    const std::filesystem::path& grouped_anchors,
    AnchorStore& store) {
    if (paths.size() != anchor_paths.size()) {
        throw std::runtime_error("chain file and chain anchor file counts differ");
    }
    const auto positions = load_anchor_positions(grouped_anchors);
    std::vector<ChainBundle> bundles;
    for (std::size_t i = 0; i < paths.size(); ++i) {
        const auto& path = paths[i];
        const std::string filename = path.filename().string();
        const std::string source =
            filename.find("refined") == std::string::npos ? "primary" : "refined";
        auto loaded = load_chains(path, source);
        attach_chain_anchors(loaded, anchor_paths[i], positions, store);
        bundles.insert(bundles.end(),
                       std::make_move_iterator(loaded.begin()),
                       std::make_move_iterator(loaded.end()));
    }
    return bundles;
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

char parse_optional_strand(const std::string& value) {
    if (value == "-" || value == "+") return value.front();
    return '.';
}

char infer_anchor_pair_strand(const AnchorPosition& reference,
                              const AnchorPosition& query) {
    if (reference.sequence.empty() || query.sequence.empty()) return '.';
    const bool forward = reference.sequence == query.sequence;
    const bool reverse = reference.sequence == reverse_complement(query.sequence);
    if (forward == reverse) return '.';
    return forward ? '+' : '-';
}

char choose_candidate_strand(const char reported,
                             const AnchorPosition& reference,
                             const AnchorPosition& query) {
    if (reported == '+' || reported == '-') return reported;
    return infer_anchor_pair_strand(reference, query);
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

ExtensionCandidate make_extension_candidate(
    const std::string& sample_a,
    const std::string& sample_b,
    const std::string& anchor_a,
    const std::string& anchor_b,
    const std::string& support_direction,
    const double score,
    const char strand,
    const AnchorPosition& position_a,
    const AnchorPosition& position_b) {
    ExtensionCandidate candidate;
    candidate.sample_a = sample_a;
    candidate.sample_b = sample_b;
    candidate.sequence_a = position_a.sequence_id;
    candidate.sequence_b = position_b.sequence_id;
    candidate.anchor_a = anchor_a;
    candidate.anchor_b = anchor_b;
    candidate.support_direction = support_direction;
    candidate.ref_start = position_a.start;
    candidate.ref_end = position_a.end;
    candidate.ref_center = position_a.center;
    candidate.query_start = position_b.start;
    candidate.query_end = position_b.end;
    candidate.query_center = position_b.center;
    candidate.score = score;
    candidate.strand = strand == '-' ? '-' : '+';
    return candidate;
}

void load_correspondence_extension_candidates(
    const std::filesystem::path& path,
    const std::unordered_map<std::string, AnchorPosition>& positions,
    std::vector<ExtensionCandidate>& candidates,
    std::unordered_set<std::string>& seen) {
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

    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= reverse_strand_column) {
            throw std::runtime_error("invalid correspondence row");
        }
        const auto& anchor_a = fields[anchor_a_column];
        const auto& anchor_b = fields[anchor_b_column];
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
        const char strand = choose_candidate_strand(
            score_and_strand.second, position_a->second, position_b->second);
        if (strand != '+' && strand != '-') continue;
        const std::string key = fields[sample_a_column] + '\x1f' +
                                fields[sample_b_column] + '\x1f' +
                                anchor_a + '\x1f' + anchor_b;
        if (!seen.insert(key).second) continue;
        candidates.push_back(make_extension_candidate(
            fields[sample_a_column], fields[sample_b_column], anchor_a, anchor_b,
            fields[support_column], score_and_strand.first,
            strand, position_a->second, position_b->second));
    }
}

void load_assignment_extension_candidates(
    const std::filesystem::path& path,
    const std::unordered_map<std::string, AnchorPosition>& positions,
    std::vector<ExtensionCandidate>& candidates,
    std::unordered_set<std::string>& seen) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open assignments: " + path.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty assignments file");
    const auto columns = header_index(split_tab(line));
    const auto target_sample_column = require_column(columns, "target_sample");
    const auto query_sample_column = require_column(columns, "query_sample");
    const auto target_anchor_column = require_column(columns, "target_anchor_id");
    const auto query_anchor_column = require_column(columns, "query_anchor_id");
    const auto status_column = require_column(columns, "status");
    const auto strand_column = require_column(columns, "assign_strand");
    const auto score_column = require_column(columns, "candidate_score");

    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= score_column) {
            throw std::runtime_error("invalid assignment row");
        }
        const auto& status = fields[status_column];
        if (status != "primary" && status != "ambiguous") continue;
        const auto& anchor_a = fields[target_anchor_column];
        const auto& anchor_b = fields[query_anchor_column];
        if (anchor_a == "." || anchor_b == ".") continue;
        const auto position_a = positions.find(anchor_a);
        const auto position_b = positions.find(anchor_b);
        if (position_a == positions.end() || position_b == positions.end()) {
            throw std::runtime_error("assignment refers to unknown anchor");
        }
        const std::string key = fields[target_sample_column] + '\x1f' +
                                fields[query_sample_column] + '\x1f' +
                                anchor_a + '\x1f' + anchor_b;
        if (!seen.insert(key).second) continue;
        const char reported_strand = parse_optional_strand(fields[strand_column]);
        const char strand = choose_candidate_strand(
            reported_strand, position_a->second, position_b->second);
        if (strand != '+' && strand != '-') continue;
        candidates.push_back(make_extension_candidate(
            fields[target_sample_column], fields[query_sample_column],
            anchor_a, anchor_b, status,
            parse_f64(fields[score_column], "candidate_score"),
            strand, position_a->second, position_b->second));
    }
}

std::vector<ExtensionCandidate> load_extension_candidates(
    const std::filesystem::path& correspondences,
    const std::filesystem::path& assignments,
    const std::filesystem::path& grouped_anchors) {
    const auto positions = load_anchor_positions(grouped_anchors);
    std::vector<ExtensionCandidate> candidates;
    std::unordered_set<std::string> seen;
    load_correspondence_extension_candidates(
        correspondences, positions, candidates, seen);
    load_assignment_extension_candidates(
        assignments, positions, candidates, seen);
    return candidates;
}

}  // namespace alignment_detail
}  // namespace vallescope2
