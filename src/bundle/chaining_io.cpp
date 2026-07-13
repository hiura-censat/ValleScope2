#include "chaining_internal.hpp"

#include "vallescope2/context/anchor_grouping.hpp"

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace vallescope2 {
namespace chaining_detail {
namespace {

const char* gap_cost_model_name(const GapCostModel model) {
    switch (model) {
        case GapCostModel::absolute:
            return "absolute";
        case GapCostModel::relative:
            return "relative";
    }
    return "unknown";
}

double emitted_both_rate(const EmittedChain& chain) {
    if (chain.anchors.empty()) return 0.0;
    return static_cast<double>(chain.both_anchor_count) /
           static_cast<double>(chain.anchors.size());
}

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

}  // namespace

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
    const auto sequence_text_column = require_column(columns, "anchor_seq");

    std::unordered_map<std::string, AnchorPosition> positions;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= end_column) throw std::runtime_error("invalid grouped anchor row");
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
        const char strand = choose_candidate_strand(
            score_and_strand.second, position_a->second, position_b->second);
        if (strand != '+' && strand != '-') continue;
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
        candidate.strand = strand;
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
        const char reported_strand = parse_optional_strand(fields[strand_column]);
        const char strand = choose_candidate_strand(
            reported_strand, position_a->second, position_b->second);
        if (strand != '+' && strand != '-') continue;
        candidate.strand = strand;
        candidates.push_back(std::move(candidate));
    }
    return candidates;
}

void write_chain(std::ostream& chains,
                 std::ostream& chain_anchors,
                 const EmittedChain& chain) {
    chains << chain.id << '\t' << chain.sample_a << '\t'
           << chain.sample_b << '\t' << chain.sequence_a << '\t'
           << chain.sequence_b << '\t' << chain.strand << '\t'
           << chain.anchors.size() << '\t' << chain.both_anchor_count << '\t'
           << emitted_both_rate(chain) << '\t'
           << chain.raw_score << '\t' << chain.score << '\t'
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
           << chain.anchors.size() << '\t' << chain.both_anchor_count << '\t'
           << emitted_both_rate(chain) << '\t'
           << chain.raw_score << '\t' << chain.score << '\t'
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

void write_metadata(const std::filesystem::path& path,
                    const ChainingParameters& parameters,
                    const ChainingResult& result) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create chain metadata");
    output << "{\n"
           << "  \"chain_predecessors\": " << parameters.predecessor_count << ",\n"
           << "  \"max_chain_gap\": " << parameters.max_chain_gap << ",\n"
           << "  \"chain_max_gap_ratio\": "
           << parameters.chain_max_gap_ratio << ",\n"
           << "  \"gap_cost_model\": \""
           << gap_cost_model_name(parameters.gap_cost_model) << "\",\n"
           << "  \"gap_weight\": " << parameters.gap_weight << ",\n"
           << "  \"gap_unit\": " << parameters.gap_unit << ",\n"
           << "  \"min_chain_anchors\": " << parameters.min_chain_anchors << ",\n"
           << "  \"min_chain_both_anchors\": "
           << parameters.min_chain_both_anchors << ",\n"
           << "  \"min_chain_score\": " << parameters.min_chain_score << ",\n"
           << "  \"chain_trim_overlap\": " << parameters.chain_trim_overlap << ",\n"
           << "  \"reciprocal_score_weight\": "
           << parameters.reciprocal_score_weight << ",\n"
           << "  \"multimap_overlap\": " << parameters.multimap_overlap << ",\n"
           << "  \"refinement_window\": " << parameters.refinement_window << ",\n"
           << "  \"refinement_min_chain_anchors\": "
           << parameters.refinement_min_chain_anchors << ",\n"
           << "  \"candidate_count\": " << result.candidate_count << ",\n"
           << "  \"raw_chain_count\": " << result.raw_chain_count << ",\n"
           << "  \"trimmed_chain_count\": " << result.trimmed_chain_count << ",\n"
           << "  \"chain_count\": " << result.chain_count << ",\n"
           << "  \"chain_anchor_count\": " << result.chain_anchor_count << ",\n"
           << "  \"refined_chain_count\": " << result.refined_chain_count << ",\n"
           << "  \"refined_chain_anchor_count\": "
           << result.refined_chain_anchor_count << "\n"
           << "}\n";
}

}  // namespace chaining_detail
}  // namespace vallescope2
