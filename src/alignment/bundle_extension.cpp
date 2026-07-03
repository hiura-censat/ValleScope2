#include "base_alignment_internal.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace vallescope2 {
namespace alignment_detail {
namespace {

struct ExtensionNode {
    const ExtensionCandidate* candidate = nullptr;
    std::uint64_t ref_center = 0;
    std::uint64_t query_center = 0;
    double score = 0.0;
    bool is_seed = false;
};

struct ExtensionState {
    double dp = -std::numeric_limits<double>::infinity();
    int predecessor = -1;
};

struct ExtensionAttempt {
    std::string classification = "unresolved_gap";
    std::string stop_reason = "no_path";
    std::vector<const ExtensionCandidate*> anchors;
    double score = 0.0;
    std::uint64_t candidate_count = 0;
    std::uint64_t copy_support_both_count = 0;
    std::string copy_support_side = ".";
    std::uint64_t ref_gap = 0;
    std::uint64_t query_gap = 0;
    std::uint64_t remaining_ref_gap = 0;
    std::uint64_t remaining_query_gap = 0;
};

std::string track_key(const ExtensionCandidate& candidate) {
    return candidate.sample_a + '\x1f' + candidate.sample_b + '\x1f' +
           candidate.sequence_a + '\x1f' + candidate.sequence_b + '\x1f' +
           candidate.strand;
}

std::string track_key(const ChainBundle& bundle) {
    return bundle.sample_a + '\x1f' + bundle.sample_b + '\x1f' +
           bundle.sequence_a + '\x1f' + bundle.sequence_b + '\x1f' +
           bundle.strand;
}

std::uint64_t positive_gap(const std::uint64_t left_end,
                           const std::uint64_t right_start) {
    return right_start > left_end ? right_start - left_end : 0;
}

std::uint64_t query_gap_between_bundles(const ChainBundle& left,
                                        const ChainBundle& right) {
    if (left.strand == '+') return positive_gap(left.query_end, right.query_start);
    return right.query_end < left.query_start ? left.query_start - right.query_end : 0;
}

bool ordered_adjacent(const ChainBundle& left, const ChainBundle& right) {
    if (!same_bundle_track(left, right)) return false;
    if (right.ref_start < left.ref_end) return false;
    if (left.strand == '+') return right.query_start >= left.query_end;
    return right.query_end <= left.query_start;
}

bool query_between(const ExtensionCandidate& candidate,
                   const ChainBundle& left,
                   const ChainBundle& right) {
    if (left.strand == '+') {
        return candidate.query_center > left.query_end &&
               candidate.query_center < right.query_start;
    }
    return candidate.query_center < left.query_start &&
           candidate.query_center > right.query_end;
}

bool candidate_between(const ExtensionCandidate& candidate,
                       const ChainBundle& left,
                       const ChainBundle& right,
                       const BaseAlignmentParameters& parameters) {
    if (track_key(candidate) != track_key(left)) return false;
    if (candidate.ref_center <= left.ref_end ||
        candidate.ref_center >= right.ref_start) {
        return false;
    }
    if (!query_between(candidate, left, right)) return false;
    if (candidate.ref_center - left.ref_end > parameters.max_chain_extension_bp ||
        right.ref_start - candidate.ref_center > parameters.max_chain_extension_bp) {
        return false;
    }
    const auto qgap = query_gap_between_bundles(left, right);
    if (qgap > parameters.max_chain_extension_bp) return false;
    return true;
}

bool predecessor_allowed(const ExtensionNode& predecessor,
                         const ExtensionNode& current,
                         const BaseAlignmentParameters& parameters) {
    if (predecessor.ref_center >= current.ref_center) return false;
    const auto dr = current.ref_center - predecessor.ref_center;
    const auto dq = current.query_center > predecessor.query_center
                        ? current.query_center - predecessor.query_center
                        : predecessor.query_center - current.query_center;
    if (dr > parameters.max_chain_extension_bp ||
        dq > parameters.max_chain_extension_bp) {
        return false;
    }
    if (dr == 0 || dq == 0) return false;
    const auto smaller_gap = std::min(dr, dq);
    const auto larger_gap = std::max(dr, dq);
    if (static_cast<double>(larger_gap) / static_cast<double>(smaller_gap) >
        parameters.chain_max_gap_ratio) {
        return false;
    }
    if (current.candidate && current.candidate->strand == '-') {
        return predecessor.query_center > current.query_center;
    }
    if (predecessor.candidate && predecessor.candidate->strand == '-') {
        return predecessor.query_center > current.query_center;
    }
    return predecessor.query_center < current.query_center;
}

double gap_cost(const ExtensionNode& predecessor,
                const ExtensionNode& current,
                const BaseAlignmentParameters& parameters) {
    const auto dr = current.ref_center - predecessor.ref_center;
    const auto dq = current.query_center > predecessor.query_center
                        ? current.query_center - predecessor.query_center
                        : predecessor.query_center - current.query_center;
    const auto delta = dr > dq ? dr - dq : dq - dr;
    if (parameters.gap_cost_model == GapCostModel::relative) {
        const double denominator =
            static_cast<double>(std::max<std::uint64_t>(1, std::min(dr, dq)));
        const double x = static_cast<double>(delta) / denominator;
        return parameters.gap_weight * x + std::sqrt(1.0 + x) - 1.0;
    }
    const double x = static_cast<double>(delta) /
                     static_cast<double>(parameters.gap_unit);
    return parameters.gap_weight * x + std::log2(1.0 + x);
}

std::vector<const ExtensionCandidate*> ordered_gap_candidates(
    const std::vector<ExtensionCandidate>& candidates,
    const ChainBundle& left,
    const ChainBundle& right,
    const BaseAlignmentParameters& parameters) {
    std::vector<const ExtensionCandidate*> selected;
    for (const auto& candidate : candidates) {
        if (candidate_between(candidate, left, right, parameters)) {
            selected.push_back(&candidate);
        }
    }
    std::sort(selected.begin(), selected.end(),
              [](const ExtensionCandidate* left_candidate,
                 const ExtensionCandidate* right_candidate) {
                  if (left_candidate->ref_center != right_candidate->ref_center) {
                      return left_candidate->ref_center < right_candidate->ref_center;
                  }
                  if (left_candidate->strand == '-') {
                      return left_candidate->query_center > right_candidate->query_center;
                  }
                  return left_candidate->query_center < right_candidate->query_center;
              });
    return selected;
}

std::uint64_t count_copy_support(
    const std::vector<ExtensionCandidate>& candidates,
    const ChainBundle& left,
    const ChainBundle& right,
    const std::string& side) {
    std::uint64_t count = 0;
    for (const auto& candidate : candidates) {
        if (track_key(candidate) != track_key(left)) continue;
        if (candidate.support_direction != "both") continue;
        if (side == "query") {
            if (left.strand == '+') {
                if (candidate.query_center > left.query_end &&
                    candidate.query_center < right.query_start) {
                    ++count;
                }
            } else if (candidate.query_center < left.query_start &&
                       candidate.query_center > right.query_end) {
                ++count;
            }
        } else if (side == "ref") {
            if (candidate.ref_center > left.ref_end &&
                candidate.ref_center < right.ref_start) {
                ++count;
            }
        }
    }
    return count;
}

ExtensionAttempt run_extension_dp(
    const ChainBundle& left,
    const ChainBundle& right,
    const std::vector<ExtensionCandidate>& all_candidates,
    const BaseAlignmentParameters& parameters) {
    ExtensionAttempt attempt;
    attempt.ref_gap = right.ref_start - left.ref_end;
    attempt.query_gap = query_gap_between_bundles(left, right);
    attempt.remaining_ref_gap = attempt.ref_gap;
    attempt.remaining_query_gap = attempt.query_gap;

    if (attempt.ref_gap > parameters.max_chain_extension_bp ||
        attempt.query_gap > parameters.max_chain_extension_bp) {
        attempt.stop_reason = "gap_too_large";
        return attempt;
    }

    if (attempt.ref_gap == 0 && attempt.query_gap > 0) {
        attempt.copy_support_side = "query";
        attempt.copy_support_both_count =
            count_copy_support(all_candidates, left, right, "query");
        if (attempt.copy_support_both_count >= parameters.min_copy_support_anchors) {
            attempt.classification = "query_dup_branch";
            attempt.stop_reason = "copy_support";
        } else {
            attempt.classification = "query_insertion";
            attempt.stop_reason = "single_axis_reached";
        }
        return attempt;
    }
    if (attempt.query_gap == 0 && attempt.ref_gap > 0) {
        attempt.copy_support_side = "ref";
        attempt.copy_support_both_count =
            count_copy_support(all_candidates, left, right, "ref");
        if (attempt.copy_support_both_count >= parameters.min_copy_support_anchors) {
            attempt.classification = "ref_dup_branch";
            attempt.stop_reason = "copy_support";
        } else {
            attempt.classification = "query_deletion";
            attempt.stop_reason = "single_axis_reached";
        }
        return attempt;
    }

    const auto selected =
        ordered_gap_candidates(all_candidates, left, right, parameters);
    attempt.candidate_count = selected.size();
    if (selected.empty()) {
        attempt.stop_reason = "no_candidate";
        return attempt;
    }

    std::vector<ExtensionNode> nodes;
    nodes.reserve(selected.size() + 1);
    ExtensionNode seed;
    seed.ref_center = left.ref_end;
    seed.query_center = left.strand == '+' ? left.query_end : left.query_start;
    seed.is_seed = true;
    nodes.push_back(seed);
    for (const auto* candidate : selected) {
        ExtensionNode node;
        node.candidate = candidate;
        node.ref_center = candidate->ref_center;
        node.query_center = candidate->query_center;
        node.score = candidate->score;
        nodes.push_back(node);
    }
    std::vector<ExtensionState> states(nodes.size());
    states.front().dp = 0.0;
    for (std::size_t i = 1; i < nodes.size(); ++i) {
        const std::size_t begin =
            i > parameters.chain_extension_predecessors
                ? i - parameters.chain_extension_predecessors
                : 0;
        for (std::size_t j = i; j-- > begin;) {
            if (!std::isfinite(states[j].dp)) continue;
            if (!predecessor_allowed(nodes[j], nodes[i], parameters)) continue;
            const double score =
                states[j].dp + nodes[i].score - gap_cost(nodes[j], nodes[i], parameters);
            if (score > states[i].dp) {
                states[i].dp = score;
                states[i].predecessor = static_cast<int>(j);
            }
        }
    }

    double best_score = -std::numeric_limits<double>::infinity();
    int best_index = -1;
    for (std::size_t i = 1; i < nodes.size(); ++i) {
        if (std::isfinite(states[i].dp) && states[i].dp > best_score) {
            best_score = states[i].dp;
            best_index = static_cast<int>(i);
        }
    }
    if (best_index < 0) {
        attempt.stop_reason = "no_reachable_candidate";
        return attempt;
    }

    std::vector<const ExtensionCandidate*> path;
    for (int index = best_index; index > 0;
         index = states[index].predecessor) {
        if (nodes[index].candidate) path.push_back(nodes[index].candidate);
        if (states[index].predecessor < 0) break;
    }
    std::reverse(path.begin(), path.end());
    if (path.size() < parameters.min_chain_extension_anchors ||
        best_score < parameters.min_chain_extension_score) {
        attempt.stop_reason = "below_threshold";
        attempt.score = best_score;
        attempt.anchors = std::move(path);
        return attempt;
    }

    const auto* last = path.back();
    attempt.classification = "partial_extension";
    attempt.stop_reason = "best_partial_path";
    attempt.anchors = std::move(path);
    attempt.score = best_score;
    attempt.remaining_ref_gap = positive_gap(last->ref_end, right.ref_start);
    if (left.strand == '+') {
        attempt.remaining_query_gap = positive_gap(last->query_end,
                                                   right.query_start);
    } else {
        attempt.remaining_query_gap =
            right.query_end < last->query_start
                ? last->query_start - right.query_end
                : 0;
    }
    return attempt;
}

void write_extension_report_header(std::ostream& output) {
    output << "extension_id\tleft_chain_id\tright_chain_id\tsample_a\tsample_b"
           << "\tsequence_a\tsequence_b\tstrand\tclassification\tstop_reason"
           << "\tn_candidate\tn_extension_anchor\textension_score"
           << "\tcopy_support_side\tcopy_support_both_count"
           << "\tref_gap\tquery_gap\tremaining_ref_gap\tremaining_query_gap"
           << "\tleft_ref_start\tleft_ref_end\tleft_query_start\tleft_query_end"
           << "\tright_ref_start\tright_ref_end\tright_query_start\tright_query_end\n";
}

void write_extension_anchor_report_header(std::ostream& output) {
    output << "extension_id\trank\tanchor_a\tanchor_b\tref_center\tquery_center"
           << "\tassign_strand\tcandidate_score\tsupport_direction\n";
}

void write_extension_report(std::ostream& output,
                            const std::uint64_t extension_id,
                            const ChainBundle& left,
                            const ChainBundle& right,
                            const ExtensionAttempt& attempt) {
    output << extension_id << '\t' << left.chain_id << '\t'
           << right.chain_id << '\t' << left.sample_a << '\t'
           << left.sample_b << '\t' << left.sequence_a << '\t'
           << left.sequence_b << '\t' << left.strand << '\t'
           << attempt.classification << '\t' << attempt.stop_reason << '\t'
           << attempt.candidate_count << '\t' << attempt.anchors.size() << '\t'
           << attempt.score << '\t' << attempt.copy_support_side << '\t'
           << attempt.copy_support_both_count << '\t' << attempt.ref_gap << '\t'
           << attempt.query_gap << '\t' << attempt.remaining_ref_gap << '\t'
           << attempt.remaining_query_gap << '\t' << left.ref_start << '\t'
           << left.ref_end << '\t' << left.query_start << '\t'
           << left.query_end << '\t' << right.ref_start << '\t'
           << right.ref_end << '\t' << right.query_start << '\t'
           << right.query_end << '\n';
}

void write_extension_anchor_report(
    std::ostream& output,
    const std::uint64_t extension_id,
    const std::vector<const ExtensionCandidate*>& anchors) {
    std::uint64_t rank = 0;
    for (const auto* anchor : anchors) {
        output << extension_id << '\t' << rank++ << '\t'
               << anchor->anchor_a << '\t' << anchor->anchor_b << '\t'
               << anchor->ref_center << '\t' << anchor->query_center << '\t'
               << anchor->strand << '\t' << anchor->score << '\t'
               << anchor->support_direction << '\n';
    }
}

void apply_partial_extension(ChainBundle& bundle,
                             const ExtensionAttempt& attempt,
                             AnchorStore& store) {
    for (const auto* candidate : attempt.anchors) {
        AnchorPair anchor;
        anchor.anchor_a = candidate->anchor_a;
        anchor.anchor_b = candidate->anchor_b;
        anchor.ref_start = candidate->ref_start;
        anchor.ref_end = candidate->ref_end;
        anchor.query_start = candidate->query_start;
        anchor.query_end = candidate->query_end;
        const auto index = store.anchors.size();
        store.anchors.push_back(std::move(anchor));
        bundle.anchor_indices.push_back(index);
    }
    const auto* last = attempt.anchors.back();
    bundle.ref_end = std::max(bundle.ref_end, last->ref_end);
    bundle.query_start = std::min(bundle.query_start, last->query_start);
    bundle.query_end = std::max(bundle.query_end, last->query_end);
    bundle.score += attempt.score;
    bundle.source = "extended";
}

void increment_classification_count(BaseAlignmentResult& result,
                                    const std::string& classification) {
    if (classification == "partial_extension") {
        ++result.extension_partial_count;
    } else if (classification == "query_dup_branch") {
        ++result.extension_query_dup_branch_count;
    } else if (classification == "ref_dup_branch") {
        ++result.extension_ref_dup_branch_count;
    } else if (classification == "query_insertion") {
        ++result.extension_query_insertion_count;
    } else if (classification == "query_deletion") {
        ++result.extension_query_deletion_count;
    } else {
        ++result.extension_unresolved_gap_count;
    }
}

}  // namespace

std::vector<ChainBundle> extend_adjacent_bundles(
    std::vector<ChainBundle> bundles,
    const std::vector<ExtensionCandidate>& candidates,
    AnchorStore& store,
    const BaseAlignmentParameters& parameters,
    const std::filesystem::path& extension_output,
    const std::filesystem::path& extension_anchor_output,
    BaseAlignmentResult& result) {
    if (!parameters.chain_extension || bundles.size() < 2) return bundles;

    std::ofstream extension_report(extension_output);
    if (!extension_report) {
        throw std::runtime_error("cannot create chain extension report");
    }
    std::ofstream extension_anchor_report(extension_anchor_output);
    if (!extension_anchor_report) {
        throw std::runtime_error("cannot create chain extension anchor report");
    }
    write_extension_report_header(extension_report);
    write_extension_anchor_report_header(extension_anchor_report);

    std::sort(bundles.begin(), bundles.end(),
              [](const ChainBundle& left, const ChainBundle& right) {
                  return std::tie(left.sample_a, left.sample_b, left.sequence_a,
                                  left.sequence_b, left.strand, left.ref_start,
                                  left.ref_end, left.query_start, left.query_end,
                                  left.source, left.chain_id) <
                         std::tie(right.sample_a, right.sample_b, right.sequence_a,
                                  right.sequence_b, right.strand, right.ref_start,
                                  right.ref_end, right.query_start,
                                  right.query_end, right.source, right.chain_id);
              });

    std::vector<ChainBundle> extended;
    std::uint64_t extension_id = 0;
    for (const auto& bundle : bundles) {
        if (extended.empty() || !ordered_adjacent(extended.back(), bundle)) {
            extended.push_back(bundle);
            continue;
        }

        auto attempt = run_extension_dp(
            extended.back(), bundle, candidates, parameters);
        write_extension_report(extension_report, extension_id,
                               extended.back(), bundle, attempt);
        write_extension_anchor_report(extension_anchor_report, extension_id,
                                      attempt.anchors);
        ++extension_id;
        ++result.extension_count;
        increment_classification_count(result, attempt.classification);

        if (attempt.classification == "partial_extension") {
            apply_partial_extension(extended.back(), attempt, store);
        }
        extended.push_back(bundle);
    }
    return extended;
}

}  // namespace alignment_detail
}  // namespace vallescope2
