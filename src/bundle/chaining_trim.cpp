#include "chaining_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace vallescope2 {
namespace chaining_detail {
namespace {

std::uint64_t overlap_length(const std::uint64_t a_start,
                             const std::uint64_t a_end,
                             const std::uint64_t b_start,
                             const std::uint64_t b_end) {
    if (a_end <= b_start || b_end <= a_start) return 0;
    return std::min(a_end, b_end) - std::max(a_start, b_start);
}

std::uint64_t span_length(const std::uint64_t start,
                          const std::uint64_t end) {
    return end > start ? end - start : 1;
}

bool same_track(const EmittedChain& left,
                const EmittedChain& right) {
    return left.sample_a == right.sample_a &&
           left.sample_b == right.sample_b &&
           left.sequence_a == right.sequence_a &&
           left.sequence_b == right.sequence_b &&
           left.strand == right.strand;
}

bool is_inside_high_chain(const Candidate& candidate,
                          const EmittedChain& high) {
    return high.ref_start <= candidate.ref_center &&
           candidate.ref_center <= high.ref_end &&
           high.query_start <= candidate.query_center &&
           candidate.query_center <= high.query_end;
}

bool is_both_supported(const Candidate& candidate) {
    return candidate.support_direction == "both";
}

double both_rate(const EmittedChain& chain) {
    if (chain.anchors.empty()) return 0.0;
    return static_cast<double>(chain.both_anchor_count) /
           static_cast<double>(chain.anchors.size());
}

double adjusted_score(const double raw_score,
                      const double rate,
                      const ChainingParameters& parameters) {
    return raw_score *
        ((1.0 - parameters.reciprocal_score_weight) +
         parameters.reciprocal_score_weight * rate);
}

bool chain_has_higher_trim_priority(const EmittedChain& left,
                                    const EmittedChain& right) {
    if (left.score != right.score) return left.score > right.score;
    const auto left_rate = both_rate(left);
    const auto right_rate = both_rate(right);
    if (left_rate != right_rate) return left_rate > right_rate;
    if (left.both_anchor_count != right.both_anchor_count) {
        return left.both_anchor_count > right.both_anchor_count;
    }
    if (left.anchors.size() != right.anchors.size()) {
        return left.anchors.size() > right.anchors.size();
    }
    return std::tie(left.sample_a, left.sample_b, left.sequence_a,
                    left.sequence_b, left.strand, left.ref_start,
                    left.ref_end, left.query_start, left.query_end,
                    left.id) <
           std::tie(right.sample_a, right.sample_b, right.sequence_a,
                    right.sequence_b, right.strand, right.ref_start,
                    right.ref_end, right.query_start, right.query_end,
                    right.id);
}

double path_gap_cost(const Candidate& predecessor,
                     const Candidate& current,
                     const ChainingParameters& parameters) {
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

double path_score(const std::vector<Candidate>& anchors,
                  const ChainingParameters& parameters) {
    if (anchors.empty()) return 0.0;
    double score = anchors.front().score;
    for (std::size_t index = 1; index < anchors.size(); ++index) {
        score += anchors[index].score -
                 path_gap_cost(anchors[index - 1], anchors[index], parameters);
    }
    return score;
}

EmittedChain make_chain_from_anchors(std::vector<Candidate> anchors,
                                     const ChainingParameters& parameters) {
    EmittedChain chain;
    if (anchors.empty()) return chain;
    const auto& first = anchors.front();
    const auto& last = anchors.back();
    chain.sample_a = first.sample_a;
    chain.sample_b = first.sample_b;
    chain.sequence_a = first.sequence_a;
    chain.sequence_b = first.sequence_b;
    chain.strand = first.strand;
    chain.ref_start = first.ref_center;
    chain.ref_end = last.ref_center;
    chain.query_start = first.query_center;
    chain.query_end = first.query_center;
    for (const auto& anchor : anchors) {
        chain.query_start = std::min(chain.query_start, anchor.query_center);
        chain.query_end = std::max(chain.query_end, anchor.query_center);
        if (anchor.support_direction == "both") ++chain.both_anchor_count;
    }
    chain.raw_score = path_score(anchors, parameters);
    chain.anchors = std::move(anchors);
    chain.score = adjusted_score(
        chain.raw_score, both_rate(chain), parameters);
    return chain;
}

bool is_one_sided_multimap(const EmittedChain& low,
                           const EmittedChain& high,
                           const ChainingParameters& parameters) {
    if (!same_track(low, high) || parameters.multimap_overlap <= 0.0) {
        return false;
    }
    const auto ref_overlap = overlap_length(
        low.ref_start, low.ref_end, high.ref_start, high.ref_end);
    const auto query_overlap = overlap_length(
        low.query_start, low.query_end, high.query_start, high.query_end);
    const double ref_ratio =
        static_cast<double>(ref_overlap) /
        static_cast<double>(std::min(span_length(low.ref_start, low.ref_end),
                                     span_length(high.ref_start, high.ref_end)));
    const double query_ratio =
        static_cast<double>(query_overlap) /
        static_cast<double>(std::min(span_length(low.query_start, low.query_end),
                                     span_length(high.query_start, high.query_end)));
    return (ref_ratio >= parameters.multimap_overlap &&
            query_ratio < parameters.chain_trim_overlap) ||
           (query_ratio >= parameters.multimap_overlap &&
            ref_ratio < parameters.chain_trim_overlap);
}

bool should_trim_against(const EmittedChain& low,
                         const EmittedChain& high,
                         const double threshold) {
    if (!same_track(low, high)) return false;
    const auto ref_overlap = overlap_length(
        low.ref_start, low.ref_end, high.ref_start, high.ref_end);
    const auto query_overlap = overlap_length(
        low.query_start, low.query_end, high.query_start, high.query_end);
    const double ref_ratio =
        static_cast<double>(ref_overlap) /
        static_cast<double>(std::min(span_length(low.ref_start, low.ref_end),
                                     span_length(high.ref_start, high.ref_end)));
    const double query_ratio =
        static_cast<double>(query_overlap) /
        static_cast<double>(std::min(span_length(low.query_start, low.query_end),
                                     span_length(high.query_start, high.query_end)));
    return ref_ratio >= threshold && query_ratio >= threshold;
}

std::vector<EmittedChain> split_after_trimming(
    const EmittedChain& low,
    const EmittedChain& high,
    const ChainingParameters& parameters) {
    std::vector<EmittedChain> fragments;
    std::vector<Candidate> current;
    for (const auto& anchor : low.anchors) {
        if (is_inside_high_chain(anchor, high) &&
            !is_both_supported(anchor)) {
            if (!current.empty()) {
                fragments.push_back(
                    make_chain_from_anchors(std::move(current), parameters));
                current.clear();
            }
        } else {
            current.push_back(anchor);
        }
    }
    if (!current.empty()) {
        fragments.push_back(make_chain_from_anchors(std::move(current), parameters));
    }
    return fragments;
}

bool chain_passes_filters(const EmittedChain& chain,
                          const ChainingParameters& parameters) {
    return chain.anchors.size() >= parameters.min_chain_anchors &&
           chain.both_anchor_count >= parameters.min_chain_both_anchors &&
           chain.score >= parameters.min_chain_score;
}

struct CoordinateInterval {
    std::int64_t start = 0;
    std::int64_t end = 0;
};

struct ContextConflict {
    bool found = false;
    std::uint64_t left_id = 0;
    std::uint64_t right_id = 0;
    std::string bracket_axis;
    std::uint64_t violation = 0;
    std::uint64_t left_gap = 0;
    std::uint64_t right_gap = 0;
    double normalized_violation = 0.0;
    bool boundary_supported = false;
};

struct ExcursionConflict {
    bool found = false;
    std::size_t candidate_index = 0;
    std::uint64_t left_id = 0;
    std::uint64_t right_id = 0;
    std::uint64_t left_jump = 0;
    std::uint64_t right_jump = 0;
    std::uint64_t net_shift = 0;
    std::uint64_t excursion = 0;
    double normalized_excursion = 0.0;
    double return_ratio = 0.0;
    double local_gain = 0.0;
};

constexpr std::uint64_t kMinExcursionBp = 2000;
constexpr double kMinNormalizedExcursion = 0.5;
constexpr double kMaxExcursionReturnRatio = 0.25;
constexpr double kSkeletonGapScale = 1000.0;
constexpr double kSkeletonGapLambda = 1600.0;

CoordinateInterval ref_interval(const EmittedChain& chain) {
    return {static_cast<std::int64_t>(chain.ref_start),
            static_cast<std::int64_t>(chain.ref_end)};
}

CoordinateInterval oriented_query_interval(const EmittedChain& chain) {
    if (chain.strand == '+') {
        return {static_cast<std::int64_t>(chain.query_start),
                static_cast<std::int64_t>(chain.query_end)};
    }
    return {-static_cast<std::int64_t>(chain.query_end),
            -static_cast<std::int64_t>(chain.query_start)};
}

bool chain_precedes(const EmittedChain& left, const EmittedChain& right) {
    const auto left_query = oriented_query_interval(left);
    const auto right_query = oriented_query_interval(right);
    return left.ref_end <= right.ref_start &&
           left_query.end <= right_query.start;
}

double skeleton_weight(const EmittedChain& chain) {
    const auto ref_span = span_length(chain.ref_start, chain.ref_end);
    const auto query_span = span_length(chain.query_start, chain.query_end);
    return std::max(0.0, chain.score) +
           static_cast<double>(std::min(ref_span, query_span));
}

std::int64_t chain_diagonal(const EmittedChain& chain) {
    const auto query = oriented_query_interval(chain);
    const auto ref_mid = static_cast<std::int64_t>(chain.ref_start) +
                         static_cast<std::int64_t>(
                             (chain.ref_end - chain.ref_start) / 2);
    const auto query_mid = query.start + (query.end - query.start) / 2;
    return query_mid - ref_mid;
}

std::uint64_t absolute_difference(const std::int64_t left,
                                  const std::int64_t right) {
    return left >= right ? static_cast<std::uint64_t>(left - right)
                         : static_cast<std::uint64_t>(right - left);
}

std::uint64_t transition_gap_difference(const EmittedChain& left,
                                        const EmittedChain& right) {
    const auto left_query = oriented_query_interval(left);
    const auto right_query = oriented_query_interval(right);
    const auto ref_gap = static_cast<std::int64_t>(right.ref_start) -
                         static_cast<std::int64_t>(left.ref_end);
    const auto query_gap = right_query.start - left_query.end;
    return absolute_difference(ref_gap, query_gap);
}

double skeleton_transition_penalty(const EmittedChain& left,
                                   const EmittedChain& right) {
    const auto difference = transition_gap_difference(left, right);
    return kSkeletonGapLambda *
           std::log1p(static_cast<double>(difference) / kSkeletonGapScale);
}

std::vector<bool> maximum_weight_skeleton(
    const std::vector<EmittedChain>& chains,
    const std::vector<std::size_t>& group,
    const std::vector<bool>& excluded) {
    std::vector<std::size_t> ordered;
    ordered.reserve(group.size());
    for (const auto index : group) {
        if (!excluded[index]) ordered.push_back(index);
    }
    std::sort(ordered.begin(), ordered.end(), [&](const auto left, const auto right) {
        if (chains[left].ref_start != chains[right].ref_start) {
            return chains[left].ref_start < chains[right].ref_start;
        }
        return oriented_query_interval(chains[left]).start <
               oriented_query_interval(chains[right]).start;
    });

    std::vector<double> score(ordered.size(), 0.0);
    std::vector<int> predecessor(ordered.size(), -1);
    std::vector<bool> skeleton(chains.size(), false);
    if (ordered.empty()) return skeleton;
    std::size_t best = 0;
    for (std::size_t index = 0; index < ordered.size(); ++index) {
        score[index] = skeleton_weight(chains[ordered[index]]);
        for (std::size_t previous = 0; previous < index; ++previous) {
            if (!chain_precedes(chains[ordered[previous]], chains[ordered[index]])) {
                continue;
            }
            const double candidate = score[previous] +
                                     skeleton_weight(chains[ordered[index]]);
            if (candidate > score[index]) {
                score[index] = candidate;
                predecessor[index] = static_cast<int>(previous);
            }
        }
        if (score[index] > score[best]) best = index;
    }

    for (int index = static_cast<int>(best); index >= 0;
         index = predecessor[static_cast<std::size_t>(index)]) {
        skeleton[ordered[static_cast<std::size_t>(index)]] = true;
    }
    return skeleton;
}

ExcursionConflict find_excursion_conflict(
    const std::vector<EmittedChain>& chains,
    const std::vector<std::size_t>& group,
    const std::vector<bool>& skeleton) {
    std::vector<std::size_t> path;
    for (const auto index : group) {
        if (skeleton[index]) path.push_back(index);
    }
    std::sort(path.begin(), path.end(), [&](const auto left, const auto right) {
        if (chains[left].ref_start != chains[right].ref_start) {
            return chains[left].ref_start < chains[right].ref_start;
        }
        return oriented_query_interval(chains[left]).start <
               oriented_query_interval(chains[right]).start;
    });

    ExcursionConflict best;
    for (std::size_t index = 1; index + 1 < path.size(); ++index) {
        const auto left_index = path[index - 1];
        const auto candidate_index = path[index];
        const auto right_index = path[index + 1];
        const auto left_diagonal = chain_diagonal(chains[left_index]);
        const auto candidate_diagonal = chain_diagonal(chains[candidate_index]);
        const auto right_diagonal = chain_diagonal(chains[right_index]);
        const auto left_delta = candidate_diagonal - left_diagonal;
        const auto right_delta = right_diagonal - candidate_diagonal;
        if (!((left_delta < 0 && right_delta > 0) ||
              (left_delta > 0 && right_delta < 0))) {
            continue;
        }

        const auto left_jump = absolute_difference(
            candidate_diagonal, left_diagonal);
        const auto right_jump = absolute_difference(
            right_diagonal, candidate_diagonal);
        const auto net_shift = absolute_difference(right_diagonal, left_diagonal);
        const auto smaller_jump = std::min(left_jump, right_jump);
        if (smaller_jump <= net_shift) continue;
        const auto excursion = smaller_jump - net_shift;
        const auto candidate_span = std::max<std::uint64_t>(
            1, std::min(span_length(chains[candidate_index].ref_start,
                                    chains[candidate_index].ref_end),
                        span_length(chains[candidate_index].query_start,
                                    chains[candidate_index].query_end)));
        const auto normalized_excursion =
            static_cast<double>(excursion) /
            static_cast<double>(candidate_span);
        const auto return_ratio = static_cast<double>(net_shift) /
                                  static_cast<double>(smaller_jump);
        if (excursion < kMinExcursionBp ||
            normalized_excursion < kMinNormalizedExcursion ||
            return_ratio > kMaxExcursionReturnRatio) {
            continue;
        }

        const auto local_gain = skeleton_weight(chains[candidate_index]) -
            skeleton_transition_penalty(chains[left_index],
                                        chains[candidate_index]) -
            skeleton_transition_penalty(chains[candidate_index],
                                        chains[right_index]) +
            skeleton_transition_penalty(chains[left_index],
                                        chains[right_index]);
        if (local_gain >= 0.0 || (best.found && local_gain >= best.local_gain)) {
            continue;
        }
        best = {true, candidate_index, chains[left_index].id,
                chains[right_index].id, left_jump, right_jump, net_shift,
                excursion, normalized_excursion, return_ratio, local_gain};
    }
    return best;
}

ContextConflict bracket_conflict(
    const std::vector<EmittedChain>& chains,
    const std::vector<std::size_t>& group,
    const std::vector<bool>& skeleton,
    const std::size_t candidate_index,
    const bool bracket_query) {
    const auto& candidate = chains[candidate_index];
    const auto candidate_ref = ref_interval(candidate);
    const auto candidate_query = oriented_query_interval(candidate);
    const auto candidate_axis = bracket_query ? candidate_query : candidate_ref;
    const auto candidate_other = bracket_query ? candidate_ref : candidate_query;

    std::size_t left_index = chains.size();
    std::size_t right_index = chains.size();
    std::int64_t left_end = std::numeric_limits<std::int64_t>::min();
    std::int64_t right_start = std::numeric_limits<std::int64_t>::max();
    for (const auto index : group) {
        if (!skeleton[index]) continue;
        const auto ref = ref_interval(chains[index]);
        const auto query = oriented_query_interval(chains[index]);
        const auto axis = bracket_query ? query : ref;
        if (axis.end <= candidate_axis.start && axis.end > left_end) {
            left_end = axis.end;
            left_index = index;
        }
        if (axis.start >= candidate_axis.end && axis.start < right_start) {
            right_start = axis.start;
            right_index = index;
        }
    }
    if (left_index == chains.size() || right_index == chains.size()) return {};

    const auto left_ref = ref_interval(chains[left_index]);
    const auto left_query = oriented_query_interval(chains[left_index]);
    const auto right_ref = ref_interval(chains[right_index]);
    const auto right_query = oriented_query_interval(chains[right_index]);
    const auto left_other = bracket_query ? left_ref : left_query;
    const auto right_other = bracket_query ? right_ref : right_query;
    if (left_other.end > right_other.start) return {};

    const auto violation = std::max<std::int64_t>(
        0, std::max(left_other.end - candidate_other.start,
                    candidate_other.end - right_other.start));
    if (violation == 0) return {};

    const auto span = std::max<std::int64_t>(
        1, std::max(candidate_ref.end - candidate_ref.start,
                    candidate_query.end - candidate_query.start));
    const auto left_gap = candidate_axis.start - left_end;
    const auto right_gap = right_start - candidate_axis.end;
    const auto boundary_slack = std::max<std::int64_t>(500, span / 10);
    ContextConflict result;
    result.found = true;
    result.left_id = chains[left_index].id;
    result.right_id = chains[right_index].id;
    result.bracket_axis = bracket_query ? "query" : "ref";
    result.violation = static_cast<std::uint64_t>(violation);
    result.left_gap = static_cast<std::uint64_t>(left_gap);
    result.right_gap = static_cast<std::uint64_t>(right_gap);
    result.normalized_violation =
        static_cast<double>(violation) / static_cast<double>(span);
    result.boundary_supported =
        left_gap <= boundary_slack && right_gap <= boundary_slack;
    return result;
}

std::string context_group_key(const EmittedChain& chain) {
    return chain.sample_a + '\x1f' + chain.sample_b + '\x1f' +
           chain.sequence_a + '\x1f' + chain.sequence_b + '\x1f' +
           chain.strand;
}

}  // namespace

std::vector<EmittedChain> trim_overlapping_chains(
    std::vector<EmittedChain> chains,
    const ChainingParameters& parameters) {
    if (parameters.chain_trim_overlap <= 0.0) {
        for (std::uint64_t index = 0; index < chains.size(); ++index) {
            chains[index].id = index;
        }
        return chains;
    }

    std::stable_sort(chains.begin(), chains.end(),
                     [](const EmittedChain& left, const EmittedChain& right) {
                         return chain_has_higher_trim_priority(left, right);
                     });

    std::vector<EmittedChain> accepted;
    for (auto& chain : chains) {
        std::vector<EmittedChain> fragments{std::move(chain)};
        for (const auto& high : accepted) {
            std::vector<EmittedChain> next_fragments;
            for (auto& fragment : fragments) {
                if (is_one_sided_multimap(fragment, high, parameters)) {
                    continue;
                }
                if (should_trim_against(fragment, high,
                                        parameters.chain_trim_overlap)) {
                    auto split = split_after_trimming(fragment, high, parameters);
                    next_fragments.insert(next_fragments.end(),
                                          std::make_move_iterator(split.begin()),
                                          std::make_move_iterator(split.end()));
                } else {
                    next_fragments.push_back(std::move(fragment));
                }
            }
            fragments = std::move(next_fragments);
            if (fragments.empty()) break;
        }
        for (auto& fragment : fragments) {
            if (chain_passes_filters(fragment, parameters)) {
                accepted.push_back(std::move(fragment));
            }
        }
    }

    std::sort(accepted.begin(), accepted.end(),
              [](const EmittedChain& left, const EmittedChain& right) {
                  const auto left_key = left.sample_a + '\x1f' + left.sample_b +
                                        '\x1f' + left.sequence_a + '\x1f' +
                                        left.sequence_b + '\x1f' + left.strand;
                  const auto right_key = right.sample_a + '\x1f' + right.sample_b +
                                         '\x1f' + right.sequence_a + '\x1f' +
                                         right.sequence_b + '\x1f' + right.strand;
                  if (left_key != right_key) return left_key < right_key;
                  if (left.ref_start != right.ref_start)
                      return left.ref_start < right.ref_start;
                  return left.query_start < right.query_start;
              });
    for (std::uint64_t index = 0; index < accepted.size(); ++index) {
        accepted[index].id = index;
    }
    return accepted;
}

std::vector<EmittedChain> filter_context_conflicting_chains(
    std::vector<EmittedChain> chains,
    const std::filesystem::path& context_output,
    std::uint64_t& filtered_count) {
    constexpr double min_normalized_violation = 0.5;
    filtered_count = 0;
    std::ofstream output(context_output);
    if (!output) throw std::runtime_error("cannot create chain context output");
    output << "chain_id\tclassification\tleft_chain_id\tright_chain_id"
              "\tbracket_axis\tviolation_bp\tnormalized_violation"
              "\tleft_boundary_gap\tright_boundary_gap"
              "\tboundary_supported\tleft_jump_bp\tright_jump_bp"
              "\tnet_shift_bp\texcursion_bp\tnormalized_excursion"
              "\treturn_ratio\tlocal_gain\n";

    std::unordered_map<std::string, std::vector<std::size_t>> groups;
    for (std::size_t index = 0; index < chains.size(); ++index) {
        groups[context_group_key(chains[index])].push_back(index);
    }

    std::vector<bool> keep(chains.size(), true);
    for (const auto& item : groups) {
        std::vector<bool> excursion_excluded(chains.size(), false);
        std::vector<ExcursionConflict> excursion_conflicts(chains.size());
        while (true) {
            const auto provisional_skeleton = maximum_weight_skeleton(
                chains, item.second, excursion_excluded);
            const auto excursion = find_excursion_conflict(
                chains, item.second, provisional_skeleton);
            if (!excursion.found) break;
            excursion_excluded[excursion.candidate_index] = true;
            excursion_conflicts[excursion.candidate_index] = excursion;
        }
        const auto skeleton = maximum_weight_skeleton(
            chains, item.second, excursion_excluded);
        for (const auto index : item.second) {
            if (excursion_excluded[index]) {
                const auto& excursion = excursion_conflicts[index];
                output << chains[index].id << "\toff_skeleton_excursion\t"
                       << excursion.left_id << '\t' << excursion.right_id
                       << "\tdiagonal\t0\t0\t0\t0\tfalse\t"
                       << excursion.left_jump << '\t'
                       << excursion.right_jump << '\t'
                       << excursion.net_shift << '\t'
                       << excursion.excursion << '\t'
                       << excursion.normalized_excursion << '\t'
                       << excursion.return_ratio << '\t'
                       << excursion.local_gain << '\n';
                keep[index] = false;
                ++filtered_count;
                continue;
            }
            if (skeleton[index]) {
                output << chains[index].id << "\tsyntenic_skeleton\t.\t.\t."
                          "\t0\t0\t0\t0\tfalse\t0\t0\t0\t0\t0\t0\t0\n";
                continue;
            }
            auto conflict = bracket_conflict(
                chains, item.second, skeleton, index, true);
            const auto ref_conflict = bracket_conflict(
                chains, item.second, skeleton, index, false);
            if (!conflict.found ||
                (ref_conflict.found &&
                 ref_conflict.normalized_violation >
                     conflict.normalized_violation)) {
                conflict = ref_conflict;
            }

            std::string classification = "off_skeleton_retained";
            if (conflict.found && conflict.boundary_supported) {
                classification = "off_skeleton_structural_candidate";
            } else if (conflict.found &&
                       conflict.normalized_violation >=
                           min_normalized_violation) {
                classification = "orphan_gap_cross_copy_rejected";
                keep[index] = false;
                ++filtered_count;
            }
            output << chains[index].id << '\t' << classification << '\t';
            if (!conflict.found) {
                output << ".\t.\t.\t0\t0\t0\t0\tfalse"
                          "\t0\t0\t0\t0\t0\t0\t0\n";
            } else {
                output << conflict.left_id << '\t' << conflict.right_id << '\t'
                       << conflict.bracket_axis << '\t' << conflict.violation << '\t'
                       << conflict.normalized_violation << '\t'
                       << conflict.left_gap << '\t' << conflict.right_gap << '\t'
                       << (conflict.boundary_supported ? "true" : "false")
                       << "\t0\t0\t0\t0\t0\t0\t0\n";
            }
        }
    }

    std::vector<EmittedChain> retained;
    retained.reserve(chains.size() - filtered_count);
    for (std::size_t index = 0; index < chains.size(); ++index) {
        if (keep[index]) retained.push_back(std::move(chains[index]));
    }
    return retained;
}

}  // namespace chaining_detail
}  // namespace vallescope2
