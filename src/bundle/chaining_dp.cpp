#include "chaining_internal.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace vallescope2 {
namespace chaining_detail {

std::string group_key(const Candidate& candidate) {
    return candidate.sample_a + '\x1f' + candidate.sample_b + '\x1f' +
           candidate.sequence_a + '\x1f' + candidate.sequence_b + '\x1f' +
           candidate.strand;
}

std::string candidate_key(const Candidate& candidate) {
    return candidate.anchor_a + '\x1f' + candidate.anchor_b;
}

namespace {

bool predecessor_allowed(const Candidate& predecessor,
                         const Candidate& current,
                         const ChainingParameters& parameters) {
    if (predecessor.ref_center >= current.ref_center) return false;
    const auto dr = current.ref_center - predecessor.ref_center;
    const auto dq = current.query_center > predecessor.query_center
                        ? current.query_center - predecessor.query_center
                        : predecessor.query_center - current.query_center;
    if (dr > parameters.max_chain_gap || dq > parameters.max_chain_gap) {
        return false;
    }
    const auto smaller_gap = std::min(dr, dq);
    const auto larger_gap = std::max(dr, dq);
    if (smaller_gap == 0 ||
        static_cast<double>(larger_gap) / static_cast<double>(smaller_gap) >
            parameters.chain_max_gap_ratio) {
        return false;
    }
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

}  // namespace

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
            if (!predecessor_allowed(candidates[j], candidates[i], parameters)) continue;
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
                                const double chain_score,
                                const ChainingParameters& parameters) {
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
    emitted.raw_score = chain_score;
    emitted.anchors.reserve(chain.size());
    for (const auto index : chain) {
        emitted.query_start = std::min(emitted.query_start, group[index].query_center);
        emitted.query_end = std::max(emitted.query_end, group[index].query_center);
        if (group[index].support_direction == "both") ++emitted.both_anchor_count;
        emitted.anchors.push_back(group[index]);
    }
    const auto both_rate = emitted.anchors.empty()
        ? 0.0
        : static_cast<double>(emitted.both_anchor_count) /
              static_cast<double>(emitted.anchors.size());
    emitted.score = emitted.raw_score *
        ((1.0 - parameters.reciprocal_score_weight) +
         parameters.reciprocal_score_weight * both_rate);
    return emitted;
}

}  // namespace chaining_detail
}  // namespace vallescope2
