#include "chaining_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
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
    }
    chain.score = path_score(anchors, parameters);
    chain.anchors = std::move(anchors);
    return chain;
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
        if (is_inside_high_chain(anchor, high)) {
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
           chain.score >= parameters.min_chain_score;
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
                         return left.score > right.score;
                     });

    std::vector<EmittedChain> accepted;
    for (auto& chain : chains) {
        std::vector<EmittedChain> fragments{std::move(chain)};
        for (const auto& high : accepted) {
            std::vector<EmittedChain> next_fragments;
            for (auto& fragment : fragments) {
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

}  // namespace chaining_detail
}  // namespace vallescope2
