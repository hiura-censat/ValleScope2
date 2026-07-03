#include "base_alignment_internal.hpp"

#include "vallescope2/context/anchor_grouping.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iterator>
#include <tuple>
#include <utility>
#include <vector>

namespace vallescope2 {
namespace alignment_detail {

std::uint32_t levenshtein_distance(const std::string& a,
                                   const std::string& b) {
    std::vector<std::uint32_t> previous(b.size() + 1);
    std::vector<std::uint32_t> current(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); ++j) {
        previous[j] = static_cast<std::uint32_t>(j);
    }
    for (std::size_t i = 1; i <= a.size(); ++i) {
        current[0] = static_cast<std::uint32_t>(i);
        for (std::size_t j = 1; j <= b.size(); ++j) {
            const auto substitution =
                previous[j - 1] + (a[i - 1] == b[j - 1] ? 0U : 1U);
            const auto insertion = current[j - 1] + 1U;
            const auto deletion = previous[j] + 1U;
            current[j] = std::min({substitution, insertion, deletion});
        }
        previous.swap(current);
    }
    return previous[b.size()];
}

double sequence_identity(const std::string& ref,
                         const std::string& query) {
    const auto denominator = std::max(ref.size(), query.size());
    if (denominator == 0) return 1.0;
    const auto distance = levenshtein_distance(ref, query);
    return 1.0 - static_cast<double>(distance) /
                     static_cast<double>(denominator);
}

struct PatchExtensionStep {
    std::uint64_t ref_consumed = 0;
    std::uint64_t query_consumed = 0;
    double identity = 0.0;
    bool success = false;
};

std::uint64_t positive_gap(const std::uint64_t left_end,
                           const std::uint64_t right_start) {
    return right_start > left_end ? right_start - left_end : 0;
}

bool same_bundle_track(const ChainBundle& left,
                       const ChainBundle& right) {
    return left.sample_a == right.sample_a &&
           left.sample_b == right.sample_b &&
           left.sequence_a == right.sequence_a &&
           left.sequence_b == right.sequence_b &&
           left.strand == right.strand;
}

std::uint64_t interval_span(const std::uint64_t start,
                            const std::uint64_t end) {
    return end > start ? end - start : 0;
}

std::uint64_t interval_overlap(const std::uint64_t a_start,
                               const std::uint64_t a_end,
                               const std::uint64_t b_start,
                               const std::uint64_t b_end) {
    const auto start = std::max(a_start, b_start);
    const auto end = std::min(a_end, b_end);
    return end > start ? end - start : 0;
}

double overlap_fraction_of_smaller(const std::uint64_t a_start,
                                   const std::uint64_t a_end,
                                   const std::uint64_t b_start,
                                   const std::uint64_t b_end) {
    const auto smaller = std::min(interval_span(a_start, a_end),
                                  interval_span(b_start, b_end));
    if (smaller == 0) return 0.0;
    return static_cast<double>(interval_overlap(a_start, a_end, b_start, b_end)) /
           static_cast<double>(smaller);
}

bool needs_bundle_trim(const ChainBundle& bundle,
                       const ChainBundle& blocker,
                       const double trim_overlap) {
    if (!same_bundle_track(bundle, blocker)) return false;
    const auto ref_overlap = overlap_fraction_of_smaller(
        bundle.ref_start, bundle.ref_end, blocker.ref_start, blocker.ref_end);
    const auto query_overlap = overlap_fraction_of_smaller(
        bundle.query_start, bundle.query_end,
        blocker.query_start, blocker.query_end);
    return ref_overlap >= trim_overlap && query_overlap >= trim_overlap;
}

void add_trimmed_fragment(std::vector<ChainBundle>& fragments,
                          const ChainBundle& bundle,
                          const std::uint64_t ref_start,
                          const std::uint64_t ref_end,
                          const std::uint64_t query_start,
                          const std::uint64_t query_end,
                          const std::uint32_t minimum_span) {
    if (interval_span(ref_start, ref_end) < minimum_span ||
        interval_span(query_start, query_end) < minimum_span) {
        return;
    }
    ChainBundle fragment = bundle;
    fragment.ref_start = ref_start;
    fragment.ref_end = ref_end;
    fragment.query_start = query_start;
    fragment.query_end = query_end;
    ++fragment.trim_count;
    fragments.push_back(std::move(fragment));
}

std::vector<ChainBundle> trim_bundle_against_blocker(
    const ChainBundle& bundle,
    const ChainBundle& blocker,
    const std::uint32_t minimum_span) {
    std::vector<ChainBundle> fragments;
    if (bundle.strand == '+') {
        add_trimmed_fragment(
            fragments, bundle,
            bundle.ref_start, std::min(bundle.ref_end, blocker.ref_start),
            bundle.query_start, std::min(bundle.query_end, blocker.query_start),
            minimum_span);
        add_trimmed_fragment(
            fragments, bundle,
            std::max(bundle.ref_start, blocker.ref_end), bundle.ref_end,
            std::max(bundle.query_start, blocker.query_end), bundle.query_end,
            minimum_span);
    } else {
        add_trimmed_fragment(
            fragments, bundle,
            bundle.ref_start, std::min(bundle.ref_end, blocker.ref_start),
            std::max(bundle.query_start, blocker.query_end), bundle.query_end,
            minimum_span);
        add_trimmed_fragment(
            fragments, bundle,
            std::max(bundle.ref_start, blocker.ref_end), bundle.ref_end,
            bundle.query_start, std::min(bundle.query_end, blocker.query_start),
            minimum_span);
    }
    return fragments;
}

std::uint64_t bundle_span_sum(const ChainBundle& bundle) {
    return interval_span(bundle.ref_start, bundle.ref_end) +
           interval_span(bundle.query_start, bundle.query_end);
}

std::vector<ChainBundle> final_trim_bundles(
    std::vector<ChainBundle> bundles,
    const BaseAlignmentParameters& parameters,
    std::uint64_t& trim_count) {
    trim_count = 0;
    if (parameters.bundle_trim_overlap <= 0.0 || bundles.size() < 2) {
        std::sort(bundles.begin(), bundles.end(),
                  [](const ChainBundle& left, const ChainBundle& right) {
                      return std::tie(left.sample_a, left.sample_b,
                                      left.sequence_a, left.sequence_b,
                                      left.strand, left.ref_start, left.ref_end,
                                      left.query_start, left.query_end,
                                      left.source, left.chain_id) <
                             std::tie(right.sample_a, right.sample_b,
                                      right.sequence_a, right.sequence_b,
                                      right.strand, right.ref_start,
                                      right.ref_end, right.query_start,
                                      right.query_end, right.source,
                                      right.chain_id);
                  });
        return bundles;
    }

    std::sort(bundles.begin(), bundles.end(),
              [](const ChainBundle& left, const ChainBundle& right) {
                  if (left.score != right.score) return left.score > right.score;
                  if (bundle_span_sum(left) != bundle_span_sum(right)) {
                      return bundle_span_sum(left) > bundle_span_sum(right);
                  }
                  return std::tie(left.sample_a, left.sample_b, left.sequence_a,
                                  left.sequence_b, left.strand, left.ref_start,
                                  left.query_start, left.source, left.chain_id) <
                         std::tie(right.sample_a, right.sample_b, right.sequence_a,
                                  right.sequence_b, right.strand, right.ref_start,
                                  right.query_start, right.source, right.chain_id);
              });

    std::vector<ChainBundle> accepted;
    for (const auto& bundle : bundles) {
        std::vector<ChainBundle> fragments{bundle};
        bool trimmed = false;
        for (const auto& blocker : accepted) {
            std::vector<ChainBundle> next_fragments;
            for (const auto& fragment : fragments) {
                if (!needs_bundle_trim(fragment, blocker,
                                       parameters.bundle_trim_overlap)) {
                    next_fragments.push_back(fragment);
                    continue;
                }
                trimmed = true;
                auto split = trim_bundle_against_blocker(
                    fragment, blocker, parameters.anchor_length);
                next_fragments.insert(next_fragments.end(),
                                      std::make_move_iterator(split.begin()),
                                      std::make_move_iterator(split.end()));
            }
            fragments = std::move(next_fragments);
            if (fragments.empty()) break;
        }
        if (trimmed) ++trim_count;
        accepted.insert(accepted.end(),
                        std::make_move_iterator(fragments.begin()),
                        std::make_move_iterator(fragments.end()));
    }

    std::sort(accepted.begin(), accepted.end(),
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
    return accepted;
}

std::uint64_t query_gap_between(const ChainBundle& left,
                                const ChainBundle& right) {
    if (left.strand == '+') return positive_gap(left.query_end, right.query_start);
    return right.query_end < left.query_start ? left.query_start - right.query_end : 0;
}

std::string fetch_patch_query(faidx_t* index,
                              const std::string& sequence_id,
                              const Interval& interval,
                              const char strand) {
    auto sequence = fetch_interval(index, sequence_id, interval);
    if (strand == '-') sequence = reverse_complement(sequence);
    return sequence;
}

std::string reversed(std::string sequence) {
    std::reverse(sequence.begin(), sequence.end());
    return sequence;
}

PatchExtensionStep trusted_wfa_prefix(const std::string& ref,
                                      const std::string& query,
                                      const BaseAlignmentParameters& parameters) {
    PatchExtensionStep best;
    if (ref.empty() || query.empty()) return best;

    const auto alignment = align_segment(ref, query, parameters);
    std::uint64_t ref_consumed = 0;
    std::uint64_t query_consumed = 0;
    std::uint64_t matches = 0;
    std::uint64_t block_length = 0;
    std::size_t ref_offset = 0;
    std::size_t query_offset = 0;

    auto consider_prefix = [&]() {
        if (ref_consumed == 0 || query_consumed == 0 ||
            block_length < parameters.anchor_length) {
            return;
        }
        const double identity =
            static_cast<double>(matches) / static_cast<double>(block_length);
        if (identity >= parameters.min_patch_identity) {
            best.ref_consumed = ref_consumed;
            best.query_consumed = query_consumed;
            best.identity = identity;
            best.success = true;
        }
    };

    std::size_t begin = 0;
    while (begin < alignment.cigar.size()) {
        std::size_t end = begin;
        while (end < alignment.cigar.size() &&
               std::isdigit(static_cast<unsigned char>(alignment.cigar[end]))) {
            ++end;
        }
        if (end == begin || end >= alignment.cigar.size()) break;
        const auto length = parse_u64(
            alignment.cigar.substr(begin, end - begin), "cigar_length");
        const char op = alignment.cigar[end];

        if (op == '=' || op == 'X' || op == 'M') {
            for (std::uint64_t i = 0; i < length; ++i) {
                if (ref_offset + i < ref.size() &&
                    query_offset + i < query.size() &&
                    ref[ref_offset + i] == query[query_offset + i]) {
                    ++matches;
                }
            }
            ref_consumed += length;
            query_consumed += length;
            block_length += length;
            ref_offset += length;
            query_offset += length;
        } else if (op == 'I') {
            query_consumed += length;
            block_length += length;
            query_offset += length;
        } else if (op == 'D') {
            ref_consumed += length;
            block_length += length;
            ref_offset += length;
        } else {
            break;
        }
        consider_prefix();
        begin = end + 1;
    }
    return best;
}

PatchExtensionStep extend_patch_window_wfa(
    faidx_t* index,
    const ChainBundle& bundle,
    const Interval& ref_interval,
    const Interval& query_interval,
    const bool reverse_direction,
    const BaseAlignmentParameters& parameters) {
    auto ref = fetch_interval(index, bundle.sequence_a, ref_interval);
    auto query = fetch_patch_query(
        index, bundle.sequence_b, query_interval, bundle.strand);
    if (reverse_direction) {
        ref = reversed(std::move(ref));
        query = reversed(std::move(query));
    }
    return trusted_wfa_prefix(ref, query, parameters);
}

bool can_patch_gap(faidx_t* index,
                   const ChainBundle& left,
                   const ChainBundle& right,
                   const BaseAlignmentParameters& parameters) {
    if (!same_bundle_track(left, right)) return false;
    if (right.ref_start < left.ref_end) return false;
    if (left.strand == '+' && right.query_start < left.query_end) return false;
    if (left.strand == '-' && right.query_end > left.query_start) return false;

    const auto ref_gap = right.ref_start - left.ref_end;
    const auto query_gap = query_gap_between(left, right);
    if (ref_gap > parameters.max_patch_gap_bp ||
        query_gap > parameters.max_patch_gap_bp) {
        return false;
    }
    if (ref_gap == 0 && query_gap == 0) return true;

    std::uint64_t left_ref_extension = 0;
    std::uint64_t left_query_extension = 0;
    std::uint64_t right_ref_extension = 0;
    std::uint64_t right_query_extension = 0;

    auto ref_gap_is_closed = [&]() {
        return left_ref_extension + right_ref_extension >= ref_gap;
    };
    auto query_gap_is_closed = [&]() {
        return left_query_extension + right_query_extension >= query_gap;
    };
    auto gap_is_closed = [&]() {
        return ref_gap_is_closed() && query_gap_is_closed();
    };

    while (!gap_is_closed()) {
        const auto remaining_ref =
            ref_gap - std::min(ref_gap, left_ref_extension + right_ref_extension);
        const auto remaining_query =
            query_gap - std::min(query_gap,
                                 left_query_extension + right_query_extension);
        if (remaining_ref == 0 || remaining_query == 0) return false;

        const auto left_ref_window =
            std::min<std::uint64_t>(parameters.patch_window_bp, remaining_ref);
        const auto left_query_window =
            std::min<std::uint64_t>(parameters.patch_window_bp, remaining_query);
        Interval left_ref{left.ref_end + left_ref_extension,
                          left.ref_end + left_ref_extension + left_ref_window};
        Interval left_query;
        if (left.strand == '+') {
            left_query = {left.query_end + left_query_extension,
                          left.query_end + left_query_extension + left_query_window};
        } else {
            left_query = {left.query_start - left_query_extension -
                              left_query_window,
                          left.query_start - left_query_extension};
        }
        const auto left_step = extend_patch_window_wfa(
            index, left, left_ref, left_query, false, parameters);
        if (!left_step.success) {
            return false;
        }
        left_ref_extension += left_step.ref_consumed;
        left_query_extension += left_step.query_consumed;
        if (gap_is_closed()) return true;
        if (ref_gap_is_closed() != query_gap_is_closed()) return false;

        const auto after_left_ref =
            ref_gap - std::min(ref_gap, left_ref_extension + right_ref_extension);
        const auto after_left_query =
            query_gap - std::min(query_gap,
                                 left_query_extension + right_query_extension);
        const auto right_ref_window =
            std::min<std::uint64_t>(parameters.patch_window_bp, after_left_ref);
        const auto right_query_window =
            std::min<std::uint64_t>(parameters.patch_window_bp, after_left_query);
        Interval right_ref{right.ref_start - right_ref_extension -
                               right_ref_window,
                           right.ref_start - right_ref_extension};
        Interval right_query;
        if (right.strand == '+') {
            right_query = {right.query_start - right_query_extension -
                               right_query_window,
                           right.query_start - right_query_extension};
        } else {
            right_query = {right.query_end + right_query_extension,
                           right.query_end + right_query_extension +
                               right_query_window};
        }
        const auto right_step = extend_patch_window_wfa(
            index, right, right_ref, right_query, true, parameters);
        if (!right_step.success) {
            return false;
        }
        right_ref_extension += right_step.ref_consumed;
        right_query_extension += right_step.query_consumed;
        if (ref_gap_is_closed() != query_gap_is_closed()) return false;
    }
    return true;
}

std::vector<ChainBundle> patch_adjacent_bundles(
    std::vector<ChainBundle> bundles,
    faidx_t* index,
    const BaseAlignmentParameters& parameters,
    std::uint64_t& patch_count) {
    patch_count = 0;
    std::sort(bundles.begin(), bundles.end(),
              [](const ChainBundle& left, const ChainBundle& right) {
                  return std::tie(left.sample_a, left.sample_b, left.sequence_a,
                                  left.sequence_b, left.strand, left.ref_start,
                                  left.ref_end, left.query_start,
                                  left.query_end, left.source, left.chain_id) <
                         std::tie(right.sample_a, right.sample_b, right.sequence_a,
                                  right.sequence_b, right.strand, right.ref_start,
                                  right.ref_end, right.query_start,
                                  right.query_end, right.source, right.chain_id);
              });

    std::vector<ChainBundle> patched;
    for (const auto& bundle : bundles) {
        if (patched.empty() ||
            !can_patch_gap(index, patched.back(), bundle, parameters)) {
            patched.push_back(bundle);
            continue;
        }
        auto& merged = patched.back();
        merged.ref_start = std::min(merged.ref_start, bundle.ref_start);
        merged.ref_end = std::max(merged.ref_end, bundle.ref_end);
        merged.query_start = std::min(merged.query_start, bundle.query_start);
        merged.query_end = std::max(merged.query_end, bundle.query_end);
        merged.score += bundle.score;
        merged.patch_count += bundle.patch_count + 1;
        merged.trim_count += bundle.trim_count;
        merged.anchor_indices.insert(merged.anchor_indices.end(),
                                     bundle.anchor_indices.begin(),
                                     bundle.anchor_indices.end());
        merged.source = "patched";
        ++patch_count;
    }
    return patched;
}

}  // namespace alignment_detail
}  // namespace vallescope2
