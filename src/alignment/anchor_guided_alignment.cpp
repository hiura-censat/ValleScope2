#include "base_alignment_internal.hpp"

#include "vallescope2/context/anchor_grouping.hpp"

#include <algorithm>
#include <optional>
#include <tuple>
#include <vector>

namespace vallescope2 {
namespace alignment_detail {

struct GuidedAnchor {
    std::uint64_t ref_start = 0;
    std::uint64_t ref_end = 0;
    std::uint64_t query_oriented_start = 0;
    std::uint64_t query_oriented_end = 0;
};

Interval oriented_to_original_query_interval(const Interval& oriented,
                                             const std::uint64_t query_length,
                                             const char strand) {
    if (strand == '+') return oriented;
    return {query_length - oriented.end, query_length - oriented.start};
}

std::string fetch_oriented_query_interval(faidx_t* index,
                                          const std::string& sequence_id,
                                          const Interval& oriented,
                                          const std::uint64_t query_length,
                                          const char strand) {
    auto sequence = fetch_interval(
        index, sequence_id,
        oriented_to_original_query_interval(oriented, query_length, strand));
    if (strand == '-') sequence = reverse_complement(sequence);
    return sequence;
}

StrandSummary summarize_anchor_strand(
    faidx_t* index,
    const ChainBundle& bundle,
    const AnchorStore& store,
    const Interval& ref_interval,
    const Interval& query_interval) {
    StrandSummary summary;
    for (const auto anchor_index : bundle.anchor_indices) {
        if (anchor_index >= store.anchors.size()) continue;
        const auto& anchor = store.anchors[anchor_index];
        if (anchor.ref_start < ref_interval.start ||
            anchor.ref_end > ref_interval.end ||
            anchor.query_start < query_interval.start ||
            anchor.query_end > query_interval.end) {
            continue;
        }
        const auto ref_anchor = fetch_interval(
            index, bundle.sequence_a, {anchor.ref_start, anchor.ref_end});
        const auto query_anchor = fetch_interval(
            index, bundle.sequence_b, {anchor.query_start, anchor.query_end});
        if (ref_anchor.empty() || query_anchor.empty() ||
            ref_anchor.size() != query_anchor.size()) {
            ++summary.incompatible_count;
        } else if (ref_anchor == query_anchor) {
            ++summary.forward_count;
        } else if (ref_anchor == reverse_complement(query_anchor)) {
            ++summary.reverse_count;
        } else {
            ++summary.incompatible_count;
        }
    }
    if (summary.reverse_count > summary.forward_count) {
        summary.effective_strand = '-';
    } else if (summary.forward_count > summary.reverse_count) {
        summary.effective_strand = '+';
    } else {
        summary.effective_strand = bundle.strand;
    }
    return summary;
}

std::vector<GuidedAnchor> collect_guided_anchors(
    const ChainBundle& bundle,
    const AnchorStore& store,
    const Interval& ref_interval,
    const Interval& query_oriented_interval,
    const std::uint64_t query_length,
    const char effective_strand) {
    std::vector<GuidedAnchor> anchors;
    for (const auto anchor_index : bundle.anchor_indices) {
        if (anchor_index >= store.anchors.size()) continue;
        const auto& anchor = store.anchors[anchor_index];
        GuidedAnchor guided;
        guided.ref_start = anchor.ref_start;
        guided.ref_end = anchor.ref_end;
        if (effective_strand == '+') {
            guided.query_oriented_start = anchor.query_start;
            guided.query_oriented_end = anchor.query_end;
        } else {
            guided.query_oriented_start = query_length - anchor.query_end;
            guided.query_oriented_end = query_length - anchor.query_start;
        }
        if (guided.ref_start < ref_interval.start ||
            guided.ref_end > ref_interval.end ||
            guided.query_oriented_start < query_oriented_interval.start ||
            guided.query_oriented_end > query_oriented_interval.end ||
            guided.ref_start >= guided.ref_end ||
            guided.query_oriented_start >= guided.query_oriented_end) {
            continue;
        }
        anchors.push_back(guided);
    }
    std::sort(anchors.begin(), anchors.end(),
              [](const GuidedAnchor& left, const GuidedAnchor& right) {
                  return std::tie(left.ref_start, left.ref_end,
                                  left.query_oriented_start,
                                  left.query_oriented_end) <
                         std::tie(right.ref_start, right.ref_end,
                                  right.query_oriented_start,
                                  right.query_oriented_end);
              });

    std::vector<GuidedAnchor> ordered;
    for (const auto& anchor : anchors) {
        if (!ordered.empty() &&
            (anchor.ref_start < ordered.back().ref_end ||
             anchor.query_oriented_start < ordered.back().query_oriented_end)) {
            continue;
        }
        ordered.push_back(anchor);
    }
    return ordered;
}

std::optional<Alignment> anchor_guided_alignment(
    faidx_t* index,
    const ChainBundle& bundle,
    const AnchorStore& store,
    const Interval& ref_interval,
    const Interval& query_interval,
    const std::uint64_t query_length,
    const char effective_strand,
    const BaseAlignmentParameters& parameters) {
    Interval query_oriented_interval = query_interval;
    if (effective_strand == '-') {
        query_oriented_interval =
            {query_length - query_interval.end, query_length - query_interval.start};
    }
    const auto anchors = collect_guided_anchors(
        bundle, store, ref_interval, query_oriented_interval, query_length,
        effective_strand);
    if (anchors.empty()) return std::nullopt;

    std::vector<GuidedAnchor> validated_anchors;
    for (const auto& anchor : anchors) {
        const auto anchor_ref = fetch_interval(
            index, bundle.sequence_a, {anchor.ref_start, anchor.ref_end});
        const auto anchor_query = fetch_oriented_query_interval(
            index, bundle.sequence_b,
            {anchor.query_oriented_start, anchor.query_oriented_end},
            query_length, effective_strand);
        if (!anchor_ref.empty() && anchor_ref == anchor_query) {
            validated_anchors.push_back(anchor);
        }
    }
    if (validated_anchors.empty()) return std::nullopt;

    Alignment result;
    std::uint64_t ref_cursor = ref_interval.start;
    std::uint64_t query_cursor = query_oriented_interval.start;
    for (const auto& anchor : validated_anchors) {
        const auto ref_segment = fetch_interval(
            index, bundle.sequence_a, {ref_cursor, anchor.ref_start});
        const auto query_segment = fetch_oriented_query_interval(
            index, bundle.sequence_b, {query_cursor, anchor.query_oriented_start},
            query_length, effective_strand);
        const auto segment_alignment =
            align_segment(ref_segment, query_segment, parameters);
        append_cigar_string(result.cigar, segment_alignment.cigar);
        result.matches += segment_alignment.matches;
        result.block_length += segment_alignment.block_length;
        result.score += segment_alignment.score;

        const auto anchor_length = anchor.ref_end - anchor.ref_start;
        append_cigar_operation(result.cigar, anchor_length, '=');
        result.matches += anchor_length;
        result.block_length += anchor_length;
        result.score += static_cast<std::int64_t>(anchor_length);
        ref_cursor = anchor.ref_end;
        query_cursor = anchor.query_oriented_end;
    }

    const auto ref_tail = fetch_interval(
        index, bundle.sequence_a, {ref_cursor, ref_interval.end});
    const auto query_tail = fetch_oriented_query_interval(
        index, bundle.sequence_b, {query_cursor, query_oriented_interval.end},
        query_length, effective_strand);
    const auto tail_alignment = align_segment(ref_tail, query_tail, parameters);
    append_cigar_string(result.cigar, tail_alignment.cigar);
    result.matches += tail_alignment.matches;
    result.block_length += tail_alignment.block_length;
    result.score += tail_alignment.score;
    return result;
}

}  // namespace alignment_detail
}  // namespace vallescope2
