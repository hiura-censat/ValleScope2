#include "base_alignment_internal.hpp"

#include "bindings/cpp/WFAligner.hpp"
#include "vallescope2/context/anchor_grouping.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
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

struct PatchAttempt {
    std::string classification = "not_candidate";
    bool merge = false;
    Interval ref_interval;
    Interval query_interval;
    std::uint64_t ref_gap = 0;
    std::uint64_t query_gap = 0;
    std::uint64_t matches = 0;
    std::uint64_t block_length = 0;
    double identity = 1.0;
    char rescue_indel_op = '.';
    std::uint64_t rescue_indel_len = 0;
    std::uint64_t rescue_left_bp = 0;
    std::uint64_t rescue_right_bp = 0;
    double rescue_left_identity = 0.0;
    double rescue_right_identity = 0.0;
    std::uint64_t rescue_extra_indel_bp = 0;
    std::uint64_t rescue_long_indel_count = 0;
    std::uint64_t rescue_ref_offset = 0;
    std::uint64_t rescue_query_offset = 0;
    bool long_indel_rescue = false;
    double left_endpoint_identity = 1.0;
    double right_endpoint_identity = 1.0;
    double max_short_error_density = 0.0;
    bool low_quality_cigar = false;
    bool zdrop_checked = false;
    bool zdrop_fired = false;
    int zdrop_status = 0;
    int zdrop_prefix_status = 0;
    int zdrop_suffix_status = 0;
    std::uint64_t zdrop_ref_span = 0;
    std::uint64_t zdrop_query_span = 0;
    std::int64_t score = 0;
    std::string cigar;
};

struct PatchCigarOperation {
    std::uint64_t length = 0;
    char op = '.';
    std::uint64_t aligned_bp = 0;
    std::uint64_t matches = 0;
    std::uint64_t ref_before = 0;
    std::uint64_t query_before = 0;
};

constexpr std::uint64_t kMinRescueIndelBp = 500;
constexpr std::uint64_t kMinRescueFlankBp = 200;
constexpr double kMinRescueFlankIdentity = 0.95;
constexpr std::uint64_t kMaxRescueExtraIndelBp = 100;
constexpr std::uint64_t kPatchQualityWindowBp = 500;
constexpr double kMinPatchEndpointIdentity = 0.85;
constexpr double kMaxPatchShortErrorDensity = 0.15;
constexpr int kPatchZDrop = 20;
constexpr int kPatchZDropSteps = 1;

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

Interval query_gap_interval(const ChainBundle& left,
                            const ChainBundle& right) {
    if (left.strand == '+') return {left.query_end, right.query_start};
    return {right.query_end, left.query_start};
}

bool overlaps_opposite_strand_bundle(
    const ChainBundle& left,
    const ChainBundle& right,
    const std::vector<ChainBundle>& bundles,
    const std::uint32_t minimum_overlap) {
    const Interval ref_gap{left.ref_end, right.ref_start};
    const Interval query_gap = query_gap_interval(left, right);
    if (ref_gap.start >= ref_gap.end || query_gap.start >= query_gap.end) {
        return false;
    }

    for (const auto& bundle : bundles) {
        if (bundle.sample_a != left.sample_a ||
            bundle.sample_b != left.sample_b ||
            bundle.sequence_a != left.sequence_a ||
            bundle.sequence_b != left.sequence_b ||
            bundle.strand == left.strand) {
            continue;
        }
        if (interval_overlap(ref_gap.start, ref_gap.end,
                             bundle.ref_start, bundle.ref_end) <
            minimum_overlap) {
            continue;
        }
        if (interval_overlap(query_gap.start, query_gap.end,
                             bundle.query_start, bundle.query_end) <
            minimum_overlap) {
            continue;
        }
        return true;
    }
    return false;
}

std::string bundle_track_key(const ChainBundle& bundle) {
    return bundle.sample_a + '\x1f' + bundle.sample_b + '\x1f' +
           bundle.sequence_a + '\x1f' + bundle.sequence_b + '\x1f' +
           bundle.strand;
}

std::string fetch_patch_query(faidx_t* index,
                              const std::string& sequence_id,
                              const Interval& interval,
                              const char strand) {
    auto sequence = fetch_interval(index, sequence_id, interval);
    if (strand == '-') sequence = reverse_complement(sequence);
    return sequence;
}

PatchAttempt summarize_patch_cigar(PatchAttempt attempt,
                                   const std::string& ref,
                                   const std::string& query) {
    std::size_t ref_offset = 0;
    std::size_t query_offset = 0;
    std::vector<PatchCigarOperation> operations;
    std::size_t begin = 0;
    while (begin < attempt.cigar.size()) {
        std::size_t end = begin;
        while (end < attempt.cigar.size() &&
               std::isdigit(static_cast<unsigned char>(attempt.cigar[end]))) {
            ++end;
        }
        if (end == begin || end >= attempt.cigar.size()) {
            throw std::runtime_error("invalid patch interval CIGAR");
        }
        const auto length = parse_u64(attempt.cigar.substr(begin, end - begin),
                                      "cigar_length");
        const char op = attempt.cigar[end];
        PatchCigarOperation operation{
            length, op, 0, 0, ref_offset, query_offset};
        if (op == '=' || op == 'X' || op == 'M') {
            operation.aligned_bp = length;
            for (std::uint64_t i = 0; i < length; ++i) {
                if (ref_offset + i < ref.size() &&
                    query_offset + i < query.size() &&
                    ref[ref_offset + i] == query[query_offset + i]) {
                    ++attempt.matches;
                    ++operation.matches;
                }
            }
            attempt.block_length += length;
            ref_offset += length;
            query_offset += length;
        } else if (op == 'I') {
            attempt.block_length += length;
            query_offset += length;
        } else if (op == 'D') {
            attempt.block_length += length;
            ref_offset += length;
        } else {
            throw std::runtime_error("unsupported patch interval CIGAR operation");
        }
        operations.push_back(operation);
        begin = end + 1;
    }
    attempt.identity =
        attempt.block_length == 0
            ? 1.0
            : static_cast<double>(attempt.matches) /
                  static_cast<double>(attempt.block_length);

    auto segment_max_error_density = [](const std::vector<std::uint8_t>& errors) {
        if (errors.empty()) return 0.0;
        const auto window = std::min<std::size_t>(kPatchQualityWindowBp,
                                                  errors.size());
        std::uint64_t error_count = 0;
        for (std::size_t i = 0; i < window; ++i) error_count += errors[i];
        double maximum = static_cast<double>(error_count) /
                         static_cast<double>(window);
        for (std::size_t i = window; i < errors.size(); ++i) {
            error_count += errors[i];
            error_count -= errors[i - window];
            maximum = std::max(maximum,
                static_cast<double>(error_count) /
                static_cast<double>(window));
        }
        return maximum;
    };
    std::vector<std::uint8_t> short_error_segment;
    std::vector<std::uint8_t> left_endpoint;
    std::vector<std::uint8_t> right_endpoint;
    for (const auto& operation : operations) {
        const bool long_indel =
            (operation.op == 'I' || operation.op == 'D') &&
            operation.length >= kMinRescueIndelBp;
        if (long_indel) {
            attempt.max_short_error_density = std::max(
                attempt.max_short_error_density,
                segment_max_error_density(short_error_segment));
            short_error_segment.clear();
            continue;
        }
        for (std::uint64_t i = 0; i < operation.length; ++i) {
            const auto error = static_cast<std::uint8_t>(
                operation.op == '=' ? 0 :
                operation.op == 'M' && i < operation.matches ? 0 : 1);
            short_error_segment.push_back(error);
            if (left_endpoint.size() < kPatchQualityWindowBp) {
                left_endpoint.push_back(error);
            }
            right_endpoint.push_back(error);
            if (right_endpoint.size() > kPatchQualityWindowBp) {
                right_endpoint.erase(right_endpoint.begin());
            }
        }
    }
    attempt.max_short_error_density = std::max(
        attempt.max_short_error_density,
        segment_max_error_density(short_error_segment));
    auto endpoint_identity = [](const std::vector<std::uint8_t>& errors) {
        if (errors.empty()) return 0.0;
        const auto error_count = std::count(errors.begin(), errors.end(), 1U);
        return 1.0 - static_cast<double>(error_count) /
                         static_cast<double>(errors.size());
    };
    attempt.left_endpoint_identity = endpoint_identity(left_endpoint);
    attempt.right_endpoint_identity = endpoint_identity(right_endpoint);

    std::size_t rescue_index = operations.size();
    for (std::size_t i = 0; i < operations.size(); ++i) {
        const auto& operation = operations[i];
        if ((operation.op == 'I' || operation.op == 'D') &&
            operation.length >= kMinRescueIndelBp) {
            ++attempt.rescue_long_indel_count;
            rescue_index = i;
        }
    }
    if (attempt.rescue_long_indel_count != 1) return attempt;

    attempt.rescue_indel_op = operations[rescue_index].op;
    attempt.rescue_indel_len = operations[rescue_index].length;
    attempt.rescue_ref_offset = operations[rescue_index].ref_before;
    attempt.rescue_query_offset = operations[rescue_index].query_before;
    std::uint64_t left_matches = 0;
    std::uint64_t right_matches = 0;
    for (std::size_t i = 0; i < operations.size(); ++i) {
        const auto& operation = operations[i];
        if (i < rescue_index) {
            attempt.rescue_left_bp += operation.aligned_bp;
            left_matches += operation.matches;
        } else if (i > rescue_index) {
            attempt.rescue_right_bp += operation.aligned_bp;
            right_matches += operation.matches;
        }
        if (i != rescue_index &&
            (operation.op == 'I' || operation.op == 'D')) {
            attempt.rescue_extra_indel_bp += operation.length;
        }
    }
    attempt.rescue_left_identity =
        attempt.rescue_left_bp == 0
            ? 0.0
            : static_cast<double>(left_matches) /
                  static_cast<double>(attempt.rescue_left_bp);
    attempt.rescue_right_identity =
        attempt.rescue_right_bp == 0
            ? 0.0
            : static_cast<double>(right_matches) /
                  static_cast<double>(attempt.rescue_right_bp);
    attempt.long_indel_rescue =
        attempt.rescue_left_bp >= kMinRescueFlankBp &&
        attempt.rescue_right_bp >= kMinRescueFlankBp &&
        attempt.rescue_left_identity >= kMinRescueFlankIdentity &&
        attempt.rescue_right_identity >= kMinRescueFlankIdentity &&
        attempt.rescue_extra_indel_bp <= kMaxRescueExtraIndelBp;
    attempt.low_quality_cigar =
        attempt.left_endpoint_identity < kMinPatchEndpointIdentity ||
        attempt.right_endpoint_identity < kMinPatchEndpointIdentity ||
        attempt.max_short_error_density > kMaxPatchShortErrorDensity;
    return attempt;
}

void validate_patch_with_zdrop(PatchAttempt& attempt,
                               faidx_t* index,
                               const ChainBundle& left,
                               const ChainBundle& right,
                               const BaseAlignmentParameters& parameters) {
    const Interval ref_interval{left.ref_start, right.ref_end};
    const Interval query_interval{
        std::min(left.query_start, right.query_start),
        std::max(left.query_end, right.query_end)};
    attempt.zdrop_ref_span = interval_span(ref_interval.start, ref_interval.end);
    attempt.zdrop_query_span = interval_span(query_interval.start, query_interval.end);
    attempt.zdrop_checked = true;
    if (attempt.zdrop_ref_span > parameters.max_bundle_align_bp ||
        attempt.zdrop_query_span > parameters.max_bundle_align_bp) {
        attempt.zdrop_status = wfa::WFAligner::StatusMaxStepsReached;
        attempt.zdrop_fired = true;
        return;
    }

    const auto ref = fetch_interval(index, left.sequence_a, ref_interval);
    const auto query = fetch_patch_query(
        index, left.sequence_b, query_interval, left.strand);
    const auto max_memory =
        static_cast<std::uint64_t>(parameters.max_wfa_memory_gb) *
        1024ULL * 1024ULL * 1024ULL;
    const auto local_ref_base = attempt.ref_interval.start - ref_interval.start;
    const auto local_query_base =
        left.strand == '+'
            ? attempt.query_interval.start - query_interval.start
            : query_interval.end - attempt.query_interval.end;
    const auto primary_ref = local_ref_base + attempt.rescue_ref_offset;
    const auto primary_query = local_query_base + attempt.rescue_query_offset;
    const auto suffix_ref = primary_ref +
        (attempt.rescue_indel_op == 'D' ? attempt.rescue_indel_len : 0);
    const auto suffix_query = primary_query +
        (attempt.rescue_indel_op == 'I' ? attempt.rescue_indel_len : 0);
    if (primary_ref > ref.size() || primary_query > query.size() ||
        suffix_ref > ref.size() || suffix_query > query.size()) {
        attempt.zdrop_status = wfa::WFAligner::StatusMaxStepsReached;
        attempt.zdrop_fired = true;
        return;
    }

    auto run_zdrop = [&](const std::string& segment_ref,
                         const std::string& segment_query) {
        if (segment_ref.empty() && segment_query.empty()) {
            return static_cast<int>(wfa::WFAligner::StatusAlgCompleted);
        }
        wfa::WFAlignerGapAffine2Pieces aligner(
            19, 39, 3, 81, 1,
            wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
        aligner.setMaxMemory(max_memory, max_memory);
        aligner.setHeuristicZDrop(kPatchZDrop, kPatchZDropSteps);
        return static_cast<int>(aligner.alignEnd2End(segment_ref, segment_query));
    };

    auto prefix_ref = ref.substr(0, primary_ref);
    auto prefix_query = query.substr(0, primary_query);
    std::reverse(prefix_ref.begin(), prefix_ref.end());
    std::reverse(prefix_query.begin(), prefix_query.end());
    attempt.zdrop_prefix_status = run_zdrop(prefix_ref, prefix_query);
    attempt.zdrop_suffix_status = run_zdrop(
        ref.substr(suffix_ref), query.substr(suffix_query));
    attempt.zdrop_status =
        attempt.zdrop_prefix_status != wfa::WFAligner::StatusAlgCompleted
            ? attempt.zdrop_prefix_status
            : attempt.zdrop_suffix_status;
    attempt.zdrop_fired =
        attempt.zdrop_prefix_status != wfa::WFAligner::StatusAlgCompleted ||
        attempt.zdrop_suffix_status != wfa::WFAligner::StatusAlgCompleted;
}

PatchAttempt align_patch_interval(faidx_t* index,
                                  const ChainBundle& left,
                                  const ChainBundle& right,
                                  const BaseAlignmentParameters& parameters) {
    PatchAttempt attempt;
    attempt.ref_gap = right.ref_start - left.ref_end;
    attempt.query_gap = query_gap_between(left, right);

    const auto ref_length = sequence_length(index, left.sequence_a);
    const auto query_length = sequence_length(index, left.sequence_b);
    const auto flank = static_cast<std::uint64_t>(parameters.patch_flank_bp);
    attempt.ref_interval = {
        left.ref_end > flank ? left.ref_end - flank : 0,
        std::min(ref_length, right.ref_start + flank)};
    if (left.strand == '+') {
        attempt.query_interval = {
            left.query_end > flank ? left.query_end - flank : 0,
            std::min(query_length, right.query_start + flank)};
    } else {
        attempt.query_interval = {
            right.query_end > flank ? right.query_end - flank : 0,
            std::min(query_length, left.query_start + flank)};
    }
    if (attempt.ref_interval.start >= attempt.ref_interval.end ||
        attempt.query_interval.start >= attempt.query_interval.end) {
        attempt.classification = "empty_patch_interval";
        return attempt;
    }

    const auto ref = fetch_interval(index, left.sequence_a, attempt.ref_interval);
    const auto query = fetch_patch_query(
        index, left.sequence_b, attempt.query_interval, left.strand);
    try {
        wfa::WFAlignerGapAffine2Pieces aligner(
            19, 39, 3, 81, 1,
            wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
        const auto max_memory =
            static_cast<std::uint64_t>(parameters.max_wfa_memory_gb) *
            1024ULL * 1024ULL * 1024ULL;
        aligner.setMaxMemory(max_memory, max_memory);
        const auto status = aligner.alignEnd2End(ref, query);
        if (status < 0) {
            attempt.classification = "patch_wfa_failed";
            return attempt;
        }
        attempt.cigar = aligner.getCIGAR(true);
        attempt.score = aligner.getAlignmentScore();
        attempt = summarize_patch_cigar(std::move(attempt), ref, query);
        if (attempt.identity >= 0.95) {
            attempt.classification = "patch_interval_wfa";
            attempt.merge = true;
        } else if (attempt.long_indel_rescue) {
            const bool extended_bundle =
                left.source == "extended" || right.source == "extended";
            if (attempt.low_quality_cigar || extended_bundle) {
                validate_patch_with_zdrop(
                    attempt, index, left, right, parameters);
            }
            if (attempt.zdrop_fired) {
                attempt.classification = "patch_zdrop_reject";
                attempt.merge = false;
            } else {
                attempt.classification = "patch_long_indel_rescue";
                attempt.merge = true;
            }
        } else {
            attempt.classification = "patch_interval_low_identity";
            attempt.merge = false;
        }
    } catch (const std::exception&) {
        attempt.classification = "patch_wfa_failed";
        attempt.merge = false;
    }
    return attempt;
}

PatchAttempt evaluate_patch_gap(faidx_t* index,
                                const ChainBundle& left,
                                const ChainBundle& right,
                                const std::vector<ChainBundle>& bundles,
                                const BaseAlignmentParameters& parameters,
                                std::ostream& patch_report,
                                const std::uint64_t patch_id) {
    PatchAttempt attempt;
    if (!same_bundle_track(left, right)) return attempt;
    if (right.ref_start < left.ref_end) return attempt;
    if (left.strand == '+' && right.query_start < left.query_end) return attempt;
    if (left.strand == '-' && right.query_end > left.query_start) return attempt;

    const auto ref_gap = right.ref_start - left.ref_end;
    const auto query_gap = query_gap_between(left, right);
    if (ref_gap > parameters.max_patch_gap_bp ||
        query_gap > parameters.max_patch_gap_bp) {
        return attempt;
    }
    if (overlaps_opposite_strand_bundle(
            left, right, bundles, parameters.anchor_length)) {
        attempt.classification = "opposite_strand_overlap";
        attempt.ref_gap = ref_gap;
        attempt.query_gap = query_gap;
        return attempt;
    }

    attempt = align_patch_interval(index, left, right, parameters);
    patch_report
        << patch_id << '\t' << attempt.classification << '\t'
        << left.chain_id << '\t' << right.chain_id << '\t'
        << left.sample_a << '\t' << left.sample_b << '\t'
        << left.sequence_a << '\t' << left.sequence_b << '\t'
        << left.strand << '\t' << ref_gap << '\t' << query_gap << '\t'
        << attempt.ref_interval.start << '\t' << attempt.ref_interval.end << '\t'
        << attempt.query_interval.start << '\t' << attempt.query_interval.end << '\t'
        << (attempt.ref_interval.end - attempt.ref_interval.start) << '\t'
        << (attempt.query_interval.end - attempt.query_interval.start) << '\t'
        << attempt.score << '\t' << attempt.identity << '\t'
        << attempt.rescue_indel_op << '\t' << attempt.rescue_indel_len << '\t'
        << attempt.rescue_left_bp << '\t' << attempt.rescue_right_bp << '\t'
        << attempt.rescue_left_identity << '\t'
        << attempt.rescue_right_identity << '\t'
        << attempt.rescue_extra_indel_bp << '\t'
        << attempt.rescue_long_indel_count << '\t'
        << attempt.rescue_ref_offset << '\t'
        << attempt.rescue_query_offset << '\t'
        << attempt.left_endpoint_identity << '\t'
        << attempt.right_endpoint_identity << '\t'
        << attempt.max_short_error_density << '\t'
        << (attempt.low_quality_cigar ? 1 : 0) << '\t'
        << (attempt.zdrop_checked ? 1 : 0) << '\t'
        << (attempt.zdrop_fired ? 1 : 0) << '\t'
        << attempt.zdrop_status << '\t'
        << attempt.zdrop_prefix_status << '\t'
        << attempt.zdrop_suffix_status << '\t'
        << attempt.zdrop_ref_span << '\t'
        << attempt.zdrop_query_span << '\t'
        << (attempt.merge ? 1 : 0) << '\t' << attempt.cigar << '\t'
        << left.ref_start << '\t' << left.ref_end << '\t'
        << left.query_start << '\t' << left.query_end << '\t'
        << right.ref_start << '\t' << right.ref_end << '\t'
        << right.query_start << '\t' << right.query_end << '\n';
    return attempt;
}

std::vector<ChainBundle> patch_adjacent_bundles(
    std::vector<ChainBundle> bundles,
    faidx_t* index,
    const std::vector<ExtensionCandidate>&,
    const BaseAlignmentParameters& parameters,
    const std::filesystem::path& patch_extension_output,
    std::uint64_t& patch_count) {
    patch_count = 0;
    std::ofstream patch_report(patch_extension_output);
    if (!patch_report) {
        throw std::runtime_error("cannot create patch extension report");
    }
    patch_report
        << "patch_id\tstatus\tleft_chain_id"
        << "\tright_chain_id\tsample_a\tsample_b\tsequence_a\tsequence_b"
        << "\tstrand\tref_gap\tquery_gap\tref_interval_start"
        << "\tref_interval_end\tquery_interval_start\tquery_interval_end"
        << "\tref_interval_len\tquery_interval_len\twfa_score"
        << "\tidentity\trescue_op\trescue_indel_len"
        << "\trescue_left_bp\trescue_right_bp"
        << "\trescue_left_identity\trescue_right_identity"
        << "\trescue_extra_indel_bp\trescue_long_indel_count"
        << "\trescue_ref_offset\trescue_query_offset"
        << "\tleft_endpoint_identity\tright_endpoint_identity"
        << "\tmax_short_error_density\tlow_quality_cigar"
        << "\tzdrop_checked\tzdrop_fired\tzdrop_status"
        << "\tzdrop_prefix_status\tzdrop_suffix_status"
        << "\tzdrop_ref_span\tzdrop_query_span"
        << "\tmerge\tcigar"
        << "\tleft_ref_start\tleft_ref_end"
        << "\tleft_query_start\tleft_query_end\tright_ref_start"
        << "\tright_ref_end\tright_query_start\tright_query_end\n";
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
    std::uint64_t patch_id = 0;
    for (const auto& bundle : bundles) {
        ChainBundle incoming = bundle;
        if (patched.empty()) {
            patched.push_back(std::move(incoming));
            continue;
        }
        const auto attempt = evaluate_patch_gap(
            index, patched.back(), incoming, bundles, parameters,
            patch_report, patch_id++);
        if (!attempt.merge) {
            patched.push_back(std::move(incoming));
            continue;
        }
        auto& merged = patched.back();
        merged.ref_start = std::min(merged.ref_start, incoming.ref_start);
        merged.ref_end = std::max(merged.ref_end, incoming.ref_end);
        merged.query_start = std::min(merged.query_start, incoming.query_start);
        merged.query_end = std::max(merged.query_end, incoming.query_end);
        merged.score += incoming.score;
        merged.patch_count += incoming.patch_count + 1;
        merged.trim_count += incoming.trim_count;
        merged.anchor_indices.insert(merged.anchor_indices.end(),
                                     incoming.anchor_indices.begin(),
                                     incoming.anchor_indices.end());
        merged.source = attempt.classification;
        ++patch_count;
    }
    return patched;
}

}  // namespace alignment_detail
}  // namespace vallescope2
