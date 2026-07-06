#include "base_alignment_internal.hpp"

#include "bindings/cpp/WFAligner.hpp"
#include "vallescope2/context/anchor_grouping.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
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

struct PatchExtensionStep {
    std::uint64_t ref_consumed = 0;
    std::uint64_t query_consumed = 0;
    double identity = 0.0;
    std::string cigar;
    bool success = false;
};

struct PatchAttempt {
    std::string classification = "unresolved";
    bool merge = false;
    std::uint64_t left_ref_extension = 0;
    std::uint64_t left_query_extension = 0;
    std::uint64_t right_ref_extension = 0;
    std::uint64_t right_query_extension = 0;
    std::string copy_support_side = ".";
    std::uint64_t copy_support_both_count = 0;
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

std::string bundle_track_key(const ChainBundle& bundle) {
    return bundle.sample_a + '\x1f' + bundle.sample_b + '\x1f' +
           bundle.sequence_a + '\x1f' + bundle.sequence_b + '\x1f' +
           bundle.strand;
}

std::string candidate_track_key(const ExtensionCandidate& candidate) {
    return candidate.sample_a + '\x1f' + candidate.sample_b + '\x1f' +
           candidate.sequence_a + '\x1f' + candidate.sequence_b + '\x1f' +
           candidate.strand;
}

std::uint64_t count_open_axis_both_support(
    const std::vector<ExtensionCandidate>& candidates,
    const ChainBundle& left,
    const ChainBundle& right,
    const std::string& side,
    const std::uint64_t left_ref_extension,
    const std::uint64_t left_query_extension,
    const std::uint64_t right_ref_extension,
    const std::uint64_t right_query_extension) {
    const auto key = bundle_track_key(left);
    std::uint64_t count = 0;
    for (const auto& candidate : candidates) {
        if (candidate_track_key(candidate) != key) continue;
        if (candidate.support_direction != "both") continue;
        if (side == "ref") {
            const auto start = left.ref_end + left_ref_extension;
            const auto end = right.ref_start > right_ref_extension
                                 ? right.ref_start - right_ref_extension
                                 : 0;
            if (start < candidate.ref_center && candidate.ref_center < end) {
                ++count;
            }
        } else if (side == "query") {
            if (left.strand == '+') {
                const auto start = left.query_end + left_query_extension;
                const auto end = right.query_start > right_query_extension
                                     ? right.query_start - right_query_extension
                                     : 0;
                if (start < candidate.query_center &&
                    candidate.query_center < end) {
                    ++count;
                }
            } else {
                const auto start = right.query_end + right_query_extension;
                const auto end = left.query_start > left_query_extension
                                     ? left.query_start - left_query_extension
                                     : 0;
                if (start < candidate.query_center &&
                    candidate.query_center < end) {
                    ++count;
                }
            }
        }
    }
    return count;
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

    wfa::WFAlignerGapAffine aligner(
        4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    const auto max_memory =
        static_cast<std::uint64_t>(parameters.max_wfa_memory_gb) * 1024ULL *
        1024ULL * 1024ULL;
    aligner.setMaxMemory(max_memory, max_memory);
    const auto status = aligner.alignExtension(
        ref.c_str(), static_cast<int>(ref.size()), query.c_str(),
        static_cast<int>(query.size()));
    if (status < 0) return best;
    const std::string cigar = aligner.getCIGAR(true);
    best.cigar = cigar;
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
    while (begin < cigar.size()) {
        std::size_t end = begin;
        while (end < cigar.size() &&
               std::isdigit(static_cast<unsigned char>(cigar[end]))) {
            ++end;
        }
        if (end == begin || end >= cigar.size()) break;
        const auto length = parse_u64(
            cigar.substr(begin, end - begin), "cigar_length");
        const char op = cigar[end];
        if ((op == 'I' || op == 'D') &&
            length > parameters.max_patch_indel_bp) {
            break;
        }

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

PatchAttempt evaluate_patch_gap(faidx_t* index,
                                const ChainBundle& left,
                                const ChainBundle& right,
                                const std::vector<ExtensionCandidate>& candidates,
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
    if (ref_gap == 0 && query_gap == 0) {
        attempt.classification = "clean_patch";
        attempt.merge = true;
        return attempt;
    }

    auto write_step = [&](const std::uint64_t step_index,
                          const char* side,
                          const Interval& ref_interval,
                          const Interval& query_interval,
                          const PatchExtensionStep& step,
                          const char* status,
                          const std::uint64_t left_ref_extension_value,
                          const std::uint64_t left_query_extension_value,
                          const std::uint64_t right_ref_extension_value,
                          const std::uint64_t right_query_extension_value) {
        patch_report
            << patch_id << '\t' << step_index << '\t' << side << '\t'
            << status << '\t' << left.chain_id << '\t' << right.chain_id
            << '\t' << left.sample_a << '\t' << left.sample_b << '\t'
            << left.sequence_a << '\t' << left.sequence_b << '\t'
            << left.strand << '\t' << ref_gap << '\t' << query_gap << '\t'
            << ref_interval.start << '\t' << ref_interval.end << '\t'
            << query_interval.start << '\t' << query_interval.end << '\t'
            << (ref_interval.end - ref_interval.start) << '\t'
            << (query_interval.end - query_interval.start) << '\t'
            << step.ref_consumed << '\t' << step.query_consumed << '\t'
            << step.identity << '\t' << (step.success ? 1 : 0) << '\t'
            << step.cigar << '\t'
            << left_ref_extension_value << '\t'
            << left_query_extension_value << '\t'
            << right_ref_extension_value << '\t'
            << right_query_extension_value << '\t'
            << left.ref_start << '\t' << left.ref_end << '\t'
            << left.query_start << '\t' << left.query_end << '\t'
            << right.ref_start << '\t' << right.ref_end << '\t'
            << right.query_start << '\t' << right.query_end << '\n';
    };

    std::uint64_t left_ref_extension = 0;
    std::uint64_t left_query_extension = 0;
    std::uint64_t right_ref_extension = 0;
    std::uint64_t right_query_extension = 0;
    std::uint64_t step_index = 0;

    auto ref_gap_is_closed = [&]() {
        return left_ref_extension + right_ref_extension >= ref_gap;
    };
    auto query_gap_is_closed = [&]() {
        return left_query_extension + right_query_extension >= query_gap;
    };
    auto gap_is_closed = [&]() {
        return ref_gap_is_closed() && query_gap_is_closed();
    };
    auto remaining_ref_gap = [&]() {
        return ref_gap - std::min(ref_gap,
                                  left_ref_extension + right_ref_extension);
    };
    auto remaining_query_gap = [&]() {
        return query_gap -
               std::min(query_gap, left_query_extension + right_query_extension);
    };
    auto limit_step_at_bridge_breakpoint =
        [&](PatchExtensionStep& step,
            const std::uint64_t current_left_ref_extension,
            const std::uint64_t current_left_query_extension,
            const std::uint64_t current_right_ref_extension,
            const std::uint64_t current_right_query_extension) {
            if (!step.success) return false;
            const auto predicted_ref =
                current_left_ref_extension + current_right_ref_extension +
                step.ref_consumed;
            const auto predicted_query =
                current_left_query_extension + current_right_query_extension +
                step.query_consumed;
            const bool ref_would_close =
                predicted_ref >= ref_gap && predicted_query < query_gap;
            const bool query_would_close =
                predicted_query >= query_gap && predicted_ref < ref_gap;
            if (!ref_would_close && !query_would_close) return false;

            const std::string open_side = ref_would_close ? "query" : "ref";
            const auto support = count_open_axis_both_support(
                candidates, left, right, open_side,
                current_left_ref_extension, current_left_query_extension,
                current_right_ref_extension, current_right_query_extension);
            attempt.copy_support_side = open_side;
            attempt.copy_support_both_count = support;
            if (support >= parameters.min_copy_support_anchors) return false;

            const auto margin =
                static_cast<std::uint64_t>(parameters.anchor_length);
            const auto ref_limit = ref_gap > margin ? ref_gap - margin : 0;
            const auto query_limit =
                query_gap > margin ? query_gap - margin : 0;
            double fraction = 1.0;
            if (ref_would_close && step.ref_consumed > 0) {
                const auto current_ref =
                    current_left_ref_extension + current_right_ref_extension;
                const auto allowed =
                    ref_limit > current_ref ? ref_limit - current_ref : 0;
                fraction = std::min(
                    fraction,
                    static_cast<double>(allowed) /
                        static_cast<double>(step.ref_consumed));
            }
            if (query_would_close && step.query_consumed > 0) {
                const auto current_query = current_left_query_extension +
                                           current_right_query_extension;
                const auto allowed = query_limit > current_query
                                         ? query_limit - current_query
                                         : 0;
                fraction = std::min(
                    fraction,
                    static_cast<double>(allowed) /
                        static_cast<double>(step.query_consumed));
            }
            if (fraction < 0.0) fraction = 0.0;
            step.ref_consumed = static_cast<std::uint64_t>(
                std::floor(static_cast<double>(step.ref_consumed) * fraction));
            step.query_consumed = static_cast<std::uint64_t>(
                std::floor(static_cast<double>(step.query_consumed) * fraction));
            if (step.ref_consumed == 0 && step.query_consumed == 0) {
                step.success = false;
            }
            return true;
        };
    auto finalize_attempt = [&]() {
        attempt.left_ref_extension = left_ref_extension;
        attempt.left_query_extension = left_query_extension;
        attempt.right_ref_extension = right_ref_extension;
        attempt.right_query_extension = right_query_extension;
        const auto remaining_ref = remaining_ref_gap();
        const auto remaining_query = remaining_query_gap();
        if (remaining_ref == 0 && remaining_query == 0) {
            attempt.classification = "clean_patch";
            attempt.merge = true;
        } else if (remaining_query <= parameters.patch_window_bp &&
                   remaining_ref > remaining_query + parameters.max_patch_indel_bp) {
            attempt.classification = "deletion_bridge";
            attempt.merge = false;
        } else if (remaining_ref <= parameters.patch_window_bp &&
                   remaining_query > remaining_ref + parameters.max_patch_indel_bp) {
            attempt.classification = "insertion_bridge";
            attempt.merge = false;
        } else {
            attempt.classification = "unresolved";
        }
        return attempt;
    };

    while (!gap_is_closed()) {
        const auto remaining_ref = remaining_ref_gap();
        const auto remaining_query = remaining_query_gap();
        if (remaining_ref == 0 || remaining_query == 0) return finalize_attempt();

        const auto left_ref_window =
            std::min<std::uint64_t>(
                static_cast<std::uint64_t>(parameters.patch_window_bp) +
                    parameters.patch_window_slack_bp,
                remaining_ref);
        const auto left_query_window =
            std::min<std::uint64_t>(
                static_cast<std::uint64_t>(parameters.patch_window_bp) +
                    parameters.patch_window_slack_bp,
                remaining_query);
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
        auto left_step = extend_patch_window_wfa(
            index, left, left_ref, left_query, false, parameters);
        const bool left_bridge_limited = limit_step_at_bridge_breakpoint(
            left_step,
            left_ref_extension,
            left_query_extension,
            right_ref_extension,
            right_query_extension);
        write_step(step_index++, "left", left_ref, left_query, left_step,
                   left_step.success
                       ? (left_bridge_limited ? "bridge_limited" : "accepted")
                       : "failed",
                   left_ref_extension + left_step.ref_consumed,
                   left_query_extension + left_step.query_consumed,
                   right_ref_extension, right_query_extension);
        if (!left_step.success) {
            return finalize_attempt();
        }
        left_ref_extension += left_step.ref_consumed;
        left_query_extension += left_step.query_consumed;
        if (left_bridge_limited) return finalize_attempt();
        if (gap_is_closed()) return finalize_attempt();
        if (ref_gap_is_closed() != query_gap_is_closed()) {
            return finalize_attempt();
        }

        const auto after_left_ref =
            ref_gap - std::min(ref_gap, left_ref_extension + right_ref_extension);
        const auto after_left_query =
            query_gap - std::min(query_gap,
                                 left_query_extension + right_query_extension);
        const auto right_ref_window =
            std::min<std::uint64_t>(
                static_cast<std::uint64_t>(parameters.patch_window_bp) +
                    parameters.patch_window_slack_bp,
                after_left_ref);
        const auto right_query_window =
            std::min<std::uint64_t>(
                static_cast<std::uint64_t>(parameters.patch_window_bp) +
                    parameters.patch_window_slack_bp,
                after_left_query);
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
        auto right_step = extend_patch_window_wfa(
            index, right, right_ref, right_query, true, parameters);
        const bool right_bridge_limited = limit_step_at_bridge_breakpoint(
            right_step,
            left_ref_extension,
            left_query_extension,
            right_ref_extension,
            right_query_extension);
        write_step(step_index++, "right", right_ref, right_query, right_step,
                   right_step.success
                       ? (right_bridge_limited ? "bridge_limited" : "accepted")
                       : "failed",
                   left_ref_extension, left_query_extension,
                   right_ref_extension + right_step.ref_consumed,
                   right_query_extension + right_step.query_consumed);
        if (!right_step.success) {
            return finalize_attempt();
        }
        right_ref_extension += right_step.ref_consumed;
        right_query_extension += right_step.query_consumed;
        if (right_bridge_limited) return finalize_attempt();
        if (ref_gap_is_closed() != query_gap_is_closed()) {
            return finalize_attempt();
        }
    }
    return finalize_attempt();
}

void extend_left_bundle(ChainBundle& bundle, const PatchAttempt& attempt) {
    bundle.ref_end += attempt.left_ref_extension;
    if (bundle.strand == '+') {
        bundle.query_end += attempt.left_query_extension;
    } else {
        bundle.query_start =
            bundle.query_start > attempt.left_query_extension
                ? bundle.query_start - attempt.left_query_extension
                : 0;
    }
}

void extend_right_bundle(ChainBundle& bundle, const PatchAttempt& attempt) {
    bundle.ref_start =
        bundle.ref_start > attempt.right_ref_extension
            ? bundle.ref_start - attempt.right_ref_extension
            : 0;
    if (bundle.strand == '+') {
        bundle.query_start =
            bundle.query_start > attempt.right_query_extension
                ? bundle.query_start - attempt.right_query_extension
                : 0;
    } else {
        bundle.query_end += attempt.right_query_extension;
    }
}

std::vector<ChainBundle> patch_adjacent_bundles(
    std::vector<ChainBundle> bundles,
    faidx_t* index,
    const std::vector<ExtensionCandidate>& candidates,
    const BaseAlignmentParameters& parameters,
    const std::filesystem::path& patch_extension_output,
    std::uint64_t& patch_count) {
    patch_count = 0;
    std::ofstream patch_report(patch_extension_output);
    if (!patch_report) {
        throw std::runtime_error("cannot create patch extension report");
    }
    patch_report
        << "patch_id\tstep_index\tside\tstatus\tleft_chain_id"
        << "\tright_chain_id\tsample_a\tsample_b\tsequence_a\tsequence_b"
        << "\tstrand\tref_gap\tquery_gap\tref_window_start"
        << "\tref_window_end\tquery_window_start\tquery_window_end"
        << "\tref_window_len\tquery_window_len\tref_consumed"
        << "\tquery_consumed\tidentity\tsuccess\tcigar"
        << "\tleft_ref_extension"
        << "\tleft_query_extension\tright_ref_extension"
        << "\tright_query_extension\tleft_ref_start\tleft_ref_end"
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
            index, patched.back(), incoming, candidates, parameters, patch_report,
            patch_id++);
        if (!attempt.merge) {
            extend_left_bundle(patched.back(), attempt);
            extend_right_bundle(incoming, attempt);
            patched.push_back(std::move(incoming));
            continue;
        }
        extend_left_bundle(patched.back(), attempt);
        extend_right_bundle(incoming, attempt);
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
