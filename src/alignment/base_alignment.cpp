#include "vallescope2/alignment/base_alignment.hpp"

#include "base_alignment_internal.hpp"

#include <htslib/faidx.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <utility>

namespace vallescope2 {
namespace {

enum class BundleAxis { ref, query };

struct PafRecord {
    std::string query_name;
    std::uint64_t query_length = 0;
    alignment_detail::Interval query_interval;
    char strand = '+';
    std::string target_name;
    std::uint64_t target_length = 0;
    alignment_detail::Interval target_interval;
    alignment_detail::Alignment alignment;
    std::uint64_t chain_id = 0;
    std::string bundle_source;
    std::uint64_t patch_count = 0;
    std::uint64_t trim_count = 0;
    std::string alignment_mode;
    char original_strand = '+';
    std::uint64_t forward_count = 0;
    std::uint64_t reverse_count = 0;
    std::uint64_t incompatible_count = 0;
    std::string sample_a;
    std::string sample_b;
    bool needs_realign = false;
};

alignment_detail::Interval bundle_alignment_interval(
    const alignment_detail::ChainBundle& bundle,
    const alignment_detail::AnchorStore& anchor_store,
    const std::uint64_t sequence_length,
    const BundleAxis axis) {
    auto start = axis == BundleAxis::ref ? bundle.ref_start : bundle.query_start;
    auto end = axis == BundleAxis::ref ? bundle.ref_end : bundle.query_end;
    if (end < start) std::swap(start, end);

    for (const auto anchor_index : bundle.anchor_indices) {
        if (anchor_index >= anchor_store.anchors.size()) continue;
        const auto& anchor = anchor_store.anchors[anchor_index];
        const auto anchor_start =
            axis == BundleAxis::ref ? anchor.ref_start : anchor.query_start;
        const auto anchor_end =
            axis == BundleAxis::ref ? anchor.ref_end : anchor.query_end;
        start = std::min(start, anchor_start);
        end = std::max(end, anchor_end);
    }

    start = std::min(start, sequence_length);
    end = std::min(end, sequence_length);
    if (end <= start && sequence_length > 0) {
        start = std::min(start, sequence_length - 1);
        end = start + 1;
    }
    return {start, end};
}

alignment_detail::Interval query_oriented_interval(const PafRecord& record) {
    if (record.strand == '+') return record.query_interval;
    return {record.query_length - record.query_interval.end,
            record.query_length - record.query_interval.start};
}

bool has_valid_intervals(const PafRecord& record) {
    return record.target_interval.start < record.target_interval.end &&
           record.query_interval.start < record.query_interval.end;
}

void realign_paf_record(PafRecord& record,
                        faidx_t* index,
                        const BaseAlignmentParameters& parameters) {
    using namespace alignment_detail;
    const auto ref_sequence =
        fetch_interval(index, record.target_name, record.target_interval);
    const auto query_sequence = fetch_oriented_query_interval(
        index, record.query_name, query_oriented_interval(record),
        record.query_length, record.strand);
    record.alignment = align_segment(ref_sequence, query_sequence, parameters);
    record.needs_realign = false;
}

void write_paf_record(std::ostream& output, const PafRecord& record) {
    output << record.query_name << '\t' << record.query_length << '\t'
           << record.query_interval.start << '\t' << record.query_interval.end
           << '\t' << record.strand << '\t' << record.target_name << '\t'
           << record.target_length << '\t' << record.target_interval.start
           << '\t' << record.target_interval.end << '\t'
           << record.alignment.matches << '\t' << record.alignment.block_length
           << "\t60"
           << "\tcg:Z:" << record.alignment.cigar
           << "\tch:i:" << record.chain_id
           << "\tbt:Z:" << record.bundle_source
           << "\tpg:i:" << record.patch_count
           << "\ttr:i:" << record.trim_count
           << "\tam:Z:" << record.alignment_mode
           << "\tos:Z:" << record.original_strand
           << "\tes:Z:" << record.strand
           << "\taf:i:" << record.forward_count
           << "\tar:i:" << record.reverse_count
           << "\tai:i:" << record.incompatible_count
           << "\tas:i:" << record.alignment.score
           << "\tsa:Z:" << record.sample_a
           << "\tsb:Z:" << record.sample_b
           << '\n';
}

void normalize_paf_records(std::vector<PafRecord>& records,
                           faidx_t* index,
                           const BaseAlignmentParameters& parameters,
                           const std::filesystem::path& event_output,
                           BaseAlignmentResult& result) {
    using Key = std::tuple<std::string, std::string, char>;
    std::map<Key, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < records.size(); ++i) {
        groups[{records[i].target_name, records[i].query_name,
                records[i].strand}]
            .push_back(i);
    }

    std::ofstream events(event_output);
    if (!events) {
        throw std::runtime_error("cannot create PAF normalization report");
    }
    events
        << "action\treason\ttarget\tquery\tstrand\tleft_chain_id"
        << "\tright_chain_id\tleft_ref_start\tleft_ref_end"
        << "\tright_ref_start\tright_ref_end\tleft_query_start"
        << "\tleft_query_end\tright_query_start\tright_query_end"
        << "\tref_overlap\tquery_overlap\n";

    for (auto& [key, indices] : groups) {
        std::sort(indices.begin(), indices.end(),
                  [&records](const std::size_t left, const std::size_t right) {
                      const auto& a = records[left];
                      const auto& b = records[right];
                      return std::tie(a.target_interval.start,
                                      a.query_interval.start,
                                      a.target_interval.end) <
                             std::tie(b.target_interval.start,
                                      b.query_interval.start,
                                      b.target_interval.end);
                  });

        for (std::size_t i = 1; i < indices.size(); ++i) {
            auto& left = records[indices[i - 1]];
            auto& right = records[indices[i]];
            const auto left_query_oriented = query_oriented_interval(left);
            const auto right_query_oriented = query_oriented_interval(right);
            const std::uint64_t ref_overlap =
                left.target_interval.end > right.target_interval.start
                    ? left.target_interval.end - right.target_interval.start
                    : 0;
            const std::uint64_t query_overlap =
                left_query_oriented.end > right_query_oriented.start
                    ? left_query_oriented.end - right_query_oriented.start
                    : 0;
            if (ref_overlap == 0 && query_overlap == 0) continue;

            ++result.paf_normalization_event_count;
            const bool anchor_length_overlap =
                (ref_overlap == parameters.anchor_length && query_overlap == 0) ||
                (query_overlap == parameters.anchor_length && ref_overlap == 0);
            const char* reason =
                anchor_length_overlap ? "anchor_length_overlap" : "non_normalized_overlap";
            constexpr bool clip_left = true;
            const char* action =
                anchor_length_overlap ? "clip_left" : "report_only";

            events << action << '\t' << reason << '\t' << left.target_name
                   << '\t' << left.query_name << '\t' << left.strand << '\t'
                   << left.chain_id << '\t' << right.chain_id << '\t'
                   << left.target_interval.start << '\t'
                   << left.target_interval.end << '\t'
                   << right.target_interval.start << '\t'
                   << right.target_interval.end << '\t'
                   << left.query_interval.start << '\t'
                   << left.query_interval.end << '\t'
                   << right.query_interval.start << '\t'
                   << right.query_interval.end << '\t'
                   << ref_overlap << '\t' << query_overlap << '\n';

            if (!anchor_length_overlap) {
                ++result.paf_normalization_unresolved_overlap_count;
                continue;
            }

            auto& clipped = clip_left ? left : right;
            const auto old_target = clipped.target_interval;
            const auto old_query = clipped.query_interval;
            if (clip_left) {
                if (ref_overlap > 0) {
                    clipped.target_interval.end = right.target_interval.start;
                }
                if (query_overlap > 0) {
                    if (clipped.strand == '+') {
                        clipped.query_interval.end = right.query_interval.start;
                    } else {
                        clipped.query_interval.start = right.query_interval.end;
                    }
                }
            } else {
                if (ref_overlap > 0) {
                    clipped.target_interval.start = left.target_interval.end;
                }
                if (query_overlap > 0) {
                    if (clipped.strand == '+') {
                        clipped.query_interval.start = left.query_interval.end;
                    } else {
                        clipped.query_interval.end = left.query_interval.start;
                    }
                }
            }
            if (!has_valid_intervals(clipped)) {
                clipped.target_interval = old_target;
                clipped.query_interval = old_query;
                ++result.paf_normalization_unresolved_overlap_count;
                continue;
            }
            clipped.needs_realign = true;
            ++result.paf_normalization_clipped_block_count;
        }
    }

    for (auto& record : records) {
        if (!record.needs_realign) continue;
        realign_paf_record(record, index, parameters);
        ++result.paf_normalization_realign_count;
    }
}

}  // namespace

BaseAlignmentResult align_chain_bundles(
    const std::vector<std::filesystem::path>& chain_files,
    const std::vector<std::filesystem::path>& chain_anchor_files,
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& correspondences,
    const std::filesystem::path& assignments,
    const std::filesystem::path& fasta,
    const std::filesystem::path& fasta_index,
    const std::filesystem::path& paf_output,
    const std::filesystem::path& metadata_output,
    const BaseAlignmentParameters& parameters) {
    using namespace alignment_detail;

    AnchorStore anchor_store;
    const auto loaded_bundles = load_chain_files(
        chain_files, chain_anchor_files, grouped_anchors, anchor_store);
    std::vector<ExtensionCandidate> extension_candidates;
    if (parameters.chain_extension) {
        extension_candidates = load_extension_candidates(
            correspondences, assignments, grouped_anchors);
    }
    BaseAlignmentResult result;
    result.loaded_bundle_count = loaded_bundles.size();

    std::unique_ptr<faidx_t, FaiDeleter> index(
        fai_load3(fasta.c_str(), fasta_index.c_str(), nullptr, 0));
    if (!index) throw std::runtime_error("cannot load FASTA index for base alignment");

    auto bundles = extend_adjacent_bundles(
        loaded_bundles, extension_candidates, anchor_store, parameters,
        metadata_output.parent_path() / "chain_extensions.tsv",
        metadata_output.parent_path() / "chain_extension_anchors.tsv", result);
    result.post_extension_bundle_count = bundles.size();

    const auto pre_patch_bundle_count = bundles.size();
    bundles = patch_adjacent_bundles(
        std::move(bundles), index.get(), extension_candidates, parameters,
        metadata_output.parent_path() / "patch_intervals.tsv",
        result.patch_count);
    result.post_patch_bundle_count = bundles.size();
    bundles = final_trim_bundles(
        std::move(bundles), parameters, result.bundle_trim_count);
    result.bundle_count = bundles.size();
    result.patched_bundle_count =
        pre_patch_bundle_count > result.post_patch_bundle_count
            ? pre_patch_bundle_count - result.post_patch_bundle_count
            : 0;

    std::vector<PafRecord> paf_records;
    std::ofstream strand_conflicts(
        metadata_output.parent_path() / "strand_conflicts.tsv");
    if (!strand_conflicts) {
        throw std::runtime_error("cannot create strand conflict report");
    }
    strand_conflicts
        << "sample_a\tsample_b\tsequence_a\tsequence_b\tchain_id"
        << "\tsource\tassigned_strand\teffective_strand"
        << "\tforward_count\treverse_count\tincompatible_count"
        << "\tref_start\tref_end\tquery_start\tquery_end\n";
    std::ofstream skipped_bundles(
        metadata_output.parent_path() / "skipped_bundles.tsv");
    if (!skipped_bundles) {
        throw std::runtime_error("cannot create skipped bundle report");
    }
    skipped_bundles
        << "reason\tsample_a\tsample_b\tsequence_a\tsequence_b\tchain_id"
        << "\tsource\tassigned_strand\teffective_strand"
        << "\tforward_count\treverse_count\tincompatible_count"
        << "\tref_start\tref_end\tquery_start\tquery_end\n";

    for (const auto& bundle : bundles) {
        const auto ref_length = sequence_length(index.get(), bundle.sequence_a);
        const auto query_length = sequence_length(index.get(), bundle.sequence_b);
        const auto ref_interval = bundle_alignment_interval(
            bundle, anchor_store, ref_length, BundleAxis::ref);
        const auto query_interval = bundle_alignment_interval(
            bundle, anchor_store, query_length, BundleAxis::query);
        const auto ref_span = ref_interval.end - ref_interval.start;
        const auto query_span = query_interval.end - query_interval.start;
        if (ref_span > parameters.max_bundle_align_bp ||
            query_span > parameters.max_bundle_align_bp) {
            skipped_bundles
                << "too_large\t" << bundle.sample_a << '\t'
                << bundle.sample_b << '\t' << bundle.sequence_a << '\t'
                << bundle.sequence_b << '\t' << bundle.chain_id << '\t'
                << bundle.source << '\t' << bundle.strand << "\t.\t.\t.\t.\t"
                << bundle.ref_start << '\t' << bundle.ref_end << '\t'
                << bundle.query_start << '\t' << bundle.query_end << '\n';
            ++result.skipped_bundle_count;
            continue;
        }

        const auto strand_summary = summarize_anchor_strand(
            index.get(), bundle, anchor_store, ref_interval, query_interval);
        result.forward_anchor_match_count += strand_summary.forward_count;
        result.reverse_anchor_match_count += strand_summary.reverse_count;
        result.incompatible_anchor_match_count +=
            strand_summary.incompatible_count;
        if (strand_summary.effective_strand != bundle.strand) {
            ++result.strand_flipped_bundle_count;
        }
        if (strand_summary.forward_count > 0 &&
            strand_summary.reverse_count > 0) {
            ++result.strand_conflict_bundle_count;
            strand_conflicts
                << bundle.sample_a << '\t' << bundle.sample_b << '\t'
                << bundle.sequence_a << '\t' << bundle.sequence_b << '\t'
                << bundle.chain_id << '\t' << bundle.source << '\t'
                << bundle.strand << '\t'
                << strand_summary.effective_strand << '\t'
                << strand_summary.forward_count << '\t'
                << strand_summary.reverse_count << '\t'
                << strand_summary.incompatible_count << '\t'
                << bundle.ref_start << '\t' << bundle.ref_end << '\t'
                << bundle.query_start << '\t' << bundle.query_end << '\n';
        }

        Alignment alignment;
        std::string alignment_mode = "anchor_guided";
        try {
            const auto guided = anchor_guided_alignment(
                index.get(), bundle, anchor_store, ref_interval, query_interval,
                query_length, strand_summary.effective_strand, parameters);
            if (!guided) {
                ++result.no_valid_anchor_bundle_count;
                Interval query_oriented_interval = query_interval;
                if (strand_summary.effective_strand == '-') {
                    query_oriented_interval =
                        {query_length - query_interval.end,
                         query_length - query_interval.start};
                }
                const auto ref_sequence = fetch_interval(
                    index.get(), bundle.sequence_a, ref_interval);
                const auto query_sequence = fetch_oriented_query_interval(
                    index.get(), bundle.sequence_b, query_oriented_interval,
                    query_length, strand_summary.effective_strand);
                alignment = align_segment(ref_sequence, query_sequence,
                                          parameters);
                alignment_mode = "whole_bundle";
                ++result.whole_bundle_alignment_count;
            } else {
                alignment = *guided;
                ++result.anchor_guided_alignment_count;
            }
        } catch (const std::exception&) {
            skipped_bundles
                << "segment_alignment_failed\t" << bundle.sample_a << '\t'
                << bundle.sample_b << '\t' << bundle.sequence_a << '\t'
                << bundle.sequence_b << '\t' << bundle.chain_id << '\t'
                << bundle.source << '\t' << bundle.strand << '\t'
                << strand_summary.effective_strand << '\t'
                << strand_summary.forward_count << '\t'
                << strand_summary.reverse_count << '\t'
                << strand_summary.incompatible_count << '\t'
                << bundle.ref_start << '\t' << bundle.ref_end << '\t'
                << bundle.query_start << '\t' << bundle.query_end << '\n';
            ++result.segment_alignment_failed_count;
            ++result.skipped_bundle_count;
            continue;
        }

        paf_records.push_back(
            PafRecord{bundle.sequence_b,
                      query_length,
                      query_interval,
                      strand_summary.effective_strand,
                      bundle.sequence_a,
                      ref_length,
                      ref_interval,
                      alignment,
                      bundle.chain_id,
                      bundle.source,
                      bundle.patch_count,
                      bundle.trim_count,
                      alignment_mode,
                      bundle.strand,
                      strand_summary.forward_count,
                      strand_summary.reverse_count,
                      strand_summary.incompatible_count,
                      bundle.sample_a,
                      bundle.sample_b});
        ++result.aligned_bundle_count;
    }

    normalize_paf_records(
        paf_records, index.get(), parameters,
        metadata_output.parent_path() / "paf_normalization.tsv", result);

    std::ofstream paf(paf_output);
    if (!paf) throw std::runtime_error("cannot create bundle alignment PAF");
    for (const auto& record : paf_records) {
        write_paf_record(paf, record);
    }
    write_metadata(metadata_output, parameters, result, chain_files, fasta);
    return result;
}

}  // namespace vallescope2
