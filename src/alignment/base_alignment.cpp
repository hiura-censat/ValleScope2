#include "vallescope2/alignment/base_alignment.hpp"

#include "base_alignment_internal.hpp"

#include <htslib/faidx.h>

#include <fstream>
#include <memory>
#include <stdexcept>
#include <utility>

namespace vallescope2 {

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
        metadata_output.parent_path() / "patch_extensions.tsv",
        result.patch_count);
    result.post_patch_bundle_count = bundles.size();
    bundles = final_trim_bundles(
        std::move(bundles), parameters, result.bundle_trim_count);
    result.bundle_count = bundles.size();
    result.patched_bundle_count =
        pre_patch_bundle_count > result.post_patch_bundle_count
            ? pre_patch_bundle_count - result.post_patch_bundle_count
            : 0;

    std::ofstream paf(paf_output);
    if (!paf) throw std::runtime_error("cannot create bundle alignment PAF");
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
        const auto ref_interval = expand_interval(
            bundle.ref_start, bundle.ref_end, ref_length, parameters.anchor_length);
        const auto query_interval = expand_interval(
            bundle.query_start, bundle.query_end, query_length, parameters.anchor_length);
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

        paf << bundle.sequence_b << '\t' << query_length << '\t'
            << query_interval.start << '\t' << query_interval.end << '\t'
            << strand_summary.effective_strand << '\t' << bundle.sequence_a << '\t'
            << ref_length << '\t' << ref_interval.start << '\t'
            << ref_interval.end << '\t' << alignment.matches << '\t'
            << alignment.block_length << "\t60"
            << "\tcg:Z:" << alignment.cigar
            << "\tch:i:" << bundle.chain_id
            << "\tbt:Z:" << bundle.source
            << "\tpg:i:" << bundle.patch_count
            << "\ttr:i:" << bundle.trim_count
            << "\tam:Z:" << alignment_mode
            << "\tos:Z:" << bundle.strand
            << "\tes:Z:" << strand_summary.effective_strand
            << "\taf:i:" << strand_summary.forward_count
            << "\tar:i:" << strand_summary.reverse_count
            << "\tai:i:" << strand_summary.incompatible_count
            << "\tas:i:" << alignment.score
            << "\tsa:Z:" << bundle.sample_a
            << "\tsb:Z:" << bundle.sample_b
            << '\n';
        ++result.aligned_bundle_count;
    }
    write_metadata(metadata_output, parameters, result, chain_files, fasta);
    return result;
}

}  // namespace vallescope2
