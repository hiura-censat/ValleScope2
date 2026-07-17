#include "base_alignment_internal.hpp"

#include <fstream>
#include <stdexcept>

namespace vallescope2 {
namespace alignment_detail {

void write_metadata(const std::filesystem::path& path,
                    const BaseAlignmentParameters& parameters,
                    const BaseAlignmentResult& result,
                    const std::vector<std::filesystem::path>& chain_files,
                    const std::filesystem::path& fasta) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create base alignment metadata");
    output << "{\n"
           << "  \"input_chains\": [";
    for (std::size_t i = 0; i < chain_files.size(); ++i) {
        if (i) output << ", ";
        output << "\"" << chain_files[i].string() << "\"";
    }
    output << "],\n"
           << "  \"fasta\": \"" << fasta.string() << "\",\n"
           << "  \"alignment_mode\": \"anchor_guided\",\n"
           << "  \"aligner\": \"WFA2-lib\",\n"
           << "  \"anchor_length\": " << parameters.anchor_length << ",\n"
           << "  \"max_bundle_align_bp\": " << parameters.max_bundle_align_bp << ",\n"
           << "  \"max_fallback_cells\": " << parameters.max_fallback_cells << ",\n"
           << "  \"max_patch_gap_bp\": " << parameters.max_patch_gap_bp << ",\n"
           << "  \"patch_flank_bp\": " << parameters.patch_flank_bp << ",\n"
           << "  \"min_patch_identity\": " << parameters.min_patch_identity << ",\n"
           << "  \"min_patch_long_indel_bp\": "
           << parameters.min_patch_long_indel_bp << ",\n"
           << "  \"min_patch_rescue_flank_bp\": "
           << parameters.min_patch_rescue_flank_bp << ",\n"
           << "  \"min_patch_rescue_flank_identity\": "
           << parameters.min_patch_rescue_flank_identity << ",\n"
           << "  \"max_patch_rescue_extra_indel_bp\": "
           << parameters.max_patch_rescue_extra_indel_bp << ",\n"
           << "  \"patch_quality_window_bp\": "
           << parameters.patch_quality_window_bp << ",\n"
           << "  \"min_patch_endpoint_identity\": "
           << parameters.min_patch_endpoint_identity << ",\n"
           << "  \"max_patch_short_error_density\": "
           << parameters.max_patch_short_error_density << ",\n"
           << "  \"min_patch_multi_segment_identity\": "
           << parameters.min_patch_multi_segment_identity << ",\n"
           << "  \"max_patch_multi_short_indel_bp\": "
           << parameters.max_patch_multi_short_indel_bp << ",\n"
           << "  \"max_patch_multi_short_error_density\": "
           << parameters.max_patch_multi_short_error_density << ",\n"
           << "  \"patch_multi_event_allow_deletions\": "
           << (parameters.patch_multi_event_allow_deletions ? "true" : "false")
           << ",\n"
           << "  \"max_wfa_memory_gb\": " << parameters.max_wfa_memory_gb << ",\n"
           << "  \"bundle_trim_overlap\": " << parameters.bundle_trim_overlap << ",\n"
           << "  \"chain_extension\": "
           << (parameters.chain_extension ? "true" : "false") << ",\n"
           << "  \"chain_extension_predecessors\": "
           << parameters.chain_extension_predecessors << ",\n"
           << "  \"max_chain_extension_bp\": "
           << parameters.max_chain_extension_bp << ",\n"
           << "  \"min_chain_extension_anchors\": "
           << parameters.min_chain_extension_anchors << ",\n"
           << "  \"min_chain_extension_score\": "
           << parameters.min_chain_extension_score << ",\n"
           << "  \"min_copy_support_anchors\": "
           << parameters.min_copy_support_anchors << ",\n"
           << "  \"loaded_bundle_count\": " << result.loaded_bundle_count << ",\n"
           << "  \"post_extension_bundle_count\": "
           << result.post_extension_bundle_count << ",\n"
           << "  \"extension_count\": " << result.extension_count << ",\n"
           << "  \"extension_partial_count\": "
           << result.extension_partial_count << ",\n"
           << "  \"extension_query_dup_branch_count\": "
           << result.extension_query_dup_branch_count << ",\n"
           << "  \"extension_ref_dup_branch_count\": "
           << result.extension_ref_dup_branch_count << ",\n"
           << "  \"extension_query_insertion_count\": "
           << result.extension_query_insertion_count << ",\n"
           << "  \"extension_query_deletion_count\": "
           << result.extension_query_deletion_count << ",\n"
           << "  \"extension_unresolved_gap_count\": "
           << result.extension_unresolved_gap_count << ",\n"
           << "  \"post_patch_bundle_count\": " << result.post_patch_bundle_count << ",\n"
           << "  \"bundle_trim_count\": " << result.bundle_trim_count << ",\n"
           << "  \"bundle_count\": " << result.bundle_count << ",\n"
           << "  \"patched_bundle_count\": " << result.patched_bundle_count << ",\n"
           << "  \"patch_count\": " << result.patch_count << ",\n"
           << "  \"aligned_bundle_count\": " << result.aligned_bundle_count << ",\n"
           << "  \"anchor_guided_alignment_count\": "
           << result.anchor_guided_alignment_count << ",\n"
           << "  \"whole_bundle_alignment_count\": "
           << result.whole_bundle_alignment_count << ",\n"
           << "  \"strand_flipped_bundle_count\": "
           << result.strand_flipped_bundle_count << ",\n"
           << "  \"strand_conflict_bundle_count\": "
           << result.strand_conflict_bundle_count << ",\n"
           << "  \"forward_anchor_match_count\": "
           << result.forward_anchor_match_count << ",\n"
           << "  \"reverse_anchor_match_count\": "
           << result.reverse_anchor_match_count << ",\n"
           << "  \"incompatible_anchor_match_count\": "
           << result.incompatible_anchor_match_count << ",\n"
           << "  \"no_valid_anchor_bundle_count\": "
           << result.no_valid_anchor_bundle_count << ",\n"
           << "  \"segment_alignment_failed_count\": "
           << result.segment_alignment_failed_count << ",\n"
           << "  \"skipped_bundle_count\": " << result.skipped_bundle_count << ",\n"
           << "  \"paf_normalization_event_count\": "
           << result.paf_normalization_event_count << ",\n"
           << "  \"paf_normalization_clipped_block_count\": "
           << result.paf_normalization_clipped_block_count << ",\n"
           << "  \"paf_normalization_unresolved_overlap_count\": "
           << result.paf_normalization_unresolved_overlap_count << ",\n"
           << "  \"paf_normalization_realign_count\": "
           << result.paf_normalization_realign_count << "\n"
           << "}\n";
}

}  // namespace alignment_detail
}  // namespace vallescope2
