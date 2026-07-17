#pragma once

#include "vallescope2/bundle/chaining.hpp"

#include <cstdint>
#include <filesystem>
#include <vector>

namespace vallescope2 {

struct BaseAlignmentParameters {
    std::uint32_t anchor_length = 50;
    std::uint32_t max_bundle_align_bp = 1000000;
    std::uint64_t max_fallback_cells = 25000000;
    std::uint32_t max_patch_gap_bp = 70000;
    std::uint32_t patch_flank_bp = 500;
    double min_patch_identity = 0.95;
    std::uint32_t min_patch_long_indel_bp = 500;
    std::uint32_t min_patch_rescue_flank_bp = 200;
    double min_patch_rescue_flank_identity = 0.95;
    std::uint32_t max_patch_rescue_extra_indel_bp = 100;
    std::uint32_t patch_quality_window_bp = 500;
    double min_patch_endpoint_identity = 0.85;
    double max_patch_short_error_density = 0.15;
    double min_patch_multi_segment_identity = 1.0;
    std::uint32_t max_patch_multi_short_indel_bp = 0;
    double max_patch_multi_short_error_density = 0.0;
    bool patch_multi_event_allow_deletions = false;
    std::uint32_t max_wfa_memory_gb = 64;
    double bundle_trim_overlap = 0.01;
    bool chain_extension = true;
    std::uint32_t chain_extension_predecessors = 50;
    std::uint32_t max_chain_extension_bp = 70000;
    std::uint32_t max_chain_gap = 300;
    double chain_max_gap_ratio = 1.05;
    GapCostModel gap_cost_model = GapCostModel::absolute;
    double gap_weight = 0.002;
    std::uint32_t gap_unit = 10;
    std::uint32_t min_chain_extension_anchors = 5;
    double min_chain_extension_score = 100.0;
    std::uint32_t min_copy_support_anchors = 1;
};

struct BaseAlignmentResult {
    std::uint64_t loaded_bundle_count = 0;
    std::uint64_t post_patch_bundle_count = 0;
    std::uint64_t post_extension_bundle_count = 0;
    std::uint64_t extension_count = 0;
    std::uint64_t extension_partial_count = 0;
    std::uint64_t extension_query_dup_branch_count = 0;
    std::uint64_t extension_ref_dup_branch_count = 0;
    std::uint64_t extension_query_insertion_count = 0;
    std::uint64_t extension_query_deletion_count = 0;
    std::uint64_t extension_unresolved_gap_count = 0;
    std::uint64_t bundle_trim_count = 0;
    std::uint64_t bundle_count = 0;
    std::uint64_t patched_bundle_count = 0;
    std::uint64_t patch_count = 0;
    std::uint64_t aligned_bundle_count = 0;
    std::uint64_t anchor_guided_alignment_count = 0;
    std::uint64_t whole_bundle_alignment_count = 0;
    std::uint64_t strand_flipped_bundle_count = 0;
    std::uint64_t strand_conflict_bundle_count = 0;
    std::uint64_t forward_anchor_match_count = 0;
    std::uint64_t reverse_anchor_match_count = 0;
    std::uint64_t incompatible_anchor_match_count = 0;
    std::uint64_t no_valid_anchor_bundle_count = 0;
    std::uint64_t segment_alignment_failed_count = 0;
    std::uint64_t skipped_bundle_count = 0;
    std::uint64_t paf_normalization_event_count = 0;
    std::uint64_t paf_normalization_clipped_block_count = 0;
    std::uint64_t paf_normalization_unresolved_overlap_count = 0;
    std::uint64_t paf_normalization_realign_count = 0;
};

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
    const BaseAlignmentParameters& parameters);

}  // namespace vallescope2
