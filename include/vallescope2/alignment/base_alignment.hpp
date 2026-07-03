#pragma once

#include <cstdint>
#include <filesystem>
#include <vector>

namespace vallescope2 {

struct BaseAlignmentParameters {
    std::uint32_t anchor_length = 50;
    std::uint32_t max_bundle_align_bp = 50000;
    std::uint64_t max_fallback_cells = 25000000;
    std::uint32_t max_patch_gap_bp = 70000;
    std::uint32_t patch_window_bp = 1000;
    double min_patch_identity = 0.85;
    std::uint32_t max_wfa_memory_gb = 64;
    double bundle_trim_overlap = 0.01;
};

struct BaseAlignmentResult {
    std::uint64_t loaded_bundle_count = 0;
    std::uint64_t post_patch_bundle_count = 0;
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
};

BaseAlignmentResult align_chain_bundles(
    const std::vector<std::filesystem::path>& chain_files,
    const std::vector<std::filesystem::path>& chain_anchor_files,
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& fasta,
    const std::filesystem::path& fasta_index,
    const std::filesystem::path& paf_output,
    const std::filesystem::path& metadata_output,
    const BaseAlignmentParameters& parameters);

}  // namespace vallescope2
