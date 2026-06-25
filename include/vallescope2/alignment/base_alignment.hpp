#pragma once

#include <cstdint>
#include <filesystem>

namespace vallescope2 {

struct BaseAlignmentParameters {
    std::uint32_t anchor_length = 50;
    std::uint32_t max_bundle_align_bp = 50000;
    std::uint64_t max_fallback_cells = 25000000;
};

struct BaseAlignmentResult {
    std::uint64_t bundle_count = 0;
    std::uint64_t aligned_bundle_count = 0;
    std::uint64_t skipped_bundle_count = 0;
};

BaseAlignmentResult align_chain_bundles(
    const std::filesystem::path& chains,
    const std::filesystem::path& fasta,
    const std::filesystem::path& fasta_index,
    const std::filesystem::path& paf_output,
    const std::filesystem::path& metadata_output,
    const BaseAlignmentParameters& parameters);

}  // namespace vallescope2
