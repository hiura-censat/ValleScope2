#pragma once

#include <cstdint>
#include <filesystem>

namespace vallescope2 {

struct ChainingParameters {
    std::uint32_t predecessor_count = 50;
    std::uint32_t max_chain_gap = 50000;
    double gap_weight = 1.0;
    std::uint32_t gap_unit = 10;
    std::uint32_t min_chain_anchors = 10;
    double min_chain_score = 0.0;
    std::uint32_t refinement_window = 50000;
    std::uint32_t refinement_min_chain_anchors = 5;
};

struct ChainingResult {
    std::uint64_t candidate_count = 0;
    std::uint64_t chain_count = 0;
    std::uint64_t chain_anchor_count = 0;
    std::uint64_t refined_chain_count = 0;
    std::uint64_t refined_chain_anchor_count = 0;
};

ChainingResult build_anchor_chains(
    const std::filesystem::path& correspondences,
    const std::filesystem::path& assignments,
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& chain_output,
    const std::filesystem::path& chain_anchor_output,
    const std::filesystem::path& metadata_output,
    const std::filesystem::path& refined_chain_output,
    const std::filesystem::path& refined_chain_anchor_output,
    const std::filesystem::path& refined_metadata_output,
    const ChainingParameters& parameters);

}  // namespace vallescope2
