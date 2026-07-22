#pragma once

#include "vallescope2/anchor/genmap_runner.hpp"
#include "vallescope2/bundle/chaining.hpp"
#include "vallescope2/correspondence/assignment.hpp"

#include <cstdint>
#include <filesystem>
#include <vector>

namespace vallescope2 {

struct ProgramOptions {
    bool combine = false;
    bool debug = false;
    bool base_align = true;
    bool chain_extension = true;
    bool patch_zdrop = false;
    bool dump_window_scores = false;
    std::uint32_t anchor_length = 50;
    std::uint32_t max_bundle_align_bp = 1000000;
    std::uint32_t max_patch_gap_bp = 70000;
    std::uint32_t patch_flank_bp = 500;
    double min_patch_identity = 0.90;
    std::uint32_t min_patch_long_indel_bp = 100;
    std::uint32_t min_patch_rescue_flank_bp = 200;
    double min_patch_rescue_flank_identity = 0.95;
    std::uint32_t max_patch_rescue_extra_indel_bp = 100;
    std::uint32_t patch_quality_window_bp = 500;
    double min_patch_endpoint_identity = 0.85;
    double max_patch_short_error_density = 0.15;
    double min_patch_multi_segment_identity = 0.99;
    std::uint32_t max_patch_multi_short_indel_bp = 100;
    double max_patch_multi_short_error_density = 0.15;
    bool patch_multi_event_allow_deletions = true;
    std::uint32_t max_wfa_memory_gb = 64;
    std::uint32_t max_chain_extension_bp = 1000000;
    std::uint32_t min_chain_extension_anchors = 5;
    double min_chain_extension_score = 100.0;
    std::uint32_t min_copy_support_anchors = 1;
    std::uint32_t min_center_distance = 50;
    std::uint32_t target_density = 15000;
    std::uint32_t distance_bin_size = 10;
    std::uint32_t context_radius_tokens = 50;
    AssignmentParameters assignment;
    ChainingParameters chaining;
    GenmapParameters genmap;
    std::vector<std::filesystem::path> input_paths;
};

ProgramOptions parse_arguments(int argc, char* argv[]);

}  // namespace vallescope2
