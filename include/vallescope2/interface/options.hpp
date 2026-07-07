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
    bool dump_window_scores = false;
    std::uint32_t anchor_length = 50;
    std::uint32_t max_bundle_align_bp = 50000;
    std::uint32_t max_patch_gap_bp = 70000;
    std::uint32_t patch_window_bp = 1000;
    std::uint32_t patch_window_slack_bp = 300;
    std::uint32_t max_patch_indel_bp = 200;
    double min_patch_identity = 0.85;
    bool joint_patch_align = true;
    std::uint32_t joint_patch_flank_bp = 500;
    std::uint32_t max_joint_patch_bp = 80000;
    std::uint32_t max_wfa_memory_gb = 64;
    std::uint32_t max_chain_extension_bp = 70000;
    std::uint32_t min_chain_extension_anchors = 1;
    double min_chain_extension_score = 0.0;
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
