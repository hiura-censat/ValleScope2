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
    bool dump_window_scores = false;
    std::uint32_t anchor_length = 50;
    std::uint32_t max_bundle_align_bp = 1000000;
    std::uint32_t max_patch_gap_bp = 70000;
    std::uint32_t patch_window_bp = 1000;
    double min_patch_identity = 0.85;
    std::uint32_t min_center_distance = 100;
    std::uint32_t target_density = 6000;
    std::uint32_t distance_bin_size = 10;
    std::uint32_t context_radius_tokens = 50;
    AssignmentParameters assignment;
    ChainingParameters chaining;
    GenmapParameters genmap;
    std::vector<std::filesystem::path> input_paths;
};

ProgramOptions parse_arguments(int argc, char* argv[]);

}  // namespace vallescope2
