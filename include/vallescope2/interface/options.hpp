#pragma once

#include "vallescope2/anchor/genmap_runner.hpp"

#include <cstdint>
#include <filesystem>
#include <vector>

namespace vallescope2 {

struct ProgramOptions {
    bool combine = false;
    bool debug = false;
    bool dump_window_scores = false;
    std::uint32_t anchor_length = 50;
    std::uint32_t min_center_distance = 100;
    std::uint32_t target_density = 6000;
    std::uint32_t distance_bin_size = 50;
    GenmapParameters genmap;
    std::vector<std::filesystem::path> input_paths;
};

ProgramOptions parse_arguments(int argc, char* argv[]);

}  // namespace vallescope2
