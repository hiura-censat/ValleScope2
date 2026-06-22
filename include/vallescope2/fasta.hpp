#pragma once

#include "vallescope2/sequence_catalog.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace vallescope2 {

std::string sample_name_from_path(const std::filesystem::path& path);

SequenceCatalog combine_fastas(
    const std::vector<std::filesystem::path>& input_paths,
    const std::vector<std::string>& roles,
    const std::filesystem::path& output_path);

}  // namespace vallescope2
