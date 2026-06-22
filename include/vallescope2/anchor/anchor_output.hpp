#pragma once

#include "vallescope2/anchor/anchor_selection.hpp"
#include "vallescope2/anchor/genmap_runner.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace vallescope2 {

struct AnchorMetadata {
    std::uint32_t anchor_length = 50;
    std::uint32_t min_center_distance = 100;
    std::uint32_t target_density = 6000;
    std::uint64_t sequence_length = 0;
    std::uint32_t genmap_k = 6;
    std::uint32_t genmap_e = 0;
    std::string genmap_command;
    std::string genmap_version;
    std::vector<std::filesystem::path> input_fastas;
    std::filesystem::path genmap_input;
};

std::string sha256_file(const std::filesystem::path& path);

void write_anchor_bed(const std::filesystem::path& path,
                      const std::vector<AnchorCandidate>& anchors,
                      const SequenceCatalog& catalog);

void write_anchor_metadata(const std::filesystem::path& path,
                           const AnchorMetadata& metadata,
                           const AnchorSelectionResult& selection);

}  // namespace vallescope2
