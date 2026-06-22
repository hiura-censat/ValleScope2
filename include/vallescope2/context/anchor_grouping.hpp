#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

namespace vallescope2 {

struct AnchorGroupingResult {
    std::uint64_t anchor_count = 0;
    std::uint64_t grouped_anchor_count = 0;
    std::uint64_t ungrouped_anchor_count = 0;
    std::uint64_t group_count = 0;
};

std::string reverse_complement(const std::string& sequence);

AnchorGroupingResult group_anchors(
    const std::filesystem::path& anchor_bed,
    const std::filesystem::path& fasta,
    const std::filesystem::path& fasta_index,
    std::uint32_t anchor_length,
    const std::filesystem::path& grouped_anchors_output,
    const std::filesystem::path& anchor_groups_output);

}  // namespace vallescope2
