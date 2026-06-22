#pragma once

#include <cstdint>
#include <filesystem>

namespace vallescope2 {

struct StructuralTokenResult {
    std::uint64_t anchor_token_count = 0;
    std::uint64_t distance_token_count = 0;
    std::uint64_t sequence_count = 0;
    std::uint64_t skipped_anchor_count = 0;

    std::uint64_t token_count() const {
        return anchor_token_count + distance_token_count;
    }
};

StructuralTokenResult build_structural_tokens(
    const std::filesystem::path& grouped_anchors,
    std::uint32_t distance_bin_size,
    const std::filesystem::path& token_output,
    const std::filesystem::path& metadata_output);

}  // namespace vallescope2
