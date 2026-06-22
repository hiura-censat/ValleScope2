#pragma once

#include <cstdint>
#include <filesystem>

namespace vallescope2 {

struct StructuralContextResult {
    std::uint64_t anchor_count = 0;
    std::uint64_t context_group_count = 0;
    std::uint64_t tmus_count = 0;
};

StructuralContextResult build_structural_contexts(
    const std::filesystem::path& structural_tokens,
    const std::filesystem::path& sequence_table,
    std::uint32_t context_radius_tokens,
    const std::filesystem::path& anchor_context_output,
    const std::filesystem::path& context_group_output,
    const std::filesystem::path& tmus_output,
    const std::filesystem::path& metadata_output);

}  // namespace vallescope2
