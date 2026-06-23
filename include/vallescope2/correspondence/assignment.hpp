#pragma once

#include <cstdint>
#include <filesystem>

namespace vallescope2 {

struct AssignmentParameters {
    std::uint32_t beta_tolerance = 30;
    std::int64_t primary_margin = 5;
    std::int64_t min_candidate_score = -10;
};

struct AssignmentResult {
    std::uint64_t ordered_pair_count = 0;
    std::uint64_t query_anchor_count = 0;
    std::uint64_t primary_count = 0;
    std::uint64_t ambiguous_query_count = 0;
    std::uint64_t ambiguous_row_count = 0;
    std::uint64_t unmatched_count = 0;
    std::uint64_t exact_context_primary_count = 0;
    std::uint64_t edit_distance_primary_count = 0;
};

AssignmentResult assign_context_correspondences(
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& anchor_contexts,
    const std::filesystem::path& sequence_table,
    const std::filesystem::path& assignment_output,
    const std::filesystem::path& metadata_output,
    const AssignmentParameters& parameters);

}  // namespace vallescope2
