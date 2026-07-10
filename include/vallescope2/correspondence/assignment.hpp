#pragma once

#include <cstdint>
#include <filesystem>

namespace vallescope2 {

enum class PairMergeMode {
    reciprocal,
    union_mode
};

struct AssignmentParameters {
    std::uint32_t min_shared_tmus = 3;
    std::int64_t primary_margin = 5;
    std::int64_t min_candidate_score = 10;
    PairMergeMode pair_merge_mode = PairMergeMode::union_mode;
};

struct AssignmentResult {
    std::uint64_t ordered_pair_count = 0;
    std::uint64_t query_anchor_count = 0;
    std::uint64_t primary_count = 0;
    std::uint64_t ambiguous_query_count = 0;
    std::uint64_t ambiguous_row_count = 0;
    std::uint64_t unmatched_count = 0;
    std::uint64_t exact_context_primary_count = 0;
    std::uint64_t shared_tmus_primary_count = 0;
    std::uint64_t shared_tmus_tie_fallback_count = 0;
    std::uint64_t edit_distance_primary_count = 0;
    std::uint64_t correspondence_count = 0;
    std::uint64_t reciprocal_correspondence_count = 0;
    std::uint64_t forward_only_correspondence_count = 0;
    std::uint64_t reverse_only_correspondence_count = 0;
    std::uint64_t discordant_correspondence_count = 0;
};

AssignmentResult assign_context_correspondences(
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& anchor_contexts,
    const std::filesystem::path& sequence_table,
    const std::filesystem::path& assignment_output,
    const std::filesystem::path& correspondence_output,
    const std::filesystem::path& correspondence_metadata_output,
    const std::filesystem::path& metadata_output,
    const AssignmentParameters& parameters);

}  // namespace vallescope2
