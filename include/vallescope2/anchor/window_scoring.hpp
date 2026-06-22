#pragma once

#include "vallescope2/sequence_catalog.hpp"

#include <cstdint>
#include <filesystem>
#include <functional>

namespace vallescope2 {

struct WindowScoringParameters {
    std::uint32_t anchor_length = 50;
    std::uint32_t kmer_length = 6;
};

struct WindowScore {
    std::uint32_t sequence_id = 0;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    double score = 0.0;
};

struct WindowScoringResult {
    std::uint64_t window_count = 0;
    std::uint64_t gap_positions = 0;
    std::uint64_t saturated_positions = 0;
};

WindowScoringResult score_windows(
    const std::filesystem::path& frequency_raw,
    const SequenceCatalog& catalog,
    const WindowScoringParameters& parameters,
    const std::function<void(const WindowScore&)>& consume);

}  // namespace vallescope2
