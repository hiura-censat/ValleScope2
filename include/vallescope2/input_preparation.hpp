#pragma once

#include "vallescope2/sequence_catalog.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace vallescope2 {

enum class ComparisonMode {
    self,
    target_query,
    multi_self,
};

struct PreparedFasta {
    std::filesystem::path fasta;
    std::filesystem::path index;
    std::filesystem::path gzi;
};

struct PreparationResult {
    ComparisonMode mode = ComparisonMode::self;
    PreparedFasta target;
    PreparedFasta query;
    std::filesystem::path sequence_table;
    SequenceCatalog catalog;
    std::uint64_t sequence_count = 0;
    std::uint64_t total_bases = 0;
};

PreparationResult prepare_inputs(
    const std::vector<std::filesystem::path>& input_paths,
    const std::filesystem::path& output_directory,
    bool combine);

}  // namespace vallescope2
