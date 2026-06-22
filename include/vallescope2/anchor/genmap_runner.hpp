#pragma once

#include "vallescope2/input_preparation.hpp"
#include "vallescope2/sequence_catalog.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace vallescope2 {

struct GenmapParameters {
    std::filesystem::path executable = "genmap";
    std::uint32_t kmer_length = 6;
    std::uint32_t errors = 0;
};

struct GenmapResult {
    std::filesystem::path input_fasta;
    std::filesystem::path index_directory;
    std::filesystem::path output_prefix;
    std::filesystem::path frequency_raw;
    std::filesystem::path log;
    std::string command;
    std::string version;
    SequenceCatalog catalog;
};

GenmapResult run_genmap(
    const std::vector<std::filesystem::path>& input_paths,
    const PreparationResult& prepared,
    const std::filesystem::path& work_directory,
    const GenmapParameters& parameters);

}  // namespace vallescope2
