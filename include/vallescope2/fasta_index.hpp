#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace vallescope2 {

struct IndexedSequence {
    std::string name;
    std::uint64_t length = 0;
};

void build_fasta_index(const std::filesystem::path& fasta_path,
                       const std::filesystem::path& index_path,
                       const std::filesystem::path& gzi_path = {});

std::vector<IndexedSequence> read_fasta_index(
    const std::filesystem::path& fasta_path,
    const std::filesystem::path& index_path,
    const std::filesystem::path& gzi_path = {});

}  // namespace vallescope2
