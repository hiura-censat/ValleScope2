#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace vallescope2 {

struct SequenceRecord {
    std::string sequence_id;
    std::string sample;
    std::string haplotype = "1";
    std::string original_id;
    std::filesystem::path source_fasta;
    std::uint64_t length = 0;
    std::string role;
};

class SequenceCatalog {
public:
    void add(SequenceRecord record);
    const SequenceRecord& at(const std::string& sequence_id) const;
    const std::vector<SequenceRecord>& records() const { return records_; }
    std::uint64_t total_bases() const { return total_bases_; }
    void write_tsv(const std::filesystem::path& path) const;

private:
    std::vector<SequenceRecord> records_;
    std::unordered_map<std::string, std::size_t> by_id_;
    std::uint64_t total_bases_ = 0;
};

}  // namespace vallescope2
