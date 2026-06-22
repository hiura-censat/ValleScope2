#include "vallescope2/sequence_catalog.hpp"

#include <fstream>
#include <stdexcept>

namespace vallescope2 {

void SequenceCatalog::add(SequenceRecord record) {
    if (record.sequence_id.empty()) {
        throw std::runtime_error("sequence ID must not be empty");
    }
    if (!by_id_.emplace(record.sequence_id, records_.size()).second) {
        throw std::runtime_error("duplicate sequence ID: " + record.sequence_id);
    }
    total_bases_ += record.length;
    records_.push_back(std::move(record));
}

const SequenceRecord& SequenceCatalog::at(
    const std::string& sequence_id) const {
    const auto found = by_id_.find(sequence_id);
    if (found == by_id_.end()) {
        throw std::runtime_error("unknown sequence ID: " + sequence_id);
    }
    return records_[found->second];
}

void SequenceCatalog::write_tsv(const std::filesystem::path& path) const {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot create sequence table: " + path.string());
    }
    output << "sequence_id\tsample\thaplotype\toriginal_id\tsource_file"
              "\tlength\trole\n";
    for (const auto& record : records_) {
        output << record.sequence_id << '\t' << record.sample << '\t'
               << record.haplotype << '\t' << record.original_id << '\t'
               << record.source_fasta.string() << '\t' << record.length << '\t'
               << record.role << '\n';
    }
}

}  // namespace vallescope2
