#include "vallescope2/fasta_index.hpp"

#include <htslib/faidx.h>

#include <memory>
#include <stdexcept>

namespace vallescope2 {
namespace {

struct FaiDeleter {
    void operator()(faidx_t* index) const { fai_destroy(index); }
};

const char* optional_path(const std::filesystem::path& path) {
    return path.empty() ? nullptr : path.c_str();
}

}  // namespace

void build_fasta_index(const std::filesystem::path& fasta_path,
                       const std::filesystem::path& index_path,
                       const std::filesystem::path& gzi_path) {
    if (fai_build3(fasta_path.c_str(), index_path.c_str(),
                   optional_path(gzi_path)) != 0) {
        throw std::runtime_error("HTSlib failed to create FASTA index: " +
                                 index_path.string());
    }
}

std::vector<IndexedSequence> read_fasta_index(
    const std::filesystem::path& fasta_path,
    const std::filesystem::path& index_path,
    const std::filesystem::path& gzi_path) {
    std::unique_ptr<faidx_t, FaiDeleter> index(
        fai_load3(fasta_path.c_str(), index_path.c_str(),
                  optional_path(gzi_path), 0));
    if (!index) {
        throw std::runtime_error("HTSlib failed to load FASTA index: " +
                                 index_path.string());
    }

    std::vector<IndexedSequence> sequences;
    const int count = faidx_nseq(index.get());
    sequences.reserve(count);
    for (int i = 0; i < count; ++i) {
        const char* name = faidx_iseq(index.get(), i);
        const hts_pos_t length = faidx_seq_len64(index.get(), name);
        if (length < 0) {
            throw std::runtime_error("HTSlib failed to read sequence length");
        }
        sequences.push_back({name, static_cast<std::uint64_t>(length)});
    }
    return sequences;
}

}  // namespace vallescope2
