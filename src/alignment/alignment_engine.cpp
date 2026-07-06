#include "base_alignment_internal.hpp"

#include "vallescope2/context/anchor_grouping.hpp"
#include "bindings/cpp/WFAligner.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <vector>

namespace vallescope2 {
namespace alignment_detail {

void FaiDeleter::operator()(faidx_t* index) const { fai_destroy(index); }

namespace {

struct SequenceDeleter {
    void operator()(char* sequence) const { free(sequence); }
};

}  // namespace

std::uint64_t sequence_length(faidx_t* index, const std::string& sequence_id) {
    const hts_pos_t length = faidx_seq_len64(index, sequence_id.c_str());
    if (length < 0) {
        throw std::runtime_error("sequence is absent from FASTA: " + sequence_id);
    }
    return static_cast<std::uint64_t>(length);
}

Interval expand_interval(const std::uint64_t start,
                         const std::uint64_t end,
                         const std::uint64_t sequence_length,
                         const std::uint32_t anchor_length) {
    const std::uint64_t flank = anchor_length / 2;
    Interval interval;
    interval.start = start > flank ? start - flank : 0;
    interval.end = std::min(sequence_length, end + flank);
    if (interval.start >= interval.end) {
        throw std::runtime_error("empty bundle interval after expansion");
    }
    return interval;
}

std::string fetch_interval(faidx_t* index,
                           const std::string& sequence_id,
                           const Interval& interval) {
    if (interval.start == interval.end) return {};
    if (interval.start > interval.end) {
        throw std::runtime_error("invalid FASTA interval: " + sequence_id);
    }
    hts_pos_t fetched_length = 0;
    std::unique_ptr<char, SequenceDeleter> fetched(
        faidx_fetch_seq64(index, sequence_id.c_str(),
                          static_cast<hts_pos_t>(interval.start),
                          static_cast<hts_pos_t>(interval.end - 1),
                          &fetched_length));
    if (!fetched ||
        fetched_length != static_cast<hts_pos_t>(interval.end - interval.start)) {
        throw std::runtime_error("failed to fetch FASTA interval: " + sequence_id);
    }
    std::string sequence(fetched.get(), static_cast<std::size_t>(fetched_length));
    std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                   [](const unsigned char value) {
                       return static_cast<char>(std::toupper(value));
                   });
    return sequence;
}

void append_cigar(std::string& cigar,
                  std::uint64_t& run_length,
                  char& run_op,
                  const char op) {
    if (run_length > 0 && run_op == op) {
        ++run_length;
        return;
    }
    if (run_length > 0) {
        cigar += std::to_string(run_length);
        cigar.push_back(run_op);
    }
    run_op = op;
    run_length = 1;
}

Alignment global_align_fallback(const std::string& ref,
                                const std::string& query,
                                const std::uint64_t max_cells) {
    const std::size_t rows = query.size() + 1;
    const std::size_t cols = ref.size() + 1;
    if (rows > 0 && cols > max_cells / rows) {
        throw std::runtime_error("fallback aligner cell limit exceeded");
    }

    constexpr std::int32_t match = 1;
    constexpr std::int32_t mismatch = -1;
    constexpr std::int32_t gap = -1;
    std::vector<std::int32_t> scores(rows * cols, 0);
    std::vector<char> trace(rows * cols, 0);
    auto at = [cols](const std::size_t row, const std::size_t col) {
        return row * cols + col;
    };
    for (std::size_t row = 1; row < rows; ++row) {
        scores[at(row, 0)] = static_cast<std::int32_t>(row) * gap;
        trace[at(row, 0)] = 'I';
    }
    for (std::size_t col = 1; col < cols; ++col) {
        scores[at(0, col)] = static_cast<std::int32_t>(col) * gap;
        trace[at(0, col)] = 'D';
    }
    for (std::size_t row = 1; row < rows; ++row) {
        for (std::size_t col = 1; col < cols; ++col) {
            const bool is_match = query[row - 1] == ref[col - 1];
            const auto diagonal =
                scores[at(row - 1, col - 1)] + (is_match ? match : mismatch);
            const auto insertion = scores[at(row - 1, col)] + gap;
            const auto deletion = scores[at(row, col - 1)] + gap;
            if (diagonal >= insertion && diagonal >= deletion) {
                scores[at(row, col)] = diagonal;
                trace[at(row, col)] = is_match ? '=' : 'X';
            } else if (insertion >= deletion) {
                scores[at(row, col)] = insertion;
                trace[at(row, col)] = 'I';
            } else {
                scores[at(row, col)] = deletion;
                trace[at(row, col)] = 'D';
            }
        }
    }

    std::size_t row = query.size();
    std::size_t col = ref.size();
    std::string reversed_cigar;
    std::uint64_t run_length = 0;
    char run_op = 0;
    Alignment alignment;
    alignment.score = scores[at(row, col)];
    while (row > 0 || col > 0) {
        const char op = trace[at(row, col)];
        append_cigar(reversed_cigar, run_length, run_op, op);
        if (op == '=') {
            ++alignment.matches;
            ++alignment.block_length;
            --row;
            --col;
        } else if (op == 'X') {
            ++alignment.block_length;
            --row;
            --col;
        } else if (op == 'I') {
            ++alignment.block_length;
            --row;
        } else {
            ++alignment.block_length;
            --col;
        }
    }
    if (run_length > 0) {
        reversed_cigar += std::to_string(run_length);
        reversed_cigar.push_back(run_op);
    }

    std::vector<std::pair<std::uint64_t, char>> ops;
    std::size_t begin = 0;
    while (begin < reversed_cigar.size()) {
        std::size_t end = begin;
        while (end < reversed_cigar.size() &&
               std::isdigit(static_cast<unsigned char>(reversed_cigar[end]))) {
            ++end;
        }
        ops.emplace_back(parse_u64(reversed_cigar.substr(begin, end - begin),
                                   "cigar_length"),
                         reversed_cigar[end]);
        begin = end + 1;
    }
    for (auto it = ops.rbegin(); it != ops.rend(); ++it) {
        alignment.cigar += std::to_string(it->first);
        alignment.cigar.push_back(it->second);
    }
    return alignment;
}

Alignment summarize_cigar(const std::string& cigar,
                          const std::string& ref,
                          const std::string& query,
                          const std::int64_t score) {
    Alignment alignment;
    alignment.cigar = cigar;
    alignment.score = score;
    std::size_t ref_offset = 0;
    std::size_t query_offset = 0;
    std::size_t begin = 0;
    while (begin < cigar.size()) {
        std::size_t end = begin;
        while (end < cigar.size() &&
               std::isdigit(static_cast<unsigned char>(cigar[end]))) {
            ++end;
        }
        if (end == begin || end >= cigar.size()) {
            throw std::runtime_error("invalid WFA2 CIGAR");
        }
        const auto length = parse_u64(cigar.substr(begin, end - begin),
                                      "cigar_length");
        const char op = cigar[end];
        if (op == '=' || op == 'X' || op == 'M') {
            alignment.block_length += length;
            for (std::uint64_t i = 0; i < length; ++i) {
                if (ref_offset + i < ref.size() &&
                    query_offset + i < query.size() &&
                    ref[ref_offset + i] == query[query_offset + i]) {
                    ++alignment.matches;
                }
            }
            ref_offset += length;
            query_offset += length;
        } else if (op == 'I') {
            alignment.block_length += length;
            query_offset += length;
        } else if (op == 'D') {
            alignment.block_length += length;
            ref_offset += length;
        } else {
            throw std::runtime_error("unsupported WFA2 CIGAR operation");
        }
        begin = end + 1;
    }
    return alignment;
}

std::uint64_t gibibytes_to_bytes(const std::uint32_t gibibytes) {
    return static_cast<std::uint64_t>(gibibytes) * 1024ULL * 1024ULL * 1024ULL;
}

Alignment global_align_wfa2(const std::string& ref,
                            const std::string& query,
                            const std::uint32_t max_memory_gb) {
    wfa::WFAlignerGapAffine aligner(
        4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    const auto max_memory = gibibytes_to_bytes(max_memory_gb);
    aligner.setMaxMemory(max_memory, max_memory);
    const auto status = aligner.alignEnd2End(ref, query);
    if (status < 0) {
        throw std::runtime_error("WFA2 alignment failed");
    }
    return summarize_cigar(aligner.getCIGAR(true), ref, query,
                           aligner.getAlignmentScore());
}

void append_cigar_operation(std::string& cigar,
                            const std::uint64_t length,
                            const char op) {
    if (length == 0) return;
    if (!cigar.empty() && cigar.back() == op) {
        std::size_t digit_end = cigar.size() - 1;
        std::size_t digit_begin = digit_end;
        while (digit_begin > 0 &&
               std::isdigit(static_cast<unsigned char>(cigar[digit_begin - 1]))) {
            --digit_begin;
        }
        const auto previous = parse_u64(
            cigar.substr(digit_begin, digit_end - digit_begin), "cigar_length");
        cigar.erase(digit_begin);
        cigar += std::to_string(previous + length);
        cigar.push_back(op);
        return;
    }
    cigar += std::to_string(length);
    cigar.push_back(op);
}

void append_cigar_string(std::string& destination,
                         const std::string& source) {
    std::size_t begin = 0;
    while (begin < source.size()) {
        std::size_t end = begin;
        while (end < source.size() &&
               std::isdigit(static_cast<unsigned char>(source[end]))) {
            ++end;
        }
        if (end == begin || end >= source.size()) {
            throw std::runtime_error("invalid segment CIGAR");
        }
        append_cigar_operation(
            destination,
            parse_u64(source.substr(begin, end - begin), "cigar_length"),
            source[end]);
        begin = end + 1;
    }
}

Alignment align_segment(const std::string& ref,
                        const std::string& query,
                        const BaseAlignmentParameters& parameters) {
    Alignment alignment;
    if (ref.empty() && query.empty()) return alignment;
    if (ref.empty()) {
        alignment.cigar = std::to_string(query.size()) + "I";
        alignment.block_length = query.size();
        alignment.score = -static_cast<std::int64_t>(query.size());
        return alignment;
    }
    if (query.empty()) {
        alignment.cigar = std::to_string(ref.size()) + "D";
        alignment.block_length = ref.size();
        alignment.score = -static_cast<std::int64_t>(ref.size());
        return alignment;
    }
    try {
        return global_align_wfa2(ref, query, parameters.max_wfa_memory_gb);
    } catch (const std::exception&) {
        return global_align_fallback(ref, query, parameters.max_fallback_cells);
    }
}

}  // namespace alignment_detail
}  // namespace vallescope2
