#include "vallescope2/alignment/base_alignment.hpp"

#include "vallescope2/context/anchor_grouping.hpp"

#include "bindings/cpp/WFAligner.hpp"

#include <htslib/faidx.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace vallescope2 {
namespace {

struct FaiDeleter {
    void operator()(faidx_t* index) const { fai_destroy(index); }
};

struct SequenceDeleter {
    void operator()(char* sequence) const { free(sequence); }
};

struct ChainBundle {
    std::uint64_t chain_id = 0;
    std::string sample_a;
    std::string sample_b;
    std::string sequence_a;
    std::string sequence_b;
    char strand = '+';
    std::uint64_t ref_start = 0;
    std::uint64_t ref_end = 0;
    std::uint64_t query_start = 0;
    std::uint64_t query_end = 0;
};

struct Interval {
    std::uint64_t start = 0;
    std::uint64_t end = 0;
};

struct Alignment {
    std::string cigar;
    std::uint64_t matches = 0;
    std::uint64_t block_length = 0;
    std::int64_t score = 0;
};

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            fields.push_back(line.substr(begin));
            break;
        }
        fields.push_back(line.substr(begin, tab - begin));
        begin = tab + 1;
    }
    return fields;
}

std::unordered_map<std::string, std::size_t>
header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < header.size(); ++i) index.emplace(header[i], i);
    return index;
}

std::size_t require_column(
    const std::unordered_map<std::string, std::size_t>& columns,
    const std::string& name) {
    const auto found = columns.find(name);
    if (found == columns.end()) throw std::runtime_error("missing column: " + name);
    return found->second;
}

std::uint64_t parse_u64(const std::string& value, const std::string& column) {
    try {
        std::size_t parsed = 0;
        const auto number = std::stoull(value, &parsed);
        if (value.empty() || parsed != value.size()) throw std::out_of_range("bad");
        return number;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer in " + column + ": " + value);
    }
}

std::vector<ChainBundle> load_chains(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open chains: " + path.string());
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("empty chains file");
    const auto columns = header_index(split_tab(line));
    const auto chain_column = require_column(columns, "chain_id");
    const auto sample_a_column = require_column(columns, "sample_a");
    const auto sample_b_column = require_column(columns, "sample_b");
    const auto sequence_a_column = require_column(columns, "sequence_a");
    const auto sequence_b_column = require_column(columns, "sequence_b");
    const auto strand_column = require_column(columns, "assign_strand");
    const auto ref_start_column = require_column(columns, "ref_start");
    const auto ref_end_column = require_column(columns, "ref_end");
    const auto query_start_column = require_column(columns, "query_start");
    const auto query_end_column = require_column(columns, "query_end");

    std::vector<ChainBundle> chains;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= query_end_column)
            throw std::runtime_error("invalid chain row");
        ChainBundle chain;
        chain.chain_id = parse_u64(fields[chain_column], "chain_id");
        chain.sample_a = fields[sample_a_column];
        chain.sample_b = fields[sample_b_column];
        chain.sequence_a = fields[sequence_a_column];
        chain.sequence_b = fields[sequence_b_column];
        chain.strand = fields[strand_column] == "-" ? '-' : '+';
        chain.ref_start = parse_u64(fields[ref_start_column], "ref_start");
        chain.ref_end = parse_u64(fields[ref_end_column], "ref_end");
        chain.query_start = parse_u64(fields[query_start_column], "query_start");
        chain.query_end = parse_u64(fields[query_end_column], "query_end");
        chains.push_back(std::move(chain));
    }
    return chains;
}

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

Alignment global_align_wfa2(const std::string& ref,
                            const std::string& query) {
    wfa::WFAlignerGapAffine aligner(
        4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    const auto status = aligner.alignEnd2End(query, ref);
    if (status < 0) {
        throw std::runtime_error("WFA2 alignment failed");
    }
    return summarize_cigar(aligner.getCIGAR(true), ref, query,
                           aligner.getAlignmentScore());
}

void write_metadata(const std::filesystem::path& path,
                    const BaseAlignmentParameters& parameters,
                    const BaseAlignmentResult& result,
                    const std::filesystem::path& chains,
                    const std::filesystem::path& fasta) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create base alignment metadata");
    output << "{\n"
           << "  \"input_chains\": \"" << chains.string() << "\",\n"
           << "  \"fasta\": \"" << fasta.string() << "\",\n"
           << "  \"alignment_mode\": \"global\",\n"
           << "  \"aligner\": \"WFA2-lib\",\n"
           << "  \"anchor_length\": " << parameters.anchor_length << ",\n"
           << "  \"max_bundle_align_bp\": " << parameters.max_bundle_align_bp << ",\n"
           << "  \"max_fallback_cells\": " << parameters.max_fallback_cells << ",\n"
           << "  \"bundle_count\": " << result.bundle_count << ",\n"
           << "  \"aligned_bundle_count\": " << result.aligned_bundle_count << ",\n"
           << "  \"skipped_bundle_count\": " << result.skipped_bundle_count << "\n"
           << "}\n";
}

}  // namespace

BaseAlignmentResult align_chain_bundles(
    const std::filesystem::path& chains,
    const std::filesystem::path& fasta,
    const std::filesystem::path& fasta_index,
    const std::filesystem::path& paf_output,
    const std::filesystem::path& metadata_output,
    const BaseAlignmentParameters& parameters) {
    const auto bundles = load_chains(chains);
    BaseAlignmentResult result;
    result.bundle_count = bundles.size();

    std::unique_ptr<faidx_t, FaiDeleter> index(
        fai_load3(fasta.c_str(), fasta_index.c_str(), nullptr, 0));
    if (!index) throw std::runtime_error("cannot load FASTA index for base alignment");

    std::ofstream paf(paf_output);
    if (!paf) throw std::runtime_error("cannot create bundle alignment PAF");

    for (const auto& bundle : bundles) {
        const auto ref_length = sequence_length(index.get(), bundle.sequence_a);
        const auto query_length = sequence_length(index.get(), bundle.sequence_b);
        const auto ref_interval = expand_interval(
            bundle.ref_start, bundle.ref_end, ref_length, parameters.anchor_length);
        const auto query_interval = expand_interval(
            bundle.query_start, bundle.query_end, query_length, parameters.anchor_length);
        const auto ref_span = ref_interval.end - ref_interval.start;
        const auto query_span = query_interval.end - query_interval.start;
        if (ref_span > parameters.max_bundle_align_bp ||
            query_span > parameters.max_bundle_align_bp) {
            ++result.skipped_bundle_count;
            continue;
        }

        const auto ref = fetch_interval(index.get(), bundle.sequence_a, ref_interval);
        auto query = fetch_interval(index.get(), bundle.sequence_b, query_interval);
        if (bundle.strand == '-') query = reverse_complement(query);

        Alignment alignment;
        try {
            alignment = global_align_wfa2(ref, query);
        } catch (const std::exception&) {
            try {
                alignment = global_align_fallback(
                    ref, query, parameters.max_fallback_cells);
            } catch (const std::exception&) {
                ++result.skipped_bundle_count;
                continue;
            }
        }

        paf << bundle.sequence_b << '\t' << query_length << '\t'
            << query_interval.start << '\t' << query_interval.end << '\t'
            << bundle.strand << '\t' << bundle.sequence_a << '\t'
            << ref_length << '\t' << ref_interval.start << '\t'
            << ref_interval.end << '\t' << alignment.matches << '\t'
            << alignment.block_length << "\t60"
            << "\tcg:Z:" << alignment.cigar
            << "\tch:i:" << bundle.chain_id
            << "\tas:i:" << alignment.score
            << "\tsa:Z:" << bundle.sample_a
            << "\tsb:Z:" << bundle.sample_b
            << '\n';
        ++result.aligned_bundle_count;
    }
    write_metadata(metadata_output, parameters, result, chains, fasta);
    return result;
}

}  // namespace vallescope2
