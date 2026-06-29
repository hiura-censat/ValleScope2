#include "vallescope2/alignment/base_alignment.hpp"

#include "vallescope2/context/anchor_grouping.hpp"

#include "bindings/cpp/WFAligner.hpp"

#include <htslib/faidx.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
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
    std::string source;
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
    std::uint64_t patch_count = 0;
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

std::vector<ChainBundle> load_chains(const std::filesystem::path& path,
                                     const std::string& source) {
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
        chain.source = source;
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

std::vector<ChainBundle> load_chain_files(
    const std::vector<std::filesystem::path>& paths) {
    std::vector<ChainBundle> bundles;
    for (const auto& path : paths) {
        const std::string filename = path.filename().string();
        const std::string source =
            filename.find("refined") == std::string::npos ? "primary" : "refined";
        auto loaded = load_chains(path, source);
        bundles.insert(bundles.end(),
                       std::make_move_iterator(loaded.begin()),
                       std::make_move_iterator(loaded.end()));
    }
    return bundles;
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

std::uint32_t levenshtein_distance(const std::string& a,
                                   const std::string& b) {
    std::vector<std::uint32_t> previous(b.size() + 1);
    std::vector<std::uint32_t> current(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); ++j) {
        previous[j] = static_cast<std::uint32_t>(j);
    }
    for (std::size_t i = 1; i <= a.size(); ++i) {
        current[0] = static_cast<std::uint32_t>(i);
        for (std::size_t j = 1; j <= b.size(); ++j) {
            const auto substitution =
                previous[j - 1] + (a[i - 1] == b[j - 1] ? 0U : 1U);
            const auto insertion = current[j - 1] + 1U;
            const auto deletion = previous[j] + 1U;
            current[j] = std::min({substitution, insertion, deletion});
        }
        previous.swap(current);
    }
    return previous[b.size()];
}

double sequence_identity(const std::string& ref,
                         const std::string& query) {
    const auto denominator = std::max(ref.size(), query.size());
    if (denominator == 0) return 1.0;
    const auto distance = levenshtein_distance(ref, query);
    return 1.0 - static_cast<double>(distance) /
                     static_cast<double>(denominator);
}

std::uint64_t positive_gap(const std::uint64_t left_end,
                           const std::uint64_t right_start) {
    return right_start > left_end ? right_start - left_end : 0;
}

bool same_bundle_track(const ChainBundle& left,
                       const ChainBundle& right) {
    return left.sample_a == right.sample_a &&
           left.sample_b == right.sample_b &&
           left.sequence_a == right.sequence_a &&
           left.sequence_b == right.sequence_b &&
           left.strand == right.strand;
}

std::uint64_t query_gap_between(const ChainBundle& left,
                                const ChainBundle& right) {
    if (left.strand == '+') return positive_gap(left.query_end, right.query_start);
    return right.query_end < left.query_start ? left.query_start - right.query_end : 0;
}

std::string fetch_patch_query(faidx_t* index,
                              const std::string& sequence_id,
                              const Interval& interval,
                              const char strand) {
    auto sequence = fetch_interval(index, sequence_id, interval);
    if (strand == '-') sequence = reverse_complement(sequence);
    return sequence;
}

bool extension_window_matches(faidx_t* index,
                              const ChainBundle& bundle,
                              const Interval& ref_interval,
                              const Interval& query_interval,
                              const BaseAlignmentParameters& parameters) {
    const auto ref = fetch_interval(index, bundle.sequence_a, ref_interval);
    const auto query = fetch_patch_query(
        index, bundle.sequence_b, query_interval, bundle.strand);
    return sequence_identity(ref, query) >= parameters.min_patch_identity;
}

bool can_patch_gap(faidx_t* index,
                   const ChainBundle& left,
                   const ChainBundle& right,
                   const BaseAlignmentParameters& parameters) {
    if (!same_bundle_track(left, right)) return false;
    if (right.ref_start < left.ref_end) return false;

    const auto ref_gap = right.ref_start - left.ref_end;
    const auto query_gap = query_gap_between(left, right);
    if (ref_gap > parameters.max_patch_gap_bp ||
        query_gap > parameters.max_patch_gap_bp) {
        return false;
    }
    if (ref_gap == 0 && query_gap == 0) return true;

    std::uint64_t left_ref_extension = 0;
    std::uint64_t left_query_extension = 0;
    std::uint64_t right_ref_extension = 0;
    std::uint64_t right_query_extension = 0;

    auto gap_is_closed = [&]() {
        return left_ref_extension + right_ref_extension >= ref_gap ||
               left_query_extension + right_query_extension >= query_gap;
    };

    while (!gap_is_closed()) {
        const auto remaining_ref =
            ref_gap - std::min(ref_gap, left_ref_extension + right_ref_extension);
        const auto remaining_query =
            query_gap - std::min(query_gap,
                                 left_query_extension + right_query_extension);
        if (remaining_ref == 0 || remaining_query == 0) return true;

        const auto left_ref_window =
            std::min<std::uint64_t>(parameters.patch_window_bp, remaining_ref);
        const auto left_query_window =
            std::min<std::uint64_t>(parameters.patch_window_bp, remaining_query);
        Interval left_ref{left.ref_end + left_ref_extension,
                          left.ref_end + left_ref_extension + left_ref_window};
        Interval left_query;
        if (left.strand == '+') {
            left_query = {left.query_end + left_query_extension,
                          left.query_end + left_query_extension + left_query_window};
        } else {
            left_query = {left.query_start - left_query_extension -
                              left_query_window,
                          left.query_start - left_query_extension};
        }
        if (!extension_window_matches(index, left, left_ref, left_query,
                                      parameters)) {
            return false;
        }
        left_ref_extension += left_ref_window;
        left_query_extension += left_query_window;
        if (gap_is_closed()) return true;

        const auto after_left_ref =
            ref_gap - std::min(ref_gap, left_ref_extension + right_ref_extension);
        const auto after_left_query =
            query_gap - std::min(query_gap,
                                 left_query_extension + right_query_extension);
        const auto right_ref_window =
            std::min<std::uint64_t>(parameters.patch_window_bp, after_left_ref);
        const auto right_query_window =
            std::min<std::uint64_t>(parameters.patch_window_bp, after_left_query);
        Interval right_ref{right.ref_start - right_ref_extension -
                               right_ref_window,
                           right.ref_start - right_ref_extension};
        Interval right_query;
        if (right.strand == '+') {
            right_query = {right.query_start - right_query_extension -
                               right_query_window,
                           right.query_start - right_query_extension};
        } else {
            right_query = {right.query_end + right_query_extension,
                           right.query_end + right_query_extension +
                               right_query_window};
        }
        if (!extension_window_matches(index, right, right_ref, right_query,
                                      parameters)) {
            return false;
        }
        right_ref_extension += right_ref_window;
        right_query_extension += right_query_window;
    }
    return true;
}

std::vector<ChainBundle> patch_adjacent_bundles(
    std::vector<ChainBundle> bundles,
    faidx_t* index,
    const BaseAlignmentParameters& parameters,
    std::uint64_t& patch_count) {
    patch_count = 0;
    std::sort(bundles.begin(), bundles.end(),
              [](const ChainBundle& left, const ChainBundle& right) {
                  return std::tie(left.sample_a, left.sample_b, left.sequence_a,
                                  left.sequence_b, left.strand, left.ref_start,
                                  left.ref_end, left.query_start,
                                  left.query_end, left.source, left.chain_id) <
                         std::tie(right.sample_a, right.sample_b, right.sequence_a,
                                  right.sequence_b, right.strand, right.ref_start,
                                  right.ref_end, right.query_start,
                                  right.query_end, right.source, right.chain_id);
              });

    std::vector<ChainBundle> patched;
    for (const auto& bundle : bundles) {
        if (patched.empty() ||
            !can_patch_gap(index, patched.back(), bundle, parameters)) {
            patched.push_back(bundle);
            continue;
        }
        auto& merged = patched.back();
        merged.ref_start = std::min(merged.ref_start, bundle.ref_start);
        merged.ref_end = std::max(merged.ref_end, bundle.ref_end);
        merged.query_start = std::min(merged.query_start, bundle.query_start);
        merged.query_end = std::max(merged.query_end, bundle.query_end);
        merged.patch_count += bundle.patch_count + 1;
        merged.source = "patched";
        ++patch_count;
    }
    return patched;
}

void write_metadata(const std::filesystem::path& path,
                    const BaseAlignmentParameters& parameters,
                    const BaseAlignmentResult& result,
                    const std::vector<std::filesystem::path>& chain_files,
                    const std::filesystem::path& fasta) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create base alignment metadata");
    output << "{\n"
           << "  \"input_chains\": [";
    for (std::size_t i = 0; i < chain_files.size(); ++i) {
        if (i) output << ", ";
        output << "\"" << chain_files[i].string() << "\"";
    }
    output << "],\n"
           << "  \"fasta\": \"" << fasta.string() << "\",\n"
           << "  \"alignment_mode\": \"global\",\n"
           << "  \"aligner\": \"WFA2-lib\",\n"
           << "  \"anchor_length\": " << parameters.anchor_length << ",\n"
           << "  \"max_bundle_align_bp\": " << parameters.max_bundle_align_bp << ",\n"
           << "  \"max_fallback_cells\": " << parameters.max_fallback_cells << ",\n"
           << "  \"max_patch_gap_bp\": " << parameters.max_patch_gap_bp << ",\n"
           << "  \"patch_window_bp\": " << parameters.patch_window_bp << ",\n"
           << "  \"min_patch_identity\": " << parameters.min_patch_identity << ",\n"
           << "  \"bundle_count\": " << result.bundle_count << ",\n"
           << "  \"patched_bundle_count\": " << result.patched_bundle_count << ",\n"
           << "  \"patch_count\": " << result.patch_count << ",\n"
           << "  \"aligned_bundle_count\": " << result.aligned_bundle_count << ",\n"
           << "  \"skipped_bundle_count\": " << result.skipped_bundle_count << "\n"
           << "}\n";
}

}  // namespace

BaseAlignmentResult align_chain_bundles(
    const std::vector<std::filesystem::path>& chain_files,
    const std::filesystem::path& fasta,
    const std::filesystem::path& fasta_index,
    const std::filesystem::path& paf_output,
    const std::filesystem::path& metadata_output,
    const BaseAlignmentParameters& parameters) {
    const auto loaded_bundles = load_chain_files(chain_files);
    BaseAlignmentResult result;

    std::unique_ptr<faidx_t, FaiDeleter> index(
        fai_load3(fasta.c_str(), fasta_index.c_str(), nullptr, 0));
    if (!index) throw std::runtime_error("cannot load FASTA index for base alignment");

    auto bundles = patch_adjacent_bundles(
        loaded_bundles, index.get(), parameters, result.patch_count);
    result.bundle_count = bundles.size();
    result.patched_bundle_count =
        loaded_bundles.size() > bundles.size() ? loaded_bundles.size() - bundles.size()
                                               : 0;

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
            << "\tbt:Z:" << bundle.source
            << "\tpg:i:" << bundle.patch_count
            << "\tas:i:" << alignment.score
            << "\tsa:Z:" << bundle.sample_a
            << "\tsb:Z:" << bundle.sample_b
            << '\n';
        ++result.aligned_bundle_count;
    }
    write_metadata(metadata_output, parameters, result, chain_files, fasta);
    return result;
}

}  // namespace vallescope2
