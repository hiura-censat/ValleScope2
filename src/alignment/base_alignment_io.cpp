#include "base_alignment_internal.hpp"

#include <fstream>
#include <iterator>
#include <stdexcept>

namespace vallescope2 {
namespace alignment_detail {

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

double parse_f64(const std::string& value, const std::string& column) {
    try {
        std::size_t parsed = 0;
        const auto number = std::stod(value, &parsed);
        if (value.empty() || parsed != value.size()) throw std::out_of_range("bad");
        return number;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid number in " + column + ": " + value);
    }
}

struct AnchorPosition {
    std::string sequence_id;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
};

std::unordered_map<std::string, AnchorPosition> load_anchor_positions(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open grouped anchors: " + path.string());
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("empty grouped anchors file");
    }
    const auto columns = header_index(split_tab(line));
    const auto anchor_column = require_column(columns, "anchor_id");
    const auto sequence_column = require_column(columns, "sequence_id");
    const auto start_column = require_column(columns, "start");
    const auto end_column = require_column(columns, "end");

    std::unordered_map<std::string, AnchorPosition> positions;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= end_column) {
            throw std::runtime_error("invalid grouped anchor row");
        }
        AnchorPosition position;
        position.sequence_id = fields[sequence_column];
        position.start = parse_u64(fields[start_column], "start");
        position.end = parse_u64(fields[end_column], "end");
        positions.emplace(fields[anchor_column], std::move(position));
    }
    return positions;
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
    const auto score_column = require_column(columns, "chain_score");
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
        chain.score = parse_f64(fields[score_column], "chain_score");
        chain.ref_start = parse_u64(fields[ref_start_column], "ref_start");
        chain.ref_end = parse_u64(fields[ref_end_column], "ref_end");
        chain.query_start = parse_u64(fields[query_start_column], "query_start");
        chain.query_end = parse_u64(fields[query_end_column], "query_end");
        chains.push_back(std::move(chain));
    }
    return chains;
}

void attach_chain_anchors(
    std::vector<ChainBundle>& chains,
    const std::filesystem::path& path,
    const std::unordered_map<std::string, AnchorPosition>& positions,
    AnchorStore& store) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open chain anchors: " + path.string());
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("empty chain anchors file");
    }
    const auto columns = header_index(split_tab(line));
    const auto chain_column = require_column(columns, "chain_id");
    const auto anchor_a_column = require_column(columns, "anchor_a");
    const auto anchor_b_column = require_column(columns, "anchor_b");

    std::unordered_map<std::uint64_t, ChainBundle*> chain_by_id;
    for (auto& chain : chains) chain_by_id.emplace(chain.chain_id, &chain);

    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= anchor_b_column) {
            throw std::runtime_error("invalid chain anchor row");
        }
        const auto chain_id = parse_u64(fields[chain_column], "chain_id");
        const auto chain_found = chain_by_id.find(chain_id);
        if (chain_found == chain_by_id.end()) continue;

        const auto anchor_a = fields[anchor_a_column];
        const auto anchor_b = fields[anchor_b_column];
        const auto position_a = positions.find(anchor_a);
        const auto position_b = positions.find(anchor_b);
        if (position_a == positions.end() || position_b == positions.end()) {
            throw std::runtime_error("chain anchor is absent from grouped anchors");
        }

        AnchorPair anchor;
        anchor.anchor_a = anchor_a;
        anchor.anchor_b = anchor_b;
        anchor.ref_start = position_a->second.start;
        anchor.ref_end = position_a->second.end;
        anchor.query_start = position_b->second.start;
        anchor.query_end = position_b->second.end;
        const auto index = store.anchors.size();
        store.anchors.push_back(std::move(anchor));
        chain_found->second->anchor_indices.push_back(index);
    }
}

std::vector<ChainBundle> load_chain_files(
    const std::vector<std::filesystem::path>& paths,
    const std::vector<std::filesystem::path>& anchor_paths,
    const std::filesystem::path& grouped_anchors,
    AnchorStore& store) {
    if (paths.size() != anchor_paths.size()) {
        throw std::runtime_error("chain file and chain anchor file counts differ");
    }
    const auto positions = load_anchor_positions(grouped_anchors);
    std::vector<ChainBundle> bundles;
    for (std::size_t i = 0; i < paths.size(); ++i) {
        const auto& path = paths[i];
        const std::string filename = path.filename().string();
        const std::string source =
            filename.find("refined") == std::string::npos ? "primary" : "refined";
        auto loaded = load_chains(path, source);
        attach_chain_anchors(loaded, anchor_paths[i], positions, store);
        bundles.insert(bundles.end(),
                       std::make_move_iterator(loaded.begin()),
                       std::make_move_iterator(loaded.end()));
    }
    return bundles;
}

}  // namespace alignment_detail
}  // namespace vallescope2
