#include "assignment_internal.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace vallescope2 {
namespace {

std::unordered_map<std::string, GroupedAnchor>
load_grouped_anchors(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open grouped anchors: " +
                                 path.string());
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("empty grouped anchors file");
    }
    const auto columns = header_index(split_tab(line));
    const auto anchor_id_column = require_column(columns, "anchor_id");
    const auto sequence_column = require_column(columns, "sequence_id");
    const auto start_column = require_column(columns, "start");
    const auto end_column = require_column(columns, "end");
    const auto group_column = require_column(columns, "anchor_group_id");
    const auto groupable_column = require_column(columns, "groupable");

    std::unordered_map<std::string, GroupedAnchor> anchors;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= groupable_column) {
            throw std::runtime_error("invalid grouped anchor row");
        }
        GroupedAnchor anchor;
        anchor.sequence_id = fields[sequence_column];
        anchor.start = parse_u64(fields[start_column], "start");
        anchor.end = parse_u64(fields[end_column], "end");
        anchor.anchor_group_id = fields[group_column];
        anchor.groupable = fields[groupable_column] == "true";
        anchors.emplace(fields[anchor_id_column], std::move(anchor));
    }
    return anchors;
}

std::uint32_t intern_token(DataSet& data, const std::string& token) {
    const auto found = data.token_ids.find(token);
    if (found != data.token_ids.end()) return found->second;
    const auto id = static_cast<std::uint32_t>(data.token_ids.size());
    data.token_ids.emplace(token, id);
    return id;
}

std::uint32_t intern_tmus_key(DataSet& data, const std::string& key) {
    const auto found = data.tmus_key_ids.find(key);
    if (found != data.tmus_key_ids.end()) return found->second;
    const auto id = static_cast<std::uint32_t>(data.tmus_key_ids.size());
    data.tmus_key_ids.emplace(key, id);
    return id;
}

std::vector<std::uint32_t> parse_tokens(DataSet& data,
                                        const std::string& text) {
    std::vector<std::uint32_t> tokens;
    std::size_t begin = 0;
    while (begin <= text.size()) {
        const std::size_t separator = text.find('|', begin);
        const std::string token =
            separator == std::string::npos
                ? text.substr(begin)
                : text.substr(begin, separator - begin);
        if (!token.empty()) tokens.push_back(intern_token(data, token));
        if (separator == std::string::npos) break;
        begin = separator + 1;
    }
    return tokens;
}

std::vector<std::uint32_t> parse_tmus_keys(DataSet& data,
                                           const std::string& text) {
    std::vector<std::uint32_t> keys;
    if (text.empty() || text == ".") return keys;
    std::size_t begin = 0;
    while (begin <= text.size()) {
        const std::size_t separator = text.find(',', begin);
        const std::string key_text =
            separator == std::string::npos
                ? text.substr(begin)
                : text.substr(begin, separator - begin);
        if (!key_text.empty()) keys.push_back(intern_tmus_key(data, key_text));
        if (separator == std::string::npos) break;
        begin = separator + 1;
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

void load_sequence_table(const std::filesystem::path& path, DataSet& data) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open sequence table: " +
                                 path.string());
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("empty sequence table");
    }
    const auto columns = header_index(split_tab(line));
    const auto sample_column = require_column(columns, "sample");
    const auto role_column = require_column(columns, "role");
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= role_column) {
            throw std::runtime_error("invalid sequence table row");
        }
        data.sample_roles.emplace(fields[sample_column], fields[role_column]);
    }
}

std::string strip_alpha_suffix(const std::string& header) {
    constexpr const char* suffix = "_alpha_prime";
    constexpr std::size_t suffix_size = 12;
    if (header.size() < suffix_size ||
        header.compare(header.size() - suffix_size, suffix_size, suffix) != 0) {
        throw std::runtime_error("invalid alpha prime column: " + header);
    }
    return header.substr(0, header.size() - suffix_size);
}

}  // namespace

DataSet load_assignment_data(const std::filesystem::path& grouped_anchors,
                             const std::filesystem::path& anchor_contexts,
                             const std::filesystem::path& sequence_table) {
    DataSet data;
    load_sequence_table(sequence_table, data);
    const auto grouped = load_grouped_anchors(grouped_anchors);

    std::ifstream input(anchor_contexts);
    if (!input) {
        throw std::runtime_error("cannot open anchor contexts: " +
                                 anchor_contexts.string());
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("empty anchor contexts file");
    }
    const auto header = split_tab(line);
    const auto columns = header_index(header);
    const auto anchor_id_column = require_column(columns, "anchor_id");
    const auto sequence_column = require_column(columns, "sequence_id");
    const auto sample_column = require_column(columns, "sample");
    const auto context_key_column =
        require_column(columns, "canonical_context_key");
    const auto forward_column = require_column(columns, "forward_context_tokens");
    const auto tmus_keys_found = columns.find("tmus_keys");

    std::vector<std::size_t> alpha_columns;
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i].size() >= 12 &&
            header[i].compare(header[i].size() - 12, 12,
                              "_alpha_prime") == 0) {
            const auto sample = strip_alpha_suffix(header[i]);
            data.sample_ids.emplace(
                sample, static_cast<std::uint32_t>(data.samples.size()));
            data.samples.push_back(sample);
            alpha_columns.push_back(i);
        }
    }
    if (data.samples.empty()) {
        throw std::runtime_error("anchor contexts have no alpha prime columns");
    }
    data.anchors_by_sample.resize(data.samples.size());

    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tab(line);
        if (fields.size() <= forward_column) {
            throw std::runtime_error("invalid anchor context row");
        }
        const auto grouped_found = grouped.find(fields[anchor_id_column]);
        if (grouped_found == grouped.end()) {
            throw std::runtime_error("anchor context has no grouped anchor: " +
                                     fields[anchor_id_column]);
        }
        if (!grouped_found->second.groupable ||
            grouped_found->second.anchor_group_id == ".") {
            continue;
        }
        const auto sample_found = data.sample_ids.find(fields[sample_column]);
        if (sample_found == data.sample_ids.end()) {
            throw std::runtime_error(
                "anchor context sample has no alpha column: " +
                fields[sample_column]);
        }

        Anchor anchor;
        anchor.anchor_id = fields[anchor_id_column];
        anchor.sample = fields[sample_column];
        anchor.sequence_id = fields[sequence_column];
        anchor.start = grouped_found->second.start;
        anchor.end = grouped_found->second.end;
        anchor.anchor_group_id = grouped_found->second.anchor_group_id;
        anchor.context_key = fields[context_key_column];
        anchor.forward_tokens = parse_tokens(data, fields[forward_column]);
        if (tmus_keys_found != columns.end()) {
            if (fields.size() <= tmus_keys_found->second) {
                throw std::runtime_error("missing tMUS key value");
            }
            anchor.tmus_keys =
                parse_tmus_keys(data, fields[tmus_keys_found->second]);
        }
        anchor.alpha_prime.reserve(alpha_columns.size());
        for (const auto column : alpha_columns) {
            if (fields.size() <= column) {
                throw std::runtime_error("missing alpha prime value");
            }
            anchor.alpha_prime.push_back(parse_u64(fields[column], "alpha_prime"));
        }
        data.anchors.push_back(std::move(anchor));
        data.anchors_by_sample[sample_found->second].push_back(
            &data.anchors.back());
    }

    // Pointers into data.anchors may have been invalidated while the vector grew.
    for (auto& sample_anchors : data.anchors_by_sample) sample_anchors.clear();
    for (const auto& anchor : data.anchors) {
        const auto sample = data.sample_ids.at(anchor.sample);
        data.anchors_by_sample[sample].push_back(&anchor);
    }
    return data;
}

std::vector<Pair> make_ordered_pairs(const DataSet& data) {
    std::vector<std::uint32_t> targets;
    std::vector<std::uint32_t> queries;
    for (std::uint32_t i = 0; i < data.samples.size(); ++i) {
        const auto role_found = data.sample_roles.find(data.samples[i]);
        const auto role =
            role_found == data.sample_roles.end() ? "self" : role_found->second;
        if (role == "target") {
            targets.push_back(i);
        } else if (role == "query") {
            queries.push_back(i);
        }
    }
    std::vector<Pair> pairs;
    if (!targets.empty() || !queries.empty()) {
        for (const auto target : targets) {
            for (const auto query : queries) pairs.push_back({target, query});
        }
        return pairs;
    }
    if (data.samples.size() == 1) {
        pairs.push_back({0, 0});
        return pairs;
    }
    for (std::uint32_t target = 0; target < data.samples.size(); ++target) {
        for (std::uint32_t query = 0; query < data.samples.size(); ++query) {
            if (target != query) pairs.push_back({target, query});
        }
    }
    return pairs;
}

}  // namespace vallescope2
