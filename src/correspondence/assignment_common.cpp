#include "assignment_internal.hpp"

#include <stdexcept>

namespace vallescope2 {

std::string pair_merge_mode_name(const PairMergeMode mode) {
    return mode == PairMergeMode::reciprocal ? "reciprocal" : "union";
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t begin = 0;
    while (true) {
        const std::size_t tab = line.find('\t', begin);
        if (tab == std::string::npos) {
            fields.emplace_back(line.substr(begin));
            break;
        }
        fields.emplace_back(line.substr(begin, tab - begin));
        begin = tab + 1;
    }
    return fields;
}

std::uint64_t parse_u64(const std::string& text, const std::string& field) {
    try {
        std::size_t parsed = 0;
        const auto value = std::stoull(text, &parsed);
        if (text.empty() || parsed != text.size()) {
            throw std::out_of_range("bad integer");
        }
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid integer in " + field + ": " + text);
    }
}

std::unordered_map<std::string, std::size_t>
header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, std::size_t> index;
    for (std::size_t i = 0; i < header.size(); ++i) {
        index.emplace(header[i], i);
    }
    return index;
}

std::size_t require_column(
    const std::unordered_map<std::string, std::size_t>& columns,
    const std::string& name) {
    const auto found = columns.find(name);
    if (found == columns.end()) {
        throw std::runtime_error("missing column: " + name);
    }
    return found->second;
}

}  // namespace vallescope2
