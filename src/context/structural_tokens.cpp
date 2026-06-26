#include "vallescope2/context/structural_tokens.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

namespace vallescope2 {
namespace {

struct GroupedAnchorRow {
    std::string anchor_id;
    std::string sequence_id;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::string group_id;
    bool groupable = false;
};

std::uint64_t parse_coordinate(const std::string& value,
                               const std::uint64_t line_number) {
    try {
        std::size_t parsed = 0;
        const auto coordinate = std::stoull(value, &parsed);
        if (parsed != value.size()) throw std::invalid_argument("trailing text");
        return coordinate;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid grouped anchor coordinate at line " +
                                 std::to_string(line_number));
    }
}

GroupedAnchorRow parse_row(const std::string& line,
                           const std::uint64_t line_number) {
    std::istringstream input(line);
    GroupedAnchorRow row;
    std::string start;
    std::string end;
    std::string anchor_sequence;
    std::string canonical_sequence;
    std::string orientation;
    std::string groupable;
    if (!(input >> row.anchor_id >> row.sequence_id >> start >> end >>
          anchor_sequence >> canonical_sequence >> row.group_id >>
          orientation >> groupable)) {
        throw std::runtime_error("invalid grouped anchor row at line " +
                                 std::to_string(line_number));
    }
    row.start = parse_coordinate(start, line_number);
    row.end = parse_coordinate(end, line_number);
    if (row.start >= row.end) {
        throw std::runtime_error("invalid grouped anchor interval at line " +
                                 std::to_string(line_number));
    }
    if (groupable == "true") {
        row.groupable = true;
        if (row.group_id == ".") {
            throw std::runtime_error("groupable anchor has no group at line " +
                                     std::to_string(line_number));
        }
    } else if (groupable == "false") {
        if (row.group_id != ".") {
            throw std::runtime_error("ungroupable anchor has a group at line " +
                                     std::to_string(line_number));
        }
    } else {
        throw std::runtime_error("invalid groupable value at line " +
                                 std::to_string(line_number));
    }
    return row;
}

void write_metadata(const std::filesystem::path& path,
                    const std::uint32_t bin_size,
                    const StructuralTokenResult& result) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot create structural token metadata: " +
                                 path.string());
    }
    output << "{\n"
           << "  \"distance_definition\": \"center_to_center\",\n"
           << "  \"distance_bin_size\": " << bin_size << ",\n"
           << "  \"distance_bin_method\": \"floor\",\n"
           << "  \"n_anchor_tokens\": " << result.anchor_token_count << ",\n"
           << "  \"n_distance_tokens\": " << result.distance_token_count << ",\n"
           << "  \"n_tokens\": " << result.token_count() << ",\n"
           << "  \"n_sequences\": " << result.sequence_count << ",\n"
           << "  \"n_ungrouped_anchors_skipped\": "
           << result.skipped_anchor_count << "\n"
           << "}\n";
}

}  // namespace

StructuralTokenResult build_structural_tokens(
    const std::filesystem::path& grouped_anchors,
    const std::uint32_t distance_bin_size,
    const std::filesystem::path& token_output,
    const std::filesystem::path& metadata_output) {
    if (distance_bin_size == 0) {
        throw std::runtime_error("distance bin size must be greater than zero");
    }
    std::ifstream input(grouped_anchors);
    if (!input) {
        throw std::runtime_error("cannot open grouped anchors: " +
                                 grouped_anchors.string());
    }
    std::ofstream output(token_output);
    if (!output) {
        throw std::runtime_error("cannot create structural token output: " +
                                 token_output.string());
    }
    output << "sequence_id\ttoken_index\ttoken_type\ttoken\tanchor_id"
              "\traw_distance\tdistance_bin\n";

    StructuralTokenResult result;
    std::unordered_set<std::string> completed_sequences;
    std::string current_sequence;
    std::uint64_t previous_start = 0;
    std::uint64_t previous_center = 0;
    bool has_previous_row = false;
    bool has_previous_grouped_anchor = false;
    std::uint64_t token_index = 0;
    std::string line;
    std::uint64_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line_number == 1 && line.rfind("anchor_id\t", 0) == 0) continue;
        if (line.empty() || line.front() == '#') continue;
        const auto row = parse_row(line, line_number);

        if (row.sequence_id != current_sequence) {
            if (!current_sequence.empty()) completed_sequences.insert(current_sequence);
            if (completed_sequences.count(row.sequence_id) != 0) {
                throw std::runtime_error(
                    "grouped anchors are not contiguous by sequence at line " +
                    std::to_string(line_number));
            }
            current_sequence = row.sequence_id;
            has_previous_row = false;
            has_previous_grouped_anchor = false;
            token_index = 0;
        }
        if (has_previous_row && row.start < previous_start) {
            throw std::runtime_error(
                "grouped anchors are not coordinate-sorted at line " +
                std::to_string(line_number));
        }
        previous_start = row.start;
        has_previous_row = true;

        if (!row.groupable) {
            ++result.skipped_anchor_count;
            continue;
        }
        const std::uint64_t center = row.start + (row.end - row.start) / 2;
        if (has_previous_grouped_anchor) {
            if (center <= previous_center) {
                throw std::runtime_error(
                    "groupable anchor centers are not increasing at line " +
                    std::to_string(line_number));
            }
            const std::uint64_t distance = center - previous_center;
            const std::uint64_t distance_bin = distance / distance_bin_size;
            output << row.sequence_id << '\t' << token_index++
                   << "\tD\tD:" << distance_bin << "\t.\t" << distance
                   << '\t' << distance_bin << '\n';
            ++result.distance_token_count;
        } else {
            ++result.sequence_count;
        }
        output << row.sequence_id << '\t' << token_index++ << "\tA\tA:"
               << row.group_id << '\t' << row.anchor_id << "\t.\t.\n";
        ++result.anchor_token_count;
        previous_center = center;
        has_previous_grouped_anchor = true;
    }
    if (!output) throw std::runtime_error("failed to write structural tokens");
    write_metadata(metadata_output, distance_bin_size, result);
    return result;
}

}  // namespace vallescope2
