#include "vallescope2/context/anchor_grouping.hpp"

#include "vallescope2/common/sha256.hpp"

#include <htslib/faidx.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace vallescope2 {
namespace {

struct FaiDeleter {
    void operator()(faidx_t* index) const { fai_destroy(index); }
};

struct SequenceDeleter {
    void operator()(char* sequence) const { std::free(sequence); }
};

struct GroupSummary {
    std::string hash;
    std::string canonical_sequence;
    std::uint64_t count = 0;
};

char complement(const char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'R': return 'Y';
        case 'Y': return 'R';
        case 'S': return 'S';
        case 'W': return 'W';
        case 'K': return 'M';
        case 'M': return 'K';
        case 'B': return 'V';
        case 'V': return 'B';
        case 'D': return 'H';
        case 'H': return 'D';
        case 'N': return 'N';
        default: throw std::runtime_error("invalid nucleotide in anchor sequence");
    }
}

std::uint64_t parse_coordinate(const std::string& value,
                               const std::string& field,
                               const std::uint64_t line_number) {
    try {
        std::size_t parsed = 0;
        const auto coordinate = std::stoull(value, &parsed);
        if (parsed != value.size()) throw std::invalid_argument("trailing text");
        return coordinate;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid BED " + field + " at line " +
                                 std::to_string(line_number));
    }
}

std::string fetch_anchor(faidx_t* index,
                         const std::string& sequence_id,
                         const std::uint64_t start,
                         const std::uint64_t end,
                         const std::uint64_t line_number) {
    const hts_pos_t sequence_length = faidx_seq_len64(index, sequence_id.c_str());
    if (sequence_length < 0) {
        throw std::runtime_error("BED sequence is absent from FASTA at line " +
                                 std::to_string(line_number) + ": " + sequence_id);
    }
    if (start >= end || end > static_cast<std::uint64_t>(sequence_length) ||
        end - 1 > static_cast<std::uint64_t>(std::numeric_limits<hts_pos_t>::max())) {
        throw std::runtime_error("BED interval is outside FASTA at line " +
                                 std::to_string(line_number));
    }

    hts_pos_t fetched_length = 0;
    std::unique_ptr<char, SequenceDeleter> fetched(
        faidx_fetch_seq64(index, sequence_id.c_str(),
                          static_cast<hts_pos_t>(start),
                          static_cast<hts_pos_t>(end - 1), &fetched_length));
    if (!fetched || fetched_length != static_cast<hts_pos_t>(end - start)) {
        throw std::runtime_error("failed to extract anchor sequence at BED line " +
                                 std::to_string(line_number));
    }
    std::string sequence(fetched.get(), static_cast<std::size_t>(fetched_length));
    std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                   [](const unsigned char value) {
                       return static_cast<char>(std::toupper(value));
                   });
    return sequence;
}

}  // namespace

std::string reverse_complement(const std::string& sequence) {
    std::string result(sequence.size(), 'N');
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        result[sequence.size() - i - 1] =
            complement(static_cast<char>(std::toupper(
                static_cast<unsigned char>(sequence[i]))));
    }
    return result;
}

AnchorGroupingResult group_anchors(
    const std::filesystem::path& anchor_bed,
    const std::filesystem::path& fasta,
    const std::filesystem::path& fasta_index,
    const std::uint32_t anchor_length,
    const std::filesystem::path& grouped_anchors_output,
    const std::filesystem::path& anchor_groups_output) {
    std::ifstream bed(anchor_bed);
    if (!bed) throw std::runtime_error("cannot open anchor BED: " + anchor_bed.string());
    std::unique_ptr<faidx_t, FaiDeleter> index(
        fai_load3(fasta.c_str(), fasta_index.c_str(), nullptr, 0));
    if (!index) throw std::runtime_error("cannot load context FASTA index");
    std::ofstream anchors(grouped_anchors_output);
    if (!anchors) throw std::runtime_error("cannot create grouped anchor output");
    anchors << "anchor_id\tsequence_id\tstart\tend\tanchor_seq\tcanonical_seq"
               "\tanchor_group_id\torientation\tgroupable\n";

    AnchorGroupingResult result;
    std::unordered_set<std::string> anchor_ids;
    std::unordered_map<std::string, std::uint32_t> group_by_hash;
    std::vector<GroupSummary> groups;
    std::string line;
    std::uint64_t line_number = 0;
    while (std::getline(bed, line)) {
        ++line_number;
        if (line.empty() || line.front() == '#') continue;
        std::istringstream fields(line);
        std::string sequence_id;
        std::string start_text;
        std::string end_text;
        std::string anchor_id;
        if (!(fields >> sequence_id >> start_text >> end_text >> anchor_id)) {
            throw std::runtime_error("anchor BED needs at least four columns at line " +
                                     std::to_string(line_number));
        }
        if (!anchor_ids.insert(anchor_id).second) {
            throw std::runtime_error("duplicate anchor ID: " + anchor_id);
        }
        const std::uint64_t start = parse_coordinate(start_text, "start", line_number);
        const std::uint64_t end = parse_coordinate(end_text, "end", line_number);
        if (end < start || end - start != anchor_length) {
            throw std::runtime_error("anchor length mismatch at BED line " +
                                     std::to_string(line_number));
        }
        const std::string anchor_sequence =
            fetch_anchor(index.get(), sequence_id, start, end, line_number);
        const std::string reverse = reverse_complement(anchor_sequence);
        const std::string canonical = std::min(anchor_sequence, reverse);
        const char orientation = anchor_sequence == reverse
                                     ? '.'
                                     : (anchor_sequence == canonical ? '+' : '-');

        anchors << anchor_id << '\t' << sequence_id << '\t' << start << '\t'
                << end << '\t' << anchor_sequence << '\t' << canonical << '\t';
        if (anchor_sequence.find('N') != std::string::npos) {
            anchors << ".\t" << orientation << "\tfalse\n";
            ++result.ungrouped_anchor_count;
        } else {
            const std::string hash = sha256_text(canonical);
            auto found = group_by_hash.find(hash);
            std::uint32_t group_id = 0;
            if (found == group_by_hash.end()) {
                group_id = static_cast<std::uint32_t>(groups.size());
                group_by_hash.emplace(hash, group_id);
                groups.push_back({hash, canonical, 0});
            } else {
                group_id = found->second;
                if (groups[group_id].canonical_sequence != canonical) {
                    throw std::runtime_error("SHA-256 collision between anchor groups");
                }
            }
            ++groups[group_id].count;
            anchors << hash << '\t' << orientation << "\ttrue\n";
            ++result.grouped_anchor_count;
        }
        ++result.anchor_count;
    }
    if (!anchors) throw std::runtime_error("failed to write grouped anchors");

    std::ofstream group_output(anchor_groups_output);
    if (!group_output) throw std::runtime_error("cannot create anchor group output");
    group_output << "anchor_group_id\tcanonical_seq\tanchor_count\n";
    for (std::size_t i = 0; i < groups.size(); ++i) {
        group_output << groups[i].hash << '\t' << groups[i].canonical_sequence << '\t'
                     << groups[i].count << '\n';
    }
    result.group_count = groups.size();
    return result;
}

}  // namespace vallescope2
