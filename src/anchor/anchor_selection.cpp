#include "vallescope2/anchor/anchor_selection.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <queue>
#include <set>
#include <stdexcept>

namespace vallescope2 {
namespace {

constexpr std::size_t kCandidatesPerChunk = 1000000;

struct GroupRecord {
    std::uint64_t begin = 0;
    std::uint64_t end = 0;
    std::uint64_t center_sum = 0;
};

template<class Value>
void write_binary(std::ofstream& output, const Value& value) {
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
    if (!output) throw std::runtime_error("failed to write candidate spool");
}

template<class Value>
bool read_binary(std::ifstream& input, Value& value) {
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (input) return true;
    if (input.eof() && input.gcount() == 0) return false;
    throw std::runtime_error("truncated candidate spool");
}

bool rank_order(const AnchorCandidate& left, const AnchorCandidate& right) {
    if (left.score != right.score) return left.score < right.score;
    if (left.distance_from_group_center != right.distance_from_group_center) {
        return left.distance_from_group_center < right.distance_from_group_center;
    }
    if (left.sequence_id != right.sequence_id) {
        return left.sequence_id < right.sequence_id;
    }
    return left.start < right.start;
}

bool conflicts(const std::set<std::uint64_t>& selected,
               const std::uint64_t center,
               const std::uint64_t minimum_distance) {
    const auto right = selected.lower_bound(center);
    if (right != selected.end() && *right - center <= minimum_distance) return true;
    if (right != selected.begin()) {
        const auto left = std::prev(right);
        if (center - *left <= minimum_distance) return true;
    }
    return false;
}

}  // namespace

struct CandidateSpool::Implementation {
    Implementation(const std::filesystem::path& directory,
                   const std::uint32_t distance)
        : candidate_path(directory / "anchor_candidates.bin"),
          group_path(directory / "anchor_groups.bin"),
          minimum_distance(distance),
          candidates(candidate_path, std::ios::binary),
          groups(group_path, std::ios::binary) {
        if (!candidates || !groups) {
            throw std::runtime_error("cannot create candidate spool");
        }
    }

    void close_group() {
        if (!has_group) return;
        write_binary(groups, GroupRecord{group_begin, count,
                                         group_first_start + group_last_start});
        has_group = false;
    }

    std::filesystem::path candidate_path;
    std::filesystem::path group_path;
    std::uint32_t minimum_distance;
    std::ofstream candidates;
    std::ofstream groups;
    std::uint64_t count = 0;
    bool has_group = false;
    std::uint64_t group_begin = 0;
    std::uint64_t group_first_start = 0;
    std::uint64_t group_last_start = 0;
    AnchorCandidate previous;
    bool finished = false;
};

CandidateSpool::CandidateSpool(const std::filesystem::path& directory,
                               const std::uint32_t min_center_distance)
    : implementation_(std::make_unique<Implementation>(directory,
                                                       min_center_distance)) {}
CandidateSpool::~CandidateSpool() = default;
CandidateSpool::CandidateSpool(CandidateSpool&&) noexcept = default;
CandidateSpool& CandidateSpool::operator=(CandidateSpool&&) noexcept = default;

void CandidateSpool::add(const AnchorCandidate& candidate) {
    auto& state = *implementation_;
    if (state.finished) throw std::runtime_error("candidate spool is closed");
    const bool same_group = state.has_group &&
        candidate.sequence_id == state.previous.sequence_id &&
        candidate.score == state.previous.score &&
        candidate.start - state.previous.start <= state.minimum_distance;
    if (!same_group) {
        state.close_group();
        state.has_group = true;
        state.group_begin = state.count;
        state.group_first_start = candidate.start;
    }
    state.group_last_start = candidate.start;
    state.previous = candidate;
    write_binary(state.candidates, candidate);
    ++state.count;
}

void CandidateSpool::finish() {
    auto& state = *implementation_;
    if (state.finished) return;
    state.close_group();
    state.candidates.close();
    state.groups.close();
    state.finished = true;
}

const std::filesystem::path& CandidateSpool::candidate_path() const {
    return implementation_->candidate_path;
}
const std::filesystem::path& CandidateSpool::group_path() const {
    return implementation_->group_path;
}
std::uint64_t CandidateSpool::size() const { return implementation_->count; }

AnchorSelectionResult select_anchors_external(
    CandidateSpool& spool,
    const std::size_t sequence_count,
    const std::uint64_t sequence_length,
    const AnchorSelectionParameters& parameters) {
    if (sequence_length == 0) {
        throw std::runtime_error("sequence length must be greater than zero");
    }
    spool.finish();

    AnchorSelectionResult result;
    result.target_count = static_cast<std::uint64_t>(std::floor(
        static_cast<long double>(sequence_length) * parameters.target_density /
        1000000.0L));
    if (result.target_count == 0 || spool.size() == 0) return result;

    std::ifstream candidates(spool.candidate_path(), std::ios::binary);
    std::ifstream groups(spool.group_path(), std::ios::binary);
    if (!candidates || !groups) throw std::runtime_error("cannot read candidate spool");

    std::vector<std::filesystem::path> chunk_paths;
    std::vector<AnchorCandidate> chunk;
    chunk.reserve(kCandidatesPerChunk);
    GroupRecord group;
    bool has_group = read_binary(groups, group);
    std::uint64_t index = 0;
    AnchorCandidate candidate;
    while (read_binary(candidates, candidate)) {
        while (has_group && index >= group.end) has_group = read_binary(groups, group);
        if (!has_group || index < group.begin) {
            throw std::runtime_error("candidate group spool is inconsistent");
        }
        const std::uint64_t doubled_start = candidate.start * 2;
        candidate.distance_from_group_center =
            doubled_start > group.center_sum
                ? doubled_start - group.center_sum
                : group.center_sum - doubled_start;
        chunk.push_back(candidate);
        ++index;
        if (chunk.size() == kCandidatesPerChunk) {
            std::sort(chunk.begin(), chunk.end(), rank_order);
            const auto path = spool.candidate_path().parent_path() /
                ("anchor_chunk_" + std::to_string(chunk_paths.size()) + ".bin");
            std::ofstream output(path, std::ios::binary);
            for (const auto& value : chunk) write_binary(output, value);
            chunk_paths.push_back(path);
            chunk.clear();
        }
    }
    if (!chunk.empty()) {
        std::sort(chunk.begin(), chunk.end(), rank_order);
        const auto path = spool.candidate_path().parent_path() /
            ("anchor_chunk_" + std::to_string(chunk_paths.size()) + ".bin");
        std::ofstream output(path, std::ios::binary);
        for (const auto& value : chunk) write_binary(output, value);
        chunk_paths.push_back(path);
    }

    struct StreamValue { AnchorCandidate candidate; std::size_t stream = 0; };
    const auto later = [](const StreamValue& left, const StreamValue& right) {
        return rank_order(right.candidate, left.candidate);
    };
    std::priority_queue<StreamValue, std::vector<StreamValue>, decltype(later)>
        queue(later);
    std::vector<std::ifstream> streams;
    streams.reserve(chunk_paths.size());
    for (std::size_t i = 0; i < chunk_paths.size(); ++i) {
        streams.emplace_back(chunk_paths[i], std::ios::binary);
        AnchorCandidate first;
        if (read_binary(streams.back(), first)) queue.push({first, i});
    }

    std::vector<std::set<std::uint64_t>> selected_centers(sequence_count);
    while (!queue.empty() && result.anchors.size() < result.target_count) {
        const auto value = queue.top();
        queue.pop();
        const auto& current = value.candidate;
        if (current.sequence_id >= selected_centers.size()) {
            throw std::runtime_error("invalid anchor sequence ID");
        }
        const std::uint64_t center = current.start + current.length / 2;
        auto& centers = selected_centers[current.sequence_id];
        if (!conflicts(centers, center, parameters.min_center_distance)) {
            centers.insert(center);
            result.anchors.push_back(current);
            result.score_cutoff = current.score;
            result.has_score_cutoff = true;
        }
        AnchorCandidate next;
        if (read_binary(streams[value.stream], next)) queue.push({next, value.stream});
    }

    streams.clear();
    for (const auto& path : chunk_paths) {
        std::error_code error;
        std::filesystem::remove(path, error);
    }
    std::sort(result.anchors.begin(), result.anchors.end(),
              [](const AnchorCandidate& left, const AnchorCandidate& right) {
                  if (left.sequence_id != right.sequence_id) {
                      return left.sequence_id < right.sequence_id;
                  }
                  return left.start < right.start;
              });
    result.actual_density = static_cast<double>(result.anchors.size()) * 1000000.0 /
                            static_cast<double>(sequence_length);
    return result;
}

}  // namespace vallescope2
