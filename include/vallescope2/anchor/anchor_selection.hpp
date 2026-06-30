#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <vector>

namespace vallescope2 {

struct AnchorSelectionParameters {
    std::uint32_t min_center_distance = 50;
    std::uint32_t target_density = 15000;
};

struct AnchorCandidate {
    std::uint64_t start = 0;
    double score = 0.0;
    std::uint64_t distance_from_group_center = 0;
    std::uint32_t sequence_id = 0;
    std::uint32_t length = 0;
};

class CandidateSpool {
public:
    CandidateSpool(const std::filesystem::path& directory,
                   std::uint32_t min_center_distance);
    ~CandidateSpool();
    CandidateSpool(CandidateSpool&&) noexcept;
    CandidateSpool& operator=(CandidateSpool&&) noexcept;
    CandidateSpool(const CandidateSpool&) = delete;
    CandidateSpool& operator=(const CandidateSpool&) = delete;

    void add(const AnchorCandidate& candidate);
    void finish();
    const std::filesystem::path& candidate_path() const;
    const std::filesystem::path& group_path() const;
    std::uint64_t size() const;

private:
    struct Implementation;
    std::unique_ptr<Implementation> implementation_;
};

struct AnchorSelectionResult {
    std::vector<AnchorCandidate> anchors;
    std::uint64_t target_count = 0;
    double actual_density = 0.0;
    double score_cutoff = 0.0;
    bool has_score_cutoff = false;
};

AnchorSelectionResult select_anchors_external(
    CandidateSpool& spool,
    std::size_t sequence_count,
    std::uint64_t sequence_length,
    const AnchorSelectionParameters& parameters);

}  // namespace vallescope2
