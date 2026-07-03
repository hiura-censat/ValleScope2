#pragma once

#include "vallescope2/bundle/chaining.hpp"

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace vallescope2 {
namespace chaining_detail {

struct AnchorPosition {
    std::string sequence_id;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::uint64_t center = 0;
};

struct Candidate {
    std::string sample_a;
    std::string sample_b;
    std::string anchor_a;
    std::string anchor_b;
    std::string support_direction;
    std::string sequence_a;
    std::string sequence_b;
    std::uint64_t ref_center = 0;
    std::uint64_t query_center = 0;
    double score = 0.0;
    char strand = '+';
};

struct ChainState {
    double dp = 0.0;
    int predecessor = -1;
};

struct EmittedChain {
    std::uint64_t id = 0;
    std::string sample_a;
    std::string sample_b;
    std::string sequence_a;
    std::string sequence_b;
    char strand = '+';
    std::uint64_t ref_start = 0;
    std::uint64_t ref_end = 0;
    std::uint64_t query_start = 0;
    std::uint64_t query_end = 0;
    double score = 0.0;
    std::uint64_t both_anchor_count = 0;
    std::vector<Candidate> anchors;
};

std::unordered_map<std::string, AnchorPosition>
load_anchor_positions(const std::filesystem::path& path);

std::vector<Candidate> load_candidates(
    const std::filesystem::path& path,
    const std::unordered_map<std::string, AnchorPosition>& positions);

std::vector<Candidate> load_assignment_candidates(
    const std::filesystem::path& path,
    const std::unordered_map<std::string, AnchorPosition>& positions);

std::string group_key(const Candidate& candidate);
std::string candidate_key(const Candidate& candidate);

std::vector<std::size_t> chain_group(std::vector<Candidate>& candidates,
                                     const ChainingParameters& parameters,
                                     double& best_score);

EmittedChain make_emitted_chain(std::uint64_t chain_id,
                                const std::vector<Candidate>& group,
                                const std::vector<std::size_t>& chain,
                                double chain_score);

std::vector<EmittedChain> trim_overlapping_chains(
    std::vector<EmittedChain> chains,
    const ChainingParameters& parameters);

void write_chain(std::ostream& chains,
                 std::ostream& chain_anchors,
                 const EmittedChain& chain);

void write_refined_chain(std::ostream& chains,
                         std::ostream& chain_anchors,
                         const EmittedChain& chain,
                         const std::string& refinement_type,
                         const std::string& parent_chain_id,
                         const std::string& left_chain_id,
                         const std::string& right_chain_id);

void write_metadata(const std::filesystem::path& path,
                    const ChainingParameters& parameters,
                    const ChainingResult& result);

std::vector<EmittedChain> refine_chains(
    const std::vector<std::pair<std::string, std::vector<Candidate>>>& leftovers,
    std::vector<EmittedChain> primary_chains,
    const std::filesystem::path& refined_chain_output,
    const std::filesystem::path& refined_chain_anchor_output,
    const ChainingParameters& parameters,
    ChainingResult& result);

}  // namespace chaining_detail
}  // namespace vallescope2
