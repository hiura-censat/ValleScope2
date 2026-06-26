#pragma once

#include "vallescope2/correspondence/assignment.hpp"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace vallescope2 {

struct GroupedAnchor {
    std::string sequence_id;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::string anchor_group_id;
    bool groupable = false;
};

struct Anchor {
    std::string anchor_id;
    std::string sample;
    std::string sequence_id;
    std::uint64_t start = 0;
    std::uint64_t end = 0;
    std::string anchor_group_id;
    std::string context_key;
    std::vector<std::uint32_t> forward_tokens;
    std::vector<std::uint32_t> tmus_keys;
    std::vector<std::uint64_t> alpha_prime;
};

struct Pair {
    std::uint32_t target_sample = 0;
    std::uint32_t query_sample = 0;
};

struct Candidate {
    const Anchor* target = nullptr;
    std::uint32_t edit_distance = 0;
    char strand = '.';
    std::int64_t score = 0;
    std::uint64_t q_alpha_prime = 0;
    std::uint64_t r_alpha_prime = 0;
    std::uint32_t shared_tmus_count = 0;
};

struct DirectedPrimary {
    std::uint64_t pair_id = 0;
    std::uint32_t target_sample = 0;
    std::uint32_t query_sample = 0;
    std::string query_anchor_id;
    std::string target_anchor_id;
    std::string method;
    char strand = '.';
    std::uint32_t edit_distance = 0;
    std::int64_t score = 0;
};

struct CanonicalPrimary {
    const DirectedPrimary* record = nullptr;
    std::string anchor_a;
    std::string anchor_b;
    bool forward = true;
};

struct DataSet {
    std::vector<std::string> samples;
    std::unordered_map<std::string, std::uint32_t> sample_ids;
    std::unordered_map<std::string, std::string> sample_roles;
    std::vector<Anchor> anchors;
    std::vector<std::vector<const Anchor*>> anchors_by_sample;
    std::unordered_map<std::string, std::uint32_t> token_ids;
    std::unordered_map<std::string, std::uint32_t> tmus_key_ids;
};

std::string pair_merge_mode_name(PairMergeMode mode);

std::vector<std::string> split_tab(const std::string& line);
std::uint64_t parse_u64(const std::string& text, const std::string& field);
std::unordered_map<std::string, std::size_t>
header_index(const std::vector<std::string>& header);
std::size_t require_column(
    const std::unordered_map<std::string, std::size_t>& columns,
    const std::string& name);

DataSet load_assignment_data(const std::filesystem::path& grouped_anchors,
                             const std::filesystem::path& anchor_contexts,
                             const std::filesystem::path& sequence_table);
std::vector<Pair> make_ordered_pairs(const DataSet& data);

std::uint32_t count_shared_tmus(const std::vector<std::uint32_t>& left,
                                const std::vector<std::uint32_t>& right);
std::vector<std::uint32_t> reversed_tokens(
    const std::vector<std::uint32_t>& tokens);
Candidate make_candidate(const Anchor& query,
                         const std::vector<std::uint32_t>& reverse_query,
                         const Anchor& target,
                         std::uint64_t q_alpha_prime,
                         std::uint64_t r_alpha_prime,
                         std::uint32_t edit_distance_limit);
bool candidate_less(const Candidate& left, const Candidate& right);

void write_assignment_header(std::ofstream& output);
void write_unmatched(std::ofstream& output,
                     std::uint64_t pair_id,
                     const std::string& target_sample,
                     const std::string& query_sample,
                     const Anchor& query);
void write_assignment(std::ofstream& output,
                      std::uint64_t pair_id,
                      const std::string& target_sample,
                      const std::string& query_sample,
                      const Anchor& query,
                      const Candidate& candidate,
                      const std::string& status,
                      const std::string& method,
                      std::int64_t score_margin,
                      std::uint64_t candidate_count);
void write_metadata(const std::filesystem::path& path,
                    const AssignmentParameters& parameters,
                    const AssignmentResult& result);

void write_correspondences(const std::filesystem::path& output_path,
                           const std::filesystem::path& metadata_path,
                           const DataSet& data,
                           const std::vector<DirectedPrimary>& primaries,
                           const AssignmentParameters& parameters,
                           AssignmentResult& result);

}  // namespace vallescope2
