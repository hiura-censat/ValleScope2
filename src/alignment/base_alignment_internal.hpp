#pragma once

#include "vallescope2/alignment/base_alignment.hpp"

#include <htslib/faidx.h>

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace vallescope2 {
namespace alignment_detail {

struct FaiDeleter {
    void operator()(faidx_t* index) const;
};

struct ChainBundle {
    std::string source;
    std::uint64_t chain_id = 0;
    double score = 0.0;
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
    std::uint64_t trim_count = 0;
    std::vector<std::size_t> anchor_indices;
};

struct AnchorPair {
    std::string anchor_a;
    std::string anchor_b;
    std::uint64_t ref_start = 0;
    std::uint64_t ref_end = 0;
    std::uint64_t query_start = 0;
    std::uint64_t query_end = 0;
};

struct ExtensionCandidate {
    std::string sample_a;
    std::string sample_b;
    std::string sequence_a;
    std::string sequence_b;
    std::string anchor_a;
    std::string anchor_b;
    std::string support_direction;
    std::uint64_t ref_start = 0;
    std::uint64_t ref_end = 0;
    std::uint64_t ref_center = 0;
    std::uint64_t query_start = 0;
    std::uint64_t query_end = 0;
    std::uint64_t query_center = 0;
    double score = 0.0;
    char strand = '+';
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

struct AnchorStore {
    std::vector<AnchorPair> anchors;
};

struct StrandSummary {
    std::uint64_t forward_count = 0;
    std::uint64_t reverse_count = 0;
    std::uint64_t incompatible_count = 0;
    char effective_strand = '+';
};

std::vector<std::string> split_tab(const std::string& line);
std::unordered_map<std::string, std::size_t>
header_index(const std::vector<std::string>& header);
std::size_t require_column(
    const std::unordered_map<std::string, std::size_t>& columns,
    const std::string& name);
std::uint64_t parse_u64(const std::string& value, const std::string& column);
double parse_f64(const std::string& value, const std::string& column);

std::vector<ChainBundle> load_chain_files(
    const std::vector<std::filesystem::path>& paths,
    const std::vector<std::filesystem::path>& anchor_paths,
    const std::filesystem::path& grouped_anchors,
    AnchorStore& store);
std::vector<ExtensionCandidate> load_extension_candidates(
    const std::filesystem::path& correspondences,
    const std::filesystem::path& assignments,
    const std::filesystem::path& grouped_anchors);

std::uint64_t sequence_length(faidx_t* index, const std::string& sequence_id);
Interval expand_interval(std::uint64_t start,
                         std::uint64_t end,
                         std::uint64_t sequence_length,
                         std::uint32_t anchor_length);
std::string fetch_interval(faidx_t* index,
                           const std::string& sequence_id,
                           const Interval& interval);
Alignment align_segment(const std::string& ref,
                        const std::string& query,
                        const BaseAlignmentParameters& parameters);
void append_cigar_operation(std::string& cigar,
                            std::uint64_t length,
                            char op);
void append_cigar_string(std::string& destination,
                         const std::string& source);
Interval oriented_to_original_query_interval(const Interval& oriented,
                                             std::uint64_t query_length,
                                             char strand);
std::string fetch_oriented_query_interval(faidx_t* index,
                                          const std::string& sequence_id,
                                          const Interval& oriented,
                                          std::uint64_t query_length,
                                          char strand);

StrandSummary summarize_anchor_strand(
    faidx_t* index,
    const ChainBundle& bundle,
    const AnchorStore& store,
    const Interval& ref_interval,
    const Interval& query_interval);
std::optional<Alignment> anchor_guided_alignment(
    faidx_t* index,
    const ChainBundle& bundle,
    const AnchorStore& store,
    const Interval& ref_interval,
    const Interval& query_interval,
    std::uint64_t query_length,
    char effective_strand,
    const BaseAlignmentParameters& parameters);

bool same_bundle_track(const ChainBundle& left, const ChainBundle& right);
std::vector<ChainBundle> extend_adjacent_bundles(
    std::vector<ChainBundle> bundles,
    const std::vector<ExtensionCandidate>& candidates,
    AnchorStore& store,
    const BaseAlignmentParameters& parameters,
    const std::filesystem::path& extension_output,
    const std::filesystem::path& extension_anchor_output,
    BaseAlignmentResult& result);
std::vector<ChainBundle> final_trim_bundles(
    std::vector<ChainBundle> bundles,
    const BaseAlignmentParameters& parameters,
    std::uint64_t& trim_count);
std::vector<ChainBundle> patch_adjacent_bundles(
    std::vector<ChainBundle> bundles,
    faidx_t* index,
    const BaseAlignmentParameters& parameters,
    std::uint64_t& patch_count);

void write_metadata(const std::filesystem::path& path,
                    const BaseAlignmentParameters& parameters,
                    const BaseAlignmentResult& result,
                    const std::vector<std::filesystem::path>& chain_files,
                    const std::filesystem::path& fasta);

}  // namespace alignment_detail
}  // namespace vallescope2
