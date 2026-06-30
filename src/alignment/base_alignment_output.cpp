#include "base_alignment_internal.hpp"

#include <fstream>
#include <stdexcept>

namespace vallescope2 {
namespace alignment_detail {

void write_metadata(const std::filesystem::path& path,
                    const BaseAlignmentParameters& parameters,
                    const BaseAlignmentResult& result,
                    const std::vector<std::filesystem::path>& chain_files,
                    const std::filesystem::path& fasta) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create base alignment metadata");
    output << "{\n"
           << "  \"input_chains\": [";
    for (std::size_t i = 0; i < chain_files.size(); ++i) {
        if (i) output << ", ";
        output << "\"" << chain_files[i].string() << "\"";
    }
    output << "],\n"
           << "  \"fasta\": \"" << fasta.string() << "\",\n"
           << "  \"alignment_mode\": \"anchor_guided\",\n"
           << "  \"aligner\": \"WFA2-lib\",\n"
           << "  \"anchor_length\": " << parameters.anchor_length << ",\n"
           << "  \"max_bundle_align_bp\": " << parameters.max_bundle_align_bp << ",\n"
           << "  \"max_fallback_cells\": " << parameters.max_fallback_cells << ",\n"
           << "  \"max_patch_gap_bp\": " << parameters.max_patch_gap_bp << ",\n"
           << "  \"patch_window_bp\": " << parameters.patch_window_bp << ",\n"
           << "  \"min_patch_identity\": " << parameters.min_patch_identity << ",\n"
           << "  \"max_wfa_memory_gb\": " << parameters.max_wfa_memory_gb << ",\n"
           << "  \"bundle_trim_overlap\": " << parameters.bundle_trim_overlap << ",\n"
           << "  \"loaded_bundle_count\": " << result.loaded_bundle_count << ",\n"
           << "  \"post_patch_bundle_count\": " << result.post_patch_bundle_count << ",\n"
           << "  \"bundle_trim_count\": " << result.bundle_trim_count << ",\n"
           << "  \"bundle_count\": " << result.bundle_count << ",\n"
           << "  \"patched_bundle_count\": " << result.patched_bundle_count << ",\n"
           << "  \"patch_count\": " << result.patch_count << ",\n"
           << "  \"aligned_bundle_count\": " << result.aligned_bundle_count << ",\n"
           << "  \"anchor_guided_alignment_count\": "
           << result.anchor_guided_alignment_count << ",\n"
           << "  \"strand_flipped_bundle_count\": "
           << result.strand_flipped_bundle_count << ",\n"
           << "  \"strand_conflict_bundle_count\": "
           << result.strand_conflict_bundle_count << ",\n"
           << "  \"forward_anchor_match_count\": "
           << result.forward_anchor_match_count << ",\n"
           << "  \"reverse_anchor_match_count\": "
           << result.reverse_anchor_match_count << ",\n"
           << "  \"incompatible_anchor_match_count\": "
           << result.incompatible_anchor_match_count << ",\n"
           << "  \"no_valid_anchor_bundle_count\": "
           << result.no_valid_anchor_bundle_count << ",\n"
           << "  \"segment_alignment_failed_count\": "
           << result.segment_alignment_failed_count << ",\n"
           << "  \"skipped_bundle_count\": " << result.skipped_bundle_count << "\n"
           << "}\n";
}

}  // namespace alignment_detail
}  // namespace vallescope2
