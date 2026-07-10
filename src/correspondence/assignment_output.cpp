#include "assignment_internal.hpp"

#include <fstream>
#include <stdexcept>

namespace vallescope2 {

void write_assignment_header(std::ofstream& output) {
    output << "pair_id\ttarget_sample\tquery_sample\tquery_anchor_id"
              "\ttarget_anchor_id\tstatus\tmethod\tassign_strand"
              "\tedit_distance\tcandidate_score\tscore_margin\tn_candidates"
              "\tq_alpha_prime\tr_alpha_prime\tshared_tmus_count"
              "\tquery_anchor_group_id\ttarget_anchor_group_id"
              "\tquery_context_key\ttarget_context_key\n";
}

void write_unmatched(std::ofstream& output,
                     const std::uint64_t pair_id,
                     const std::string& target_sample,
                     const std::string& query_sample,
                     const Anchor& query) {
    output << pair_id << '\t' << target_sample << '\t' << query_sample << '\t'
           << query.anchor_id
           << "\t.\tunmatched\tunmatched\t.\t.\t.\t.\t0\t.\t.\t.\t"
           << query.anchor_group_id << "\t.\t" << query.context_key << "\t.\n";
}

void write_assignment(std::ofstream& output,
                      const std::uint64_t pair_id,
                      const std::string& target_sample,
                      const std::string& query_sample,
                      const Anchor& query,
                      const Candidate& candidate,
                      const std::string& status,
                      const std::string& method,
                      const std::int64_t score_margin,
                      const std::uint64_t candidate_count) {
    output << pair_id << '\t' << target_sample << '\t' << query_sample << '\t'
           << query.anchor_id << '\t' << candidate.target->anchor_id << '\t'
           << status << '\t' << method << '\t' << candidate.strand << '\t'
           << candidate.edit_distance << '\t' << candidate.score << '\t'
           << score_margin << '\t' << candidate_count << '\t'
           << candidate.q_alpha_prime << '\t' << candidate.r_alpha_prime << '\t'
           << candidate.shared_tmus_count << '\t' << query.anchor_group_id
           << '\t' << candidate.target->anchor_group_id << '\t'
           << query.context_key << '\t' << candidate.target->context_key << '\n';
}

void write_metadata(const std::filesystem::path& path,
                    const AssignmentParameters& parameters,
                    const AssignmentResult& result) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create assignment metadata");
    output << "{\n"
           << "  \"candidate_filter\": "
           << "\"candidate_score > min_candidate_score\",\n"
           << "  \"min_candidate_score\": " << parameters.min_candidate_score
           << ",\n"
           << "  \"min_shared_tmus\": " << parameters.min_shared_tmus
           << ",\n"
           << "  \"pair_merge_mode\": \""
           << pair_merge_mode_name(parameters.pair_merge_mode) << "\",\n"
           << "  \"primary_margin\": " << parameters.primary_margin << ",\n"
           << "  \"ordered_pair_count\": " << result.ordered_pair_count
           << ",\n"
           << "  \"query_anchor_count\": " << result.query_anchor_count
           << ",\n"
           << "  \"primary_count\": " << result.primary_count << ",\n"
           << "  \"ambiguous_query_count\": "
           << result.ambiguous_query_count << ",\n"
           << "  \"ambiguous_row_count\": " << result.ambiguous_row_count
           << ",\n"
           << "  \"unmatched_count\": " << result.unmatched_count << ",\n"
           << "  \"exact_context_primary_count\": "
           << result.exact_context_primary_count << ",\n"
           << "  \"shared_tmus_primary_count\": "
           << result.shared_tmus_primary_count << ",\n"
           << "  \"shared_tmus_tie_fallback_count\": "
           << result.shared_tmus_tie_fallback_count << ",\n"
           << "  \"edit_distance_primary_count\": "
           << result.edit_distance_primary_count << "\n"
           << "}\n";
}

}  // namespace vallescope2
