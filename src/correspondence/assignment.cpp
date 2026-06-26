#include "vallescope2/correspondence/assignment.hpp"

#include "assignment_internal.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace vallescope2 {
namespace {

using AnchorTable = std::unordered_map<std::string, std::vector<const Anchor*>>;
using TmusTable = std::unordered_map<std::uint32_t, std::vector<const Anchor*>>;
using TmusIndex = std::unordered_map<std::string, TmusTable>;

void add_primary(std::vector<DirectedPrimary>& primaries,
                 const std::uint64_t pair_index,
                 const Pair& pair,
                 const Anchor& query,
                 const Candidate& candidate,
                 const std::string& method) {
    primaries.push_back(
        {pair_index, pair.target_sample, pair.query_sample, query.anchor_id,
         candidate.target->anchor_id, method, candidate.strand,
         candidate.edit_distance, candidate.score});
}

bool try_exact_context_match(std::ofstream& output,
                             const AnchorTable& by_context,
                             const std::uint64_t pair_index,
                             const Pair& pair,
                             const std::string& target_sample_name,
                             const std::string& query_sample_name,
                             const Anchor& query,
                             std::vector<DirectedPrimary>& primaries,
                             AssignmentResult& result) {
    const auto exact_found = by_context.find(query.context_key);
    if (exact_found == by_context.end() || exact_found->second.size() != 1) {
        return false;
    }

    const auto* target = exact_found->second.front();
    if (target->anchor_id == query.anchor_id) {
        write_unmatched(output, pair_index, target_sample_name,
                        query_sample_name, query);
        ++result.unmatched_count;
        return true;
    }

    Candidate candidate;
    candidate.target = target;
    candidate.edit_distance = 0;
    candidate.strand = '.';
    candidate.q_alpha_prime = query.alpha_prime[pair.target_sample];
    candidate.r_alpha_prime = target->alpha_prime[pair.query_sample];
    candidate.score = static_cast<std::int64_t>(candidate.q_alpha_prime) +
                      static_cast<std::int64_t>(candidate.r_alpha_prime);
    candidate.shared_tmus_count =
        count_shared_tmus(query.tmus_keys, target->tmus_keys);

    if (candidate.q_alpha_prime == 0 && candidate.r_alpha_prime == 0) {
        return false;
    }

    write_assignment(output, pair_index, target_sample_name, query_sample_name,
                     query, candidate, "primary", "exact_context", 0, 1);
    add_primary(primaries, pair_index, pair, query, candidate, "exact_context");
    ++result.primary_count;
    ++result.exact_context_primary_count;
    return true;
}

bool try_shared_tmus_match(std::ofstream& output,
                           const std::vector<const Anchor*>& targets,
                           const TmusIndex& by_group_tmus,
                           const std::uint64_t pair_index,
                           const Pair& pair,
                           const std::string& target_sample_name,
                           const std::string& query_sample_name,
                           const Anchor& query,
                           const std::uint64_t q_alpha_prime,
                           const AssignmentParameters& parameters,
                           std::vector<DirectedPrimary>& primaries,
                           AssignmentResult& result) {
    if (parameters.min_shared_tmus == 0) {
        // Preserve the old all-target scan behavior for this unusual setting.
        Candidate best_shared;
        bool has_shared = false;
        bool shared_tie = false;
        for (const auto* target : targets) {
            if (target->anchor_id == query.anchor_id) continue;
            const auto r_alpha_prime = target->alpha_prime[pair.query_sample];
            if (r_alpha_prime == 0) continue;
            const auto shared =
                count_shared_tmus(query.tmus_keys, target->tmus_keys);
            if (!has_shared || shared > best_shared.shared_tmus_count) {
                has_shared = true;
                shared_tie = false;
                best_shared = Candidate{};
                best_shared.target = target;
                best_shared.edit_distance = 0;
                best_shared.strand = '.';
                best_shared.q_alpha_prime = q_alpha_prime;
                best_shared.r_alpha_prime = r_alpha_prime;
                best_shared.shared_tmus_count = shared;
                best_shared.score = static_cast<std::int64_t>(q_alpha_prime) +
                                    static_cast<std::int64_t>(r_alpha_prime);
            } else if (shared == best_shared.shared_tmus_count) {
                shared_tie = true;
            }
        }
        if (has_shared && !shared_tie) {
            write_assignment(output, pair_index, target_sample_name,
                             query_sample_name, query, best_shared, "primary",
                             "shared_tmus", 0, 1);
            add_primary(primaries, pair_index, pair, query, best_shared,
                        "shared_tmus");
            ++result.primary_count;
            ++result.shared_tmus_primary_count;
            return true;
        }
        if (has_shared && shared_tie) ++result.shared_tmus_tie_fallback_count;
        return false;
    }

    const auto group_found = by_group_tmus.find(query.anchor_group_id);
    if (group_found == by_group_tmus.end()) return false;

    std::unordered_map<const Anchor*, std::uint32_t> shared_counts;
    for (const auto& key : query.tmus_keys) {
        const auto key_found = group_found->second.find(key);
        if (key_found == group_found->second.end()) continue;
        for (const auto* target : key_found->second) ++shared_counts[target];
    }

    Candidate best_shared;
    bool has_shared = false;
    bool shared_tie = false;

    for (const auto& item : shared_counts) {
        const auto* target = item.first;
        if (target->anchor_id == query.anchor_id) continue;
        const auto r_alpha_prime = target->alpha_prime[pair.query_sample];
        if (r_alpha_prime == 0) continue;
        const auto shared = item.second;
        if (shared < parameters.min_shared_tmus) continue;

        if (!has_shared || shared > best_shared.shared_tmus_count) {
            has_shared = true;
            shared_tie = false;
            best_shared = Candidate{};
            best_shared.target = target;
            best_shared.edit_distance = 0;
            best_shared.strand = '.';
            best_shared.q_alpha_prime = q_alpha_prime;
            best_shared.r_alpha_prime = r_alpha_prime;
            best_shared.shared_tmus_count = shared;
            best_shared.score = static_cast<std::int64_t>(q_alpha_prime) +
                                static_cast<std::int64_t>(r_alpha_prime);
        } else if (shared == best_shared.shared_tmus_count) {
            shared_tie = true;
        }
    }

    if (has_shared && !shared_tie) {
        write_assignment(output, pair_index, target_sample_name,
                         query_sample_name, query, best_shared, "primary",
                         "shared_tmus", 0, 1);
        add_primary(primaries, pair_index, pair, query, best_shared,
                    "shared_tmus");
        ++result.primary_count;
        ++result.shared_tmus_primary_count;
        return true;
    }
    if (has_shared && shared_tie) ++result.shared_tmus_tie_fallback_count;
    return false;
}

std::vector<Candidate> collect_edit_distance_candidates(
    const std::vector<const Anchor*>& targets,
    const Pair& pair,
    const Anchor& query,
    const std::uint64_t q_alpha_prime,
    const AssignmentParameters& parameters) {
    std::vector<Candidate> candidates;
    const auto reverse_query = reversed_tokens(query.forward_tokens);
    for (const auto* target : targets) {
        if (target->anchor_id == query.anchor_id) continue;
        const auto r_alpha_prime = target->alpha_prime[pair.query_sample];
        if (r_alpha_prime == 0) continue;

        const auto maximum_passing_edit_distance =
            static_cast<std::int64_t>(q_alpha_prime) +
            static_cast<std::int64_t>(r_alpha_prime) -
            parameters.min_candidate_score - 1;
        if (maximum_passing_edit_distance < 0) continue;

        const auto maximum_possible_edit_distance =
            std::max(query.forward_tokens.size(), target->forward_tokens.size());
        const auto edit_distance_limit =
            static_cast<std::uint32_t>(std::min<std::uint64_t>(
                static_cast<std::uint64_t>(maximum_passing_edit_distance),
                maximum_possible_edit_distance));
        auto candidate = make_candidate(query, reverse_query, *target,
                                        q_alpha_prime, r_alpha_prime,
                                        edit_distance_limit);
        candidate.shared_tmus_count =
            count_shared_tmus(query.tmus_keys, target->tmus_keys);
        if (candidate.score > parameters.min_candidate_score) {
            candidates.push_back(candidate);
        }
    }
    return candidates;
}

void write_edit_distance_result(std::ofstream& output,
                                std::vector<Candidate>& candidates,
                                const std::uint64_t pair_index,
                                const Pair& pair,
                                const std::string& target_sample_name,
                                const std::string& query_sample_name,
                                const Anchor& query,
                                const AssignmentParameters& parameters,
                                std::vector<DirectedPrimary>& primaries,
                                AssignmentResult& result) {
    std::sort(candidates.begin(), candidates.end(), candidate_less);
    if (candidates.size() == 1) {
        write_assignment(output, pair_index, target_sample_name,
                         query_sample_name, query, candidates.front(),
                         "primary", "edit_distance", 0, 1);
        add_primary(primaries, pair_index, pair, query, candidates.front(),
                    "edit_distance");
        ++result.primary_count;
        ++result.edit_distance_primary_count;
        return;
    }

    const auto margin = candidates[0].score - candidates[1].score;
    if (margin >= parameters.primary_margin) {
        write_assignment(output, pair_index, target_sample_name,
                         query_sample_name, query, candidates.front(),
                         "primary", "edit_distance", margin,
                         candidates.size());
        add_primary(primaries, pair_index, pair, query, candidates.front(),
                    "edit_distance");
        ++result.primary_count;
        ++result.edit_distance_primary_count;
        return;
    }

    ++result.ambiguous_query_count;
    for (const auto& candidate : candidates) {
        if (candidates.front().score - candidate.score >
            parameters.primary_margin) {
            break;
        }
        write_assignment(output, pair_index, target_sample_name,
                         query_sample_name, query, candidate, "ambiguous",
                         "edit_distance", margin, candidates.size());
        ++result.ambiguous_row_count;
    }
}

}  // namespace

AssignmentResult assign_context_correspondences(
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& anchor_contexts,
    const std::filesystem::path& sequence_table,
    const std::filesystem::path& assignment_output,
    const std::filesystem::path& correspondence_output,
    const std::filesystem::path& correspondence_metadata_output,
    const std::filesystem::path& metadata_output,
    const AssignmentParameters& parameters) {
    const auto data =
        load_assignment_data(grouped_anchors, anchor_contexts, sequence_table);
    const auto pairs = make_ordered_pairs(data);

    std::ofstream output(assignment_output);
    if (!output) throw std::runtime_error("cannot create assignment output");
    write_assignment_header(output);

    AssignmentResult result;
    result.ordered_pair_count = pairs.size();
    std::vector<DirectedPrimary> primaries;

    for (std::uint64_t pair_index = 0; pair_index < pairs.size(); ++pair_index) {
        const auto pair = pairs[pair_index];
        const auto& target_sample_name = data.samples[pair.target_sample];
        const auto& query_sample_name = data.samples[pair.query_sample];

        AnchorTable by_context;
        AnchorTable by_group;
        TmusIndex by_group_tmus;
        for (const auto* target : data.anchors_by_sample[pair.target_sample]) {
            by_context[target->context_key].push_back(target);
            by_group[target->anchor_group_id].push_back(target);
            for (const auto& key : target->tmus_keys) {
                by_group_tmus[target->anchor_group_id][key].push_back(target);
            }
        }

        for (const auto* query_ptr : data.anchors_by_sample[pair.query_sample]) {
            const auto& query = *query_ptr;
            ++result.query_anchor_count;

            if (try_exact_context_match(output, by_context, pair_index, pair,
                                        target_sample_name, query_sample_name,
                                        query, primaries, result)) {
                continue;
            }

            const auto q_alpha_prime = query.alpha_prime[pair.target_sample];
            const auto group_found = by_group.find(query.anchor_group_id);
            if (q_alpha_prime > 0 && group_found != by_group.end()) {
                const auto& targets = group_found->second;
                if (try_shared_tmus_match(output, targets, by_group_tmus,
                                          pair_index, pair, target_sample_name,
                                          query_sample_name, query,
                                          q_alpha_prime, parameters, primaries,
                                          result)) {
                    continue;
                }

                auto candidates = collect_edit_distance_candidates(
                    targets, pair, query, q_alpha_prime, parameters);
                if (!candidates.empty()) {
                    write_edit_distance_result(
                        output, candidates, pair_index, pair, target_sample_name,
                        query_sample_name, query, parameters, primaries, result);
                    continue;
                }
            }

            write_unmatched(output, pair_index, target_sample_name,
                            query_sample_name, query);
            ++result.unmatched_count;
        }
    }

    write_correspondences(correspondence_output, correspondence_metadata_output,
                          data, primaries, parameters, result);
    write_metadata(metadata_output, parameters, result);
    return result;
}

}  // namespace vallescope2
