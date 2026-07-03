#include "chaining_internal.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace vallescope2 {
namespace chaining_detail {
namespace {

bool same_chain_group(const EmittedChain& chain,
                      const std::string& key) {
    Candidate candidate;
    candidate.sample_a = chain.sample_a;
    candidate.sample_b = chain.sample_b;
    candidate.sequence_a = chain.sequence_a;
    candidate.sequence_b = chain.sequence_b;
    candidate.strand = chain.strand;
    return group_key(candidate) == key;
}

std::uint64_t subtract_window(const std::uint64_t value,
                              const std::uint32_t window) {
    return value > window ? value - window : 0;
}

bool candidate_in_box(const Candidate& candidate,
                      const std::uint64_t ref_start,
                      const std::uint64_t ref_end,
                      const std::uint64_t query_start,
                      const std::uint64_t query_end,
                      const std::unordered_set<std::string>& used) {
    if (used.find(candidate_key(candidate)) != used.end()) return false;
    return ref_start <= candidate.ref_center && candidate.ref_center <= ref_end &&
           query_start <= candidate.query_center && candidate.query_center <= query_end;
}

bool emit_refined_subset(
    std::vector<Candidate> subset,
    const ChainingParameters& parameters,
    std::uint64_t& next_chain_id,
    ChainingResult& result,
    std::unordered_set<std::string>& used,
    std::vector<EmittedChain>& refined_chains,
    std::ostream& chains,
    std::ostream& chain_anchors,
    const std::string& refinement_type,
    const std::string& parent_chain_id,
    const std::string& left_chain_id,
    const std::string& right_chain_id) {
    if (subset.size() < parameters.refinement_min_chain_anchors) return false;
    double chain_score = 0.0;
    const auto chain = chain_group(subset, parameters, chain_score);
    if (chain.size() < parameters.refinement_min_chain_anchors ||
        chain_score < parameters.min_chain_score) {
        return false;
    }
    auto emitted = make_emitted_chain(next_chain_id++, subset, chain, chain_score);
    if (emitted.both_anchor_count < parameters.min_chain_both_anchors) {
        return false;
    }
    write_refined_chain(chains, chain_anchors, emitted, refinement_type,
                        parent_chain_id, left_chain_id, right_chain_id);
    for (const auto& candidate : emitted.anchors) used.insert(candidate_key(candidate));
    ++result.refined_chain_count;
    result.refined_chain_anchor_count += emitted.anchors.size();
    refined_chains.push_back(std::move(emitted));
    return true;
}

}  // namespace

std::vector<EmittedChain> refine_chains(
    const std::vector<std::pair<std::string, std::vector<Candidate>>>& leftovers,
    std::vector<EmittedChain> primary_chains,
    const std::filesystem::path& refined_chain_output,
    const std::filesystem::path& refined_chain_anchor_output,
    const ChainingParameters& parameters,
    ChainingResult& result) {
    std::ofstream chains(refined_chain_output);
    if (!chains) throw std::runtime_error("cannot create refined chain output");
    std::ofstream chain_anchors(refined_chain_anchor_output);
    if (!chain_anchors) throw std::runtime_error("cannot create refined chain anchor output");
    chains << "chain_id\trefinement_type\tparent_chain_id\tleft_chain_id"
              "\tright_chain_id\tsample_a\tsample_b\tsequence_a\tsequence_b"
              "\tassign_strand\tn_candidates\tboth_anchor_count\tboth_rate"
              "\tchain_score"
              "\tref_start\tref_end\tquery_start\tquery_end\n";
    chain_anchors << "chain_id\trank\trefinement_type\tsample_a\tsample_b"
                     "\tanchor_a\tanchor_b\tref_center\tquery_center"
                     "\tassign_strand\tcandidate_score\tsupport_direction\n";

    std::uint64_t next_chain_id = 0;
    std::unordered_set<std::string> refined_used;
    std::vector<EmittedChain> refined_chains;
    for (const auto& primary : primary_chains) {
        for (const auto& candidate : primary.anchors) {
            refined_used.insert(candidate_key(candidate));
        }
    }

    for (const auto& item : leftovers) {
        const auto& key = item.first;
        const auto& candidates = item.second;
        for (const auto& primary : primary_chains) {
            if (!same_chain_group(primary, key)) continue;
            std::vector<Candidate> subset;
            const auto ref_start = subtract_window(primary.ref_start,
                                                   parameters.refinement_window);
            const auto ref_end = primary.ref_end + parameters.refinement_window;
            const auto query_start = subtract_window(primary.query_start,
                                                     parameters.refinement_window);
            const auto query_end = primary.query_end + parameters.refinement_window;
            for (const auto& candidate : candidates) {
                if (candidate_in_box(candidate, ref_start, ref_end, query_start,
                                     query_end, refined_used)) {
                    subset.push_back(candidate);
                }
            }
            emit_refined_subset(std::move(subset), parameters, next_chain_id, result,
                                refined_used, refined_chains, chains,
                                chain_anchors, "around", std::to_string(primary.id),
                                ".", ".");
        }
    }

    std::sort(primary_chains.begin(), primary_chains.end(),
              [](const EmittedChain& left, const EmittedChain& right) {
                  const auto left_key = left.sample_a + '\x1f' + left.sample_b + '\x1f' +
                                        left.sequence_a + '\x1f' + left.sequence_b + '\x1f' +
                                        left.strand;
                  const auto right_key = right.sample_a + '\x1f' + right.sample_b + '\x1f' +
                                         right.sequence_a + '\x1f' + right.sequence_b + '\x1f' +
                                         right.strand;
                  if (left_key != right_key) return left_key < right_key;
                  return left.ref_start < right.ref_start;
              });
    for (std::size_t index = 1; index < primary_chains.size(); ++index) {
        const auto& left = primary_chains[index - 1];
        const auto& right = primary_chains[index];
        Candidate left_group;
        left_group.sample_a = left.sample_a;
        left_group.sample_b = left.sample_b;
        left_group.sequence_a = left.sequence_a;
        left_group.sequence_b = left.sequence_b;
        left_group.strand = left.strand;
        const auto key = group_key(left_group);
        if (!same_chain_group(right, key)) continue;
        if (left.ref_end > right.ref_start) continue;

        std::uint64_t query_start = 0;
        std::uint64_t query_end = 0;
        if (left.strand == '+') {
            if (left.query_end > right.query_start) continue;
            query_start = left.query_end;
            query_end = right.query_start;
        } else {
            if (right.query_end > left.query_start) continue;
            query_start = right.query_end;
            query_end = left.query_start;
        }

        const auto found = std::find_if(
            leftovers.begin(), leftovers.end(),
            [&](const auto& item) { return item.first == key; });
        if (found == leftovers.end()) continue;
        std::vector<Candidate> subset;
        for (const auto& candidate : found->second) {
            if (candidate_in_box(candidate, left.ref_end, right.ref_start,
                                 query_start, query_end, refined_used)) {
                subset.push_back(candidate);
            }
        }
        emit_refined_subset(std::move(subset), parameters, next_chain_id, result,
                            refined_used, refined_chains, chains, chain_anchors,
                            "between", ".", std::to_string(left.id),
                            std::to_string(right.id));
    }
    return refined_chains;
}

}  // namespace chaining_detail
}  // namespace vallescope2
