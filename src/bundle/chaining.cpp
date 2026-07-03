#include "chaining_internal.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace vallescope2 {
namespace {
using namespace chaining_detail;

void validate_parameters(const ChainingParameters& parameters) {
    if (parameters.predecessor_count == 0) {
        throw std::runtime_error("chain predecessor count must be greater than zero");
    }
    if (parameters.max_chain_gap == 0) {
        throw std::runtime_error("maximum chain gap must be greater than zero");
    }
    if (parameters.gap_unit == 0) {
        throw std::runtime_error("gap unit must be greater than zero");
    }
    if (parameters.chain_max_gap_ratio < 1.0) {
        throw std::runtime_error("chain max gap ratio must be at least 1");
    }
    if (parameters.min_chain_anchors == 0) {
        throw std::runtime_error("minimum chain anchor count must be greater than zero");
    }
    if (parameters.chain_trim_overlap < 0.0 ||
        parameters.chain_trim_overlap > 1.0) {
        throw std::runtime_error("chain trim overlap must be between 0 and 1");
    }
    if (parameters.refinement_min_chain_anchors == 0) {
        throw std::runtime_error("minimum refinement chain anchor count must be greater than zero");
    }
}

std::vector<std::pair<std::string, std::vector<Candidate>>> group_candidates(
    std::vector<Candidate> candidates) {
    std::unordered_map<std::string, std::vector<Candidate>> groups;
    for (auto& candidate : candidates) {
        groups[group_key(candidate)].push_back(std::move(candidate));
    }

    std::vector<std::pair<std::string, std::vector<Candidate>>> ordered_groups(
        std::make_move_iterator(groups.begin()),
        std::make_move_iterator(groups.end()));
    std::sort(ordered_groups.begin(), ordered_groups.end(),
              [](const auto& left, const auto& right) {
                  return left.first < right.first;
              });
    return ordered_groups;
}

void add_assignment_candidates(
    std::vector<std::pair<std::string, std::vector<Candidate>>>& ordered_groups,
    std::vector<Candidate> assignment_candidates) {
    std::unordered_map<std::string, std::size_t> group_index;
    for (std::size_t index = 0; index < ordered_groups.size(); ++index) {
        group_index.emplace(ordered_groups[index].first, index);
    }
    for (auto& candidate : assignment_candidates) {
        const auto key = group_key(candidate);
        const auto found = group_index.find(key);
        if (found == group_index.end()) {
            group_index.emplace(key, ordered_groups.size());
            ordered_groups.emplace_back(key, std::vector<Candidate>{});
            ordered_groups.back().second.push_back(std::move(candidate));
        } else {
            ordered_groups[found->second].second.push_back(std::move(candidate));
        }
    }
    std::sort(ordered_groups.begin(), ordered_groups.end(),
              [](const auto& left, const auto& right) {
                  return left.first < right.first;
              });
}

}  // namespace

ChainingResult build_anchor_chains(
    const std::filesystem::path& correspondences,
    const std::filesystem::path& assignments,
    const std::filesystem::path& grouped_anchors,
    const std::filesystem::path& chain_output,
    const std::filesystem::path& chain_anchor_output,
    const std::filesystem::path& metadata_output,
    const std::filesystem::path& refined_chain_output,
    const std::filesystem::path& refined_chain_anchor_output,
    const std::filesystem::path& refined_metadata_output,
    const ChainingParameters& parameters) {
    validate_parameters(parameters);

    const auto positions = load_anchor_positions(grouped_anchors);
    auto candidates = load_candidates(correspondences, positions);
    auto assignment_candidates = load_assignment_candidates(assignments, positions);

    ChainingResult result;
    result.candidate_count = candidates.size();

    std::ofstream chains(chain_output);
    if (!chains) throw std::runtime_error("cannot create chain output");
    std::ofstream chain_anchors(chain_anchor_output);
    if (!chain_anchors) throw std::runtime_error("cannot create chain anchor output");
    chains << "chain_id\tsample_a\tsample_b\tsequence_a\tsequence_b"
              "\tassign_strand\tn_candidates\tboth_anchor_count\tboth_rate"
              "\tchain_score"
              "\tref_start\tref_end\tquery_start\tquery_end\n";
    chain_anchors << "chain_id\trank\tsample_a\tsample_b\tanchor_a\tanchor_b"
                     "\tref_center\tquery_center\tassign_strand"
                     "\tcandidate_score\tsupport_direction\n";

    auto ordered_groups = group_candidates(std::move(candidates));

    std::vector<EmittedChain> primary_chains;
    for (auto& item : ordered_groups) {
        auto& group = item.second;
        std::vector<Candidate> refinement_candidates;
        while (!group.empty()) {
            double chain_score = 0.0;
            const auto chain = chain_group(group, parameters, chain_score);
            if (chain.empty()) break;
            if (chain_score < parameters.min_chain_score) break;

            std::vector<bool> used(group.size(), false);
            for (const auto index : chain) used[index] = true;
            if (chain.size() < parameters.min_chain_anchors) {
                std::vector<Candidate> remaining;
                remaining.reserve(group.size() - chain.size());
                for (std::size_t index = 0; index < group.size(); ++index) {
                    if (used[index])
                        refinement_candidates.push_back(std::move(group[index]));
                    else
                        remaining.push_back(std::move(group[index]));
                }
                group = std::move(remaining);
                continue;
            }

            auto emitted =
                make_emitted_chain(primary_chains.size(), group, chain, chain_score);
            if (emitted.both_anchor_count < parameters.min_chain_both_anchors) {
                std::vector<Candidate> remaining;
                remaining.reserve(group.size() - chain.size());
                for (std::size_t index = 0; index < group.size(); ++index) {
                    if (used[index])
                        refinement_candidates.push_back(std::move(group[index]));
                    else
                        remaining.push_back(std::move(group[index]));
                }
                group = std::move(remaining);
                continue;
            }
            primary_chains.push_back(std::move(emitted));

            std::vector<Candidate> remaining;
            remaining.reserve(group.size() - chain.size());
            for (std::size_t index = 0; index < group.size(); ++index) {
                if (!used[index]) remaining.push_back(std::move(group[index]));
            }
            group = std::move(remaining);
        }
        group.insert(group.end(),
                     std::make_move_iterator(refinement_candidates.begin()),
                     std::make_move_iterator(refinement_candidates.end()));
    }

    result.raw_chain_count = primary_chains.size();
    add_assignment_candidates(ordered_groups, std::move(assignment_candidates));
    auto refined_chains = refine_chains(
        ordered_groups, primary_chains, refined_chain_output,
        refined_chain_anchor_output, parameters, result);

    std::vector<EmittedChain> emitted_chains = std::move(primary_chains);
    emitted_chains.insert(emitted_chains.end(),
                          std::make_move_iterator(refined_chains.begin()),
                          std::make_move_iterator(refined_chains.end()));
    const auto pre_trim_chain_count = emitted_chains.size();
    emitted_chains = trim_overlapping_chains(std::move(emitted_chains), parameters);
    result.trimmed_chain_count =
        pre_trim_chain_count > emitted_chains.size()
            ? pre_trim_chain_count - emitted_chains.size()
            : 0;
    result.chain_count = emitted_chains.size();
    for (const auto& emitted : emitted_chains) {
        result.chain_anchor_count += emitted.anchors.size();
        write_chain(chains, chain_anchors, emitted);
    }

    write_metadata(metadata_output, parameters, result);
    write_metadata(refined_metadata_output, parameters, result);
    return result;
}

}  // namespace vallescope2
