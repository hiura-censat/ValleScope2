#include "assignment_internal.hpp"

#include "vallescope2/correspondence/edit_distance.hpp"

#include <algorithm>

namespace vallescope2 {

std::uint32_t count_shared_tmus(const std::vector<std::uint32_t>& left,
                                const std::vector<std::uint32_t>& right) {
    std::uint32_t count = 0;
    auto l = left.begin();
    auto r = right.begin();
    while (l != left.end() && r != right.end()) {
        if (*l == *r) {
            ++count;
            ++l;
            ++r;
        } else if (*l < *r) {
            ++l;
        } else {
            ++r;
        }
    }
    return count;
}

std::vector<std::uint32_t> reversed_tokens(
    const std::vector<std::uint32_t>& tokens) {
    return {tokens.rbegin(), tokens.rend()};
}

Candidate make_candidate(const Anchor& query,
                         const std::vector<std::uint32_t>& reverse_query,
                         const Anchor& target,
                         const std::uint64_t q_alpha_prime,
                         const std::uint64_t r_alpha_prime,
                         const std::uint32_t edit_distance_limit) {
    const auto forward_ed = bounded_token_edit_distance(
        query.forward_tokens, target.forward_tokens, edit_distance_limit);
    const auto reverse_ed = bounded_token_edit_distance(
        reverse_query, target.forward_tokens, edit_distance_limit);
    const bool use_reverse = reverse_ed < forward_ed;
    const auto edit_distance = use_reverse ? reverse_ed : forward_ed;

    Candidate candidate;
    candidate.target = &target;
    candidate.edit_distance = edit_distance;
    candidate.strand = use_reverse ? '-' : '+';
    candidate.q_alpha_prime = q_alpha_prime;
    candidate.r_alpha_prime = r_alpha_prime;
    candidate.score = static_cast<std::int64_t>(candidate.q_alpha_prime) +
                      static_cast<std::int64_t>(candidate.r_alpha_prime) -
                      static_cast<std::int64_t>(edit_distance);
    return candidate;
}

bool candidate_less(const Candidate& left, const Candidate& right) {
    if (left.score != right.score) return left.score > right.score;
    if (left.edit_distance != right.edit_distance) {
        return left.edit_distance < right.edit_distance;
    }
    return left.target->anchor_id < right.target->anchor_id;
}

}  // namespace vallescope2
