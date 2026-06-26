#include "vallescope2/correspondence/edit_distance.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void expect_distance(const std::string& name,
                     const std::vector<std::uint32_t>& left,
                     const std::vector<std::uint32_t>& right,
                     const std::uint32_t limit,
                     const std::uint32_t expected) {
    const auto observed =
        vallescope2::bounded_token_edit_distance(left, right, limit);
    if (observed != expected) {
        std::cerr << name << ": expected " << expected << ", observed "
                  << observed << '\n';
        std::exit(1);
    }
}

void expect_strand_choice() {
    const std::vector<std::uint32_t> query{1, 2, 3, 4};
    const std::vector<std::uint32_t> target{4, 3, 2, 1};
    auto reversed_query = query;
    std::reverse(reversed_query.begin(), reversed_query.end());

    const auto forward =
        vallescope2::bounded_token_edit_distance(query, target, 4);
    const auto reverse =
        vallescope2::bounded_token_edit_distance(reversed_query, target, 4);
    if (!(reverse == 0 && forward > reverse)) {
        std::cerr << "strand choice: expected reverse ED 0 and smaller than "
                     "forward ED; observed forward="
                  << forward << ", reverse=" << reverse << '\n';
        std::exit(1);
    }
}

void expect_anchor_distance_pair_indel_resynchronizes() {
    // A structural context is an alternating A/D token stream.  If one
    // anchor-distance pair is inserted, the edit distance should be 2 and the
    // suffix should realign.  A substitution-only interpretation would make
    // most downstream tokens mismatch.
    const std::vector<std::uint32_t> query{
        101, 201, 102, 202, 103, 203, 104, 204, 105, 205};
    const std::vector<std::uint32_t> target{
        101, 201, 999, 299, 102, 202, 103, 203, 104, 204, 105, 205};
    expect_distance("A-D pair insertion resynchronizes", query, target, 12, 2);
    expect_distance("A-D pair insertion over tight limit", query, target, 1, 2);
}

void expect_anchor_distance_pair_deletion_resynchronizes() {
    const std::vector<std::uint32_t> query{
        101, 201, 999, 299, 102, 202, 103, 203, 104, 204, 105, 205};
    const std::vector<std::uint32_t> target{
        101, 201, 102, 202, 103, 203, 104, 204, 105, 205};
    expect_distance("A-D pair deletion resynchronizes", query, target, 12, 2);
}

void expect_long_shift_resynchronizes() {
    std::vector<std::uint32_t> query;
    std::vector<std::uint32_t> target;
    for (std::uint32_t i = 0; i < 30; ++i) {
        query.push_back(1000 + i);  // A-like token
        query.push_back(2000 + i);  // D-like token
    }
    target = query;
    target.insert(target.begin() + 20, {9001, 9002, 9003, 9004});
    expect_distance("long suffix after insertion resynchronizes",
                    query, target, 10, 4);
}

void expect_repeated_prefix_pair_insertion_resynchronizes() {
    // Repeated A1/D1 units can make a naive visual alignment ambiguous.
    // The optimal edit distance should still treat the extra two A1/D1 pairs
    // as four inserted tokens and then resynchronize at A2/D2.
    const std::vector<std::uint32_t> query{
        101, 201, 102, 202, 103, 203, 104, 204, 105, 205};
    const std::vector<std::uint32_t> target{
        101, 201, 101, 201, 101, 201, 102, 202,
        103, 203, 104, 204, 105, 205};
    expect_distance("repeated prefix A-D pair insertion resynchronizes",
                    query, target, 10, 4);
}

void expect_tandem_duplication_resynchronizes() {
    // A local tandem duplication of A3/D3/A4/D4 should cost four inserted
    // tokens, not force all downstream A/D tokens into substitutions.
    const std::vector<std::uint32_t> query{
        101, 201, 102, 202, 103, 203, 104, 204, 105, 205};
    const std::vector<std::uint32_t> target{
        101, 201, 102, 202, 103, 203, 104, 204,
        103, 203, 104, 204, 105, 205};
    expect_distance("tandem duplication resynchronizes", query, target, 10, 4);
}

void expect_phase_disrupted_context_has_higher_distance() {
    // This intentionally includes an odd token insertion, creating adjacent
    // A-like tokens.  Since the A/D phase is disrupted, the edit distance
    // should be higher than the clean A-D pair insertion case but still
    // finite and bounded.
    const std::vector<std::uint32_t> query{
        101, 201, 102, 202, 103, 203, 104, 204, 105, 205};
    const std::vector<std::uint32_t> target{
        101, 201, 102, 202, 103, 203, 104,
        103, 203, 104, 204, 105, 205};
    expect_distance("phase-disrupted local duplication", query, target, 10, 3);
}

void expect_repetitive_ambiguous_indel_still_finds_minimum() {
    // A more repetitive case where several A/D pairs are identical.  There are
    // many equally good alignments, but the minimum distance should still be
    // the number of extra tokens.
    const std::vector<std::uint32_t> query{
        101, 201, 101, 201, 102, 202, 102, 202, 103, 203};
    const std::vector<std::uint32_t> target{
        101, 201, 101, 201, 101, 201,
        102, 202, 102, 202, 103, 203};
    expect_distance("ambiguous repetitive insertion finds minimum",
                    query, target, 10, 2);
}

}  // namespace

int main() {
    expect_distance("identical", {1, 2, 3}, {1, 2, 3}, 3, 0);
    expect_distance("substitution", {1, 2, 3}, {1, 9, 3}, 3, 1);
    expect_distance("insertion", {1, 2, 3}, {1, 2, 8, 3}, 3, 1);
    expect_distance("deletion", {1, 2, 8, 3}, {1, 2, 3}, 3, 1);
    expect_distance("mixed edits", {1, 2, 3, 4}, {1, 9, 3, 8, 4}, 5, 2);
    expect_distance("empty left", {}, {1, 2}, 2, 2);
    expect_distance("empty right", {1, 2}, {}, 2, 2);

    expect_distance("length difference over limit", {1, 2, 3}, {1}, 1, 2);
    expect_distance("distance over limit", {1, 2, 3}, {9, 8, 7}, 2, 3);
    expect_distance("distance equals limit", {1, 2, 3}, {9, 8, 3}, 2, 2);

    expect_strand_choice();
    expect_anchor_distance_pair_indel_resynchronizes();
    expect_anchor_distance_pair_deletion_resynchronizes();
    expect_long_shift_resynchronizes();
    expect_repeated_prefix_pair_insertion_resynchronizes();
    expect_tandem_duplication_resynchronizes();
    expect_phase_disrupted_context_has_higher_distance();
    expect_repetitive_ambiguous_indel_still_finds_minimum();
    return 0;
}
