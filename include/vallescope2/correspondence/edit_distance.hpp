#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <vector>

namespace vallescope2 {

inline std::uint32_t bounded_token_edit_distance(
    const std::vector<std::uint32_t>& left,
    const std::vector<std::uint32_t>& right,
    const std::uint32_t limit) {
    const std::uint32_t fail = limit + 1;
    const auto n = static_cast<std::int64_t>(left.size());
    const auto m = static_cast<std::int64_t>(right.size());
    if (std::llabs(n - m) > static_cast<long long>(limit)) return fail;
    if (left.empty()) return static_cast<std::uint32_t>(right.size());
    if (right.empty()) return static_cast<std::uint32_t>(left.size());

    const std::uint32_t inf = fail;
    std::vector<std::uint32_t> previous(right.size() + 1, inf);
    std::vector<std::uint32_t> current(right.size() + 1, inf);
    const std::uint32_t first_end =
        std::min<std::uint32_t>(static_cast<std::uint32_t>(right.size()), limit);
    for (std::uint32_t j = 0; j <= first_end; ++j) previous[j] = j;

    for (std::uint32_t i = 1; i <= left.size(); ++i) {
        std::fill(current.begin(), current.end(), inf);
        const std::uint32_t begin = i > limit ? i - limit : 1;
        const std::uint32_t end =
            std::min<std::uint32_t>(static_cast<std::uint32_t>(right.size()),
                                    i + limit);
        if (begin == 1) current[0] = i <= limit ? i : inf;
        std::uint32_t row_min = inf;
        for (std::uint32_t j = begin; j <= end; ++j) {
            const std::uint32_t substitution =
                previous[j - 1] + (left[i - 1] == right[j - 1] ? 0 : 1);
            const std::uint32_t deletion = previous[j] + 1;
            const std::uint32_t insertion = current[j - 1] + 1;
            current[j] = std::min({substitution, deletion, insertion, inf});
            row_min = std::min(row_min, current[j]);
        }
        if (row_min > limit) return fail;
        previous.swap(current);
    }
    return previous[right.size()] <= limit ? previous[right.size()] : fail;
}

}  // namespace vallescope2
