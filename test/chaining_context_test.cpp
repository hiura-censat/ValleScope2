#include "chaining_internal.hpp"

#include <chrono>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using vallescope2::chaining_detail::EmittedChain;

EmittedChain make_chain(const std::uint64_t id,
                        const std::uint64_t ref_start,
                        const std::uint64_t ref_end,
                        const std::uint64_t query_start,
                        const std::uint64_t query_end,
                        const double score) {
    EmittedChain chain;
    chain.id = id;
    chain.sample_a = "sampleA";
    chain.sample_b = "sampleB";
    chain.sequence_a = "sampleA";
    chain.sequence_b = "sampleB";
    chain.strand = '+';
    chain.ref_start = ref_start;
    chain.ref_end = ref_end;
    chain.query_start = query_start;
    chain.query_end = query_end;
    chain.score = score;
    return chain;
}

void require(const bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    const auto directory = std::filesystem::temp_directory_path() /
        ("vallescope2-chain-context-test-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        std::filesystem::create_directories(directory);
        const auto left = make_chain(10, 0, 100, 0, 100, 100.0);
        const auto right = make_chain(11, 200, 300, 1000, 1100, 100.0);

        std::uint64_t filtered = 0;
        auto rejected = vallescope2::chaining_detail::
            filter_context_conflicting_chains(
                {left, right, make_chain(12, 1000, 1100, 110, 200, 10.0)},
                directory / "rejected.tsv", filtered);
        require(filtered == 1 && rejected.size() == 2,
                "orphan-gap cross-copy was not rejected");

        auto rescued = vallescope2::chaining_detail::
            filter_context_conflicting_chains(
                {left, right, make_chain(13, 1200, 1300, 110, 990, 10.0)},
                directory / "rescued.tsv", filtered);
        require(filtered == 0 && rescued.size() == 3,
                "boundary-supported structural candidate was rejected");

        std::filesystem::remove_all(directory);
        return 0;
    } catch (...) {
        std::filesystem::remove_all(directory);
        throw;
    }
}
