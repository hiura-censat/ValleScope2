#include "chaining_internal.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
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

std::string read_text(const std::filesystem::path& path) {
    std::ifstream input(path);
    return {std::istreambuf_iterator<char>(input),
            std::istreambuf_iterator<char>()};
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

        const auto excursion_path = directory / "excursion.tsv";
        auto excursion_filtered = vallescope2::chaining_detail::
            filter_context_conflicting_chains(
                {make_chain(20, 0, 100, 34000, 34100, 100.0),
                 make_chain(21, 14000, 15000, 35400, 36400, 100.0),
                 make_chain(22, 24000, 25000, 56000, 57000, 100.0)},
                excursion_path, filtered);
        require(filtered == 1 && excursion_filtered.size() == 2,
                "short diagonal excursion was not removed from the skeleton");
        require(read_text(excursion_path).find("off_skeleton_excursion") !=
                    std::string::npos,
                "diagonal excursion was not recorded");

        auto deletion_retained = vallescope2::chaining_detail::
            filter_context_conflicting_chains(
                {make_chain(30, 0, 100, 0, 100, 100.0),
                 make_chain(31, 9000, 10000, 1000, 2000, 100.0),
                 make_chain(32, 11000, 12000, 3000, 4000, 100.0)},
                directory / "deletion.tsv", filtered);
        require(filtered == 0 && deletion_retained.size() == 3,
                "persistent one-way deletion shift was rejected");

        std::filesystem::remove_all(directory);
        return 0;
    } catch (...) {
        std::filesystem::remove_all(directory);
        throw;
    }
}
