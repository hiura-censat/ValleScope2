#include "vallescope2/anchor/window_scoring.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace vallescope2 {
namespace {

constexpr std::uint16_t kMaximumFrequency = 65535;

const std::array<double, 65536>& log_frequency_table() {
    static const auto table = [] {
        std::array<double, 65536> values{};
        for (std::size_t i = 0; i < values.size(); ++i) {
            values[i] = std::log1p(static_cast<double>(i));
        }
        return values;
    }();
    return table;
}

}  // namespace

WindowScoringResult score_windows(
    const std::filesystem::path& frequency_raw,
    const SequenceCatalog& catalog,
    const WindowScoringParameters& parameters,
    const std::function<void(const WindowScore&)>& consume) {
    if (parameters.anchor_length < parameters.kmer_length) {
        throw std::runtime_error(
            "anchor length must be greater than or equal to k-mer length");
    }
    std::ifstream input(frequency_raw, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot open GenMap raw frequency file: " +
                                 frequency_raw.string());
    }

    const std::uint64_t values_per_window =
        parameters.anchor_length - parameters.kmer_length + 1;
    std::vector<double> ring(values_per_window);
    const auto& log_values = log_frequency_table();
    WindowScoringResult result;

    for (std::size_t sequence_id = 0;
         sequence_id < catalog.records().size(); ++sequence_id) {
        const auto& sequence = catalog.records()[sequence_id];
        std::size_t ring_size = 0;
        std::size_t ring_next = 0;
        long double running_sum = 0.0;
        const std::uint64_t valid_positions =
            sequence.length >= parameters.kmer_length
                ? sequence.length - parameters.kmer_length + 1
                : 0;

        for (std::uint64_t position = 0; position < sequence.length; ++position) {
            std::uint16_t frequency = 0;
            input.read(reinterpret_cast<char*>(&frequency), sizeof(frequency));
            if (!input) {
                throw std::runtime_error(
                    "GenMap raw frequency file ended unexpectedly");
            }
            if (position >= valid_positions) continue;
            if (frequency == 0) {
                ring_size = 0;
                ring_next = 0;
                running_sum = 0.0;
                ++result.gap_positions;
                continue;
            }
            if (frequency == kMaximumFrequency) {
                ++result.saturated_positions;
            }

            const double value = log_values[frequency];
            if (ring_size < ring.size()) {
                ring[ring_next] = value;
                running_sum += value;
                ++ring_size;
            } else {
                running_sum -= ring[ring_next];
                ring[ring_next] = value;
                running_sum += value;
            }
            ring_next = (ring_next + 1) % ring.size();

            if (ring_size == ring.size()) {
                const std::uint64_t start = position + 1 - values_per_window;
                consume({static_cast<std::uint32_t>(sequence_id),
                         start,
                         start + parameters.anchor_length,
                         static_cast<double>(running_sum / values_per_window)});
                ++result.window_count;
            }
        }
    }

    char extra = 0;
    if (input.read(&extra, 1)) {
        throw std::runtime_error("GenMap raw frequency file is larger than expected");
    }
    return result;
}

}  // namespace vallescope2
