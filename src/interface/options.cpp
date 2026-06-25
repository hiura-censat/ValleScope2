#include "vallescope2/interface/options.hpp"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#ifndef VALLESCOPE2_VERSION
#define VALLESCOPE2_VERSION "unknown"
#endif

namespace vallescope2 {
namespace {

void print_usage(std::ostream& output) {
    output << "Usage: vallescope2 [options] <input.fa>\n"
           << "       vallescope2 [options] <target.fa> <query.fa>\n"
           << "       vallescope2 [options] --combine <input1.fa> <input2.fa> [...]\n\n"
           << "Options:\n"
           << "  --combine  Combine all input FASTAs for all-vs-all self comparison\n"
           << "  -K, --kmer-length INT  GenMap k-mer length [6]\n"
           << "  -E, --kmer-errors INT  GenMap errors [0]\n"
           << "  -al, --anchor-length INT  Anchor/window length [50]\n"
           << "  -md, --min-center-distance INT  Minimum center distance [100]\n"
           << "  -td, --target-density INT  Target anchors per Mb [6000]\n"
           << "  -db, --distance-bin-size INT  Distance token bin size [10]\n"
           << "  -cr, --context-radius-tokens INT  Context radius in tokens [50]\n"
           << "  -bt, --beta-tolerance INT  Legacy context ED tolerance [30]\n"
           << "  -mcs, --min-candidate-score INT  Minimum candidate score [-10]\n"
           << "  -pm, --primary-margin INT  Primary assignment score margin [5]\n"
           << "  --pair-merge-mode MODE  reciprocal or union [union]\n"
           << "  -t, --threads INT  Number of GenMap threads [20]\n"
           << "  --genmap PATH  GenMap executable [genmap]\n"
           << "  --debug  Keep intermediate FASTA, index, and GenMap files\n"
           << "  --dump-window-scores  Write window_scores.tsv (very large)\n"
           << "  -h, --help  Show this help message\n"
           << "  -v, --version  Show version information\n";
}

std::uint32_t parse_unsigned(const int argc, char* argv[], int& index,
                             const std::string& option) {
    if (++index >= argc) throw std::runtime_error(option + " requires an integer value");
    const std::string value = argv[index];
    try {
        std::size_t parsed = 0;
        const unsigned long number = std::stoul(value, &parsed);
        if (value.empty() || value.front() == '-' || parsed != value.size() ||
            number > std::numeric_limits<std::uint32_t>::max()) {
            throw std::out_of_range("invalid integer");
        }
        return static_cast<std::uint32_t>(number);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
}

std::int64_t parse_signed(const int argc, char* argv[], int& index,
                          const std::string& option) {
    if (++index >= argc) throw std::runtime_error(option + " requires an integer value");
    const std::string value = argv[index];
    try {
        std::size_t parsed = 0;
        const long long number = std::stoll(value, &parsed);
        if (value.empty() || parsed != value.size()) {
            throw std::out_of_range("invalid integer");
        }
        return static_cast<std::int64_t>(number);
    } catch (const std::exception&) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
}

PairMergeMode parse_pair_merge_mode(const int argc, char* argv[], int& index,
                                    const std::string& option) {
    if (++index >= argc) throw std::runtime_error(option + " requires a mode");
    const std::string value = argv[index];
    if (value == "reciprocal") return PairMergeMode::reciprocal;
    if (value == "union") return PairMergeMode::union_mode;
    throw std::runtime_error("invalid value for " + option + ": " + value);
}

}  // namespace

ProgramOptions parse_arguments(const int argc, char* argv[]) {
    ProgramOptions options;
    const auto local_genmap = std::filesystem::current_path() / ".tools/genmap/bin/genmap";
    if (std::filesystem::is_regular_file(local_genmap)) options.genmap.executable = local_genmap;

    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--combine") options.combine = true;
        else if (argument == "--debug") options.debug = true;
        else if (argument == "--dump-window-scores") options.dump_window_scores = true;
        else if (argument == "-K" || argument == "--kmer-length")
            options.genmap.kmer_length = parse_unsigned(argc, argv, index, argument);
        else if (argument == "-E" || argument == "--kmer-errors")
            options.genmap.errors = parse_unsigned(argc, argv, index, argument);
        else if (argument == "-al" || argument == "--anchor-length")
            options.anchor_length = parse_unsigned(argc, argv, index, argument);
        else if (argument == "-md" || argument == "--min-center-distance")
            options.min_center_distance = parse_unsigned(argc, argv, index, argument);
        else if (argument == "-td" || argument == "--target-density")
            options.target_density = parse_unsigned(argc, argv, index, argument);
        else if (argument == "-db" || argument == "--distance-bin-size")
            options.distance_bin_size = parse_unsigned(argc, argv, index, argument);
        else if (argument == "-cr" || argument == "--context-radius-tokens")
            options.context_radius_tokens = parse_unsigned(argc, argv, index, argument);
        else if (argument == "-bt" || argument == "--beta-tolerance")
            options.assignment.beta_tolerance = parse_unsigned(argc, argv, index, argument);
        else if (argument == "-mcs" || argument == "--min-candidate-score")
            options.assignment.min_candidate_score =
                parse_signed(argc, argv, index, argument);
        else if (argument == "-pm" || argument == "--primary-margin")
            options.assignment.primary_margin = parse_unsigned(argc, argv, index, argument);
        else if (argument == "--pair-merge-mode")
            options.assignment.pair_merge_mode =
                parse_pair_merge_mode(argc, argv, index, argument);
        else if (argument == "-t" || argument == "--threads")
            options.genmap.threads = parse_unsigned(argc, argv, index, argument);
        else if (argument == "--genmap") {
            if (++index >= argc) throw std::runtime_error("--genmap requires a path");
            options.genmap.executable = argv[index];
        } else if (argument == "-h" || argument == "--help") {
            print_usage(std::cout);
            std::exit(0);
        } else if (argument == "-v" || argument == "--version") {
            std::cout << "vallescope2 " << VALLESCOPE2_VERSION << '\n';
            std::exit(0);
        } else if (!argument.empty() && argument.front() == '-') {
            throw std::runtime_error("unknown option: " + argument);
        } else options.input_paths.emplace_back(argument);
    }
    if (options.input_paths.empty()) throw std::runtime_error("at least one FASTA file is required");
    if (options.combine && options.input_paths.size() < 2)
        throw std::runtime_error("--combine requires at least two FASTA files");
    if (!options.combine && options.input_paths.size() > 2)
        throw std::runtime_error("more than two FASTA files require the --combine option");
    if (options.genmap.kmer_length == 0)
        throw std::runtime_error("GenMap k-mer length must be greater than zero");
    if (options.genmap.threads == 0)
        throw std::runtime_error("thread count must be greater than zero");
    if (options.distance_bin_size == 0)
        throw std::runtime_error("distance bin size must be greater than zero");
    if (options.context_radius_tokens >
        (std::numeric_limits<std::uint32_t>::max() - 1) / 2) {
        throw std::runtime_error("context radius is too large");
    }
    if (options.assignment.beta_tolerance ==
        std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("beta tolerance is too large");
    }
    if (options.anchor_length < options.genmap.kmer_length)
        throw std::runtime_error("anchor length must be greater than or equal to k-mer length");
    return options;
}

}  // namespace vallescope2
