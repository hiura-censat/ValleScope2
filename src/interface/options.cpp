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
           << "  -md, --min-center-distance INT  Minimum center distance [50]\n"
           << "  -td, --target-density INT  Target anchors per Mb [15000]\n"
           << "  -db, --distance-bin-size INT  Distance token bin size [10]\n"
           << "  -cr, --context-radius-tokens INT  Context radius in tokens [50]\n"
           << "  --min-shared-tmus INT  Shared exact tMUS fast-path threshold [3]\n"
           << "  -mcs, --min-candidate-score INT  Minimum candidate score [10]\n"
           << "  -pm, --primary-margin INT  Primary assignment score margin [5]\n"
           << "  --pair-merge-mode MODE  reciprocal or union [union]\n"
           << "  --chain-predecessors INT  Chaining predecessor search limit [50]\n"
           << "  --max-chain-gap INT  Maximum ref/query gap between chained anchors [300]\n"
           << "  --chain-max-gap-ratio FLOAT  Maximum max(dr,dq)/min(dr,dq) [1.05]\n"
           << "  --gap-cost-model MODE  absolute or relative [absolute]\n"
           << "  --gap-weight FLOAT  Linear chaining gap weight [0.002]\n"
           << "  --min-chain-anchors INT  Minimum anchors per emitted chain [10]\n"
           << "  --min-chain-both-anchors INT  Minimum both-supported anchors per emitted chain [1]\n"
           << "  --min-chain-score FLOAT  Minimum score per emitted chain [0]\n"
           << "  --chain-trim-overlap FLOAT  Ref/query overlap for trimming lower-score chains/bundles [0.01]\n"
           << "  --chain-reciprocal-weight FLOAT  Both-rate weight in adjusted chain score [0.8]\n"
           << "  --chain-multimap-overlap FLOAT  One-axis overlap for competing-chain suppression [0.5]\n"
           << "  --refinement-window INT  Refinement interval extension in bp [100000]\n"
           << "  --refinement-min-chain-anchors INT  Minimum anchors per refined chain [20]\n"
           << "  --base-align  Write base-level bundle PAF to stdout with cg:Z [default]\n"
           << "  --no-base-align  Stop after chaining and do not write final PAF\n"
           << "  --chain-extension  Extend adjacent chains using correspondence DP before gap patching [default]\n"
           << "  --no-chain-extension  Disable chain extension before gap patching\n"
           << "  --max-chain-extension-bp INT  Maximum bundle gap considered for chain extension [70000]\n"
           << "  --min-chain-extension-anchors INT  Minimum anchors for partial chain extension [5]\n"
           << "  --min-chain-extension-score FLOAT  Minimum DP score for partial chain extension [100]\n"
           << "  --min-copy-support-anchors INT  both-supported anchors needed to classify a residual as dup branch [1]\n"
           << "  --max-bundle-align-bp INT  Skip bundle alignments longer than this [1000000]\n"
           << "  --max-patch-gap-bp INT  Maximum adjacent bundle gap for patching [70000]\n"
           << "  --patch-flank-bp INT  Trusted flank size for patch interval WFA [500]\n"
           << "  --max-wfa-memory-gb INT  WFA2 memory limit in GB [64]\n"
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

GapCostModel parse_gap_cost_model(const int argc, char* argv[], int& index,
                                  const std::string& option) {
    if (++index >= argc) throw std::runtime_error(option + " requires a mode");
    const std::string value = argv[index];
    if (value == "absolute") return GapCostModel::absolute;
    if (value == "relative") return GapCostModel::relative;
    throw std::runtime_error("invalid value for " + option + ": " + value);
}

double parse_double(const int argc, char* argv[], int& index,
                    const std::string& option) {
    if (++index >= argc) throw std::runtime_error(option + " requires a value");
    const std::string value = argv[index];
    try {
        std::size_t parsed = 0;
        const double number = std::stod(value, &parsed);
        if (value.empty() || parsed != value.size()) {
            throw std::out_of_range("invalid number");
        }
        return number;
    } catch (const std::exception&) {
        throw std::runtime_error("invalid value for " + option + ": " + value);
    }
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
        else if (argument == "--base-align") options.base_align = true;
        else if (argument == "--no-base-align") options.base_align = false;
        else if (argument == "--chain-extension") options.chain_extension = true;
        else if (argument == "--no-chain-extension") options.chain_extension = false;
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
        else if (argument == "--min-shared-tmus")
            options.assignment.min_shared_tmus =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "-mcs" || argument == "--min-candidate-score")
            options.assignment.min_candidate_score =
                parse_signed(argc, argv, index, argument);
        else if (argument == "-pm" || argument == "--primary-margin")
            options.assignment.primary_margin = parse_unsigned(argc, argv, index, argument);
        else if (argument == "--pair-merge-mode")
            options.assignment.pair_merge_mode =
                parse_pair_merge_mode(argc, argv, index, argument);
        else if (argument == "--chain-predecessors")
            options.chaining.predecessor_count =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--max-chain-gap")
            options.chaining.max_chain_gap =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--chain-max-gap-ratio")
            options.chaining.chain_max_gap_ratio =
                parse_double(argc, argv, index, argument);
        else if (argument == "--gap-cost-model")
            options.chaining.gap_cost_model =
                parse_gap_cost_model(argc, argv, index, argument);
        else if (argument == "--gap-weight")
            options.chaining.gap_weight =
                parse_double(argc, argv, index, argument);
        else if (argument == "--min-chain-anchors")
            options.chaining.min_chain_anchors =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--min-chain-both-anchors")
            options.chaining.min_chain_both_anchors =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--min-chain-score")
            options.chaining.min_chain_score =
                parse_double(argc, argv, index, argument);
        else if (argument == "--chain-trim-overlap")
            options.chaining.chain_trim_overlap =
                parse_double(argc, argv, index, argument);
        else if (argument == "--chain-reciprocal-weight")
            options.chaining.reciprocal_score_weight =
                parse_double(argc, argv, index, argument);
        else if (argument == "--chain-multimap-overlap")
            options.chaining.multimap_overlap =
                parse_double(argc, argv, index, argument);
        else if (argument == "--refinement-window")
            options.chaining.refinement_window =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--refinement-min-chain-anchors")
            options.chaining.refinement_min_chain_anchors =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--max-chain-extension-bp")
            options.max_chain_extension_bp =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--min-chain-extension-anchors")
            options.min_chain_extension_anchors =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--min-chain-extension-score")
            options.min_chain_extension_score =
                parse_double(argc, argv, index, argument);
        else if (argument == "--min-copy-support-anchors")
            options.min_copy_support_anchors =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--max-bundle-align-bp")
            options.max_bundle_align_bp =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--max-patch-gap-bp")
            options.max_patch_gap_bp =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--patch-flank-bp")
            options.patch_flank_bp =
                parse_unsigned(argc, argv, index, argument);
        else if (argument == "--max-wfa-memory-gb")
            options.max_wfa_memory_gb =
                parse_unsigned(argc, argv, index, argument);
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
    options.chaining.gap_unit = 10;
    if (options.chaining.predecessor_count == 0)
        throw std::runtime_error("chain predecessor count must be greater than zero");
    if (options.chaining.max_chain_gap == 0)
        throw std::runtime_error("maximum chain gap must be greater than zero");
    if (options.chaining.gap_weight < 0.0)
        throw std::runtime_error("gap weight must be non-negative");
    if (options.chaining.chain_max_gap_ratio < 1.0)
        throw std::runtime_error("chain max gap ratio must be at least 1");
    if (options.chaining.min_chain_anchors == 0)
        throw std::runtime_error("minimum chain anchor count must be greater than zero");
    if (options.chaining.chain_trim_overlap < 0.0 ||
        options.chaining.chain_trim_overlap > 1.0)
        throw std::runtime_error("chain trim overlap must be between 0 and 1");
    if (options.chaining.reciprocal_score_weight < 0.0 ||
        options.chaining.reciprocal_score_weight > 1.0)
        throw std::runtime_error(
            "chain reciprocal weight must be between 0 and 1");
    if (options.chaining.multimap_overlap < 0.0 ||
        options.chaining.multimap_overlap > 1.0)
        throw std::runtime_error(
            "chain multimap overlap must be between 0 and 1");
    if (options.chaining.refinement_min_chain_anchors == 0)
        throw std::runtime_error("minimum refinement chain anchor count must be greater than zero");
    if (options.max_bundle_align_bp == 0)
        throw std::runtime_error("maximum bundle alignment length must be greater than zero");
    if (options.max_chain_extension_bp == 0)
        throw std::runtime_error("maximum chain extension length must be greater than zero");
    if (options.min_chain_extension_anchors == 0)
        throw std::runtime_error("minimum chain extension anchor count must be greater than zero");
    if (options.patch_flank_bp == 0)
        throw std::runtime_error("patch flank size must be greater than zero");
    if (options.max_wfa_memory_gb == 0)
        throw std::runtime_error("maximum WFA memory must be greater than zero");
    if (options.context_radius_tokens >
        (std::numeric_limits<std::uint32_t>::max() - 1) / 2) {
        throw std::runtime_error("context radius is too large");
    }
    if (options.anchor_length < options.genmap.kmer_length)
        throw std::runtime_error("anchor length must be greater than or equal to k-mer length");
    return options;
}

}  // namespace vallescope2
