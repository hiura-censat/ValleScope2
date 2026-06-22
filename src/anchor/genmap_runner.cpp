#include "vallescope2/anchor/genmap_runner.hpp"

#include "vallescope2/common/process.hpp"
#include "vallescope2/fasta.hpp"

#include <fstream>
#include <stdexcept>

namespace vallescope2 {

GenmapResult run_genmap(
    const std::vector<std::filesystem::path>& input_paths,
    const PreparationResult& prepared,
    const std::filesystem::path& work_directory,
    const GenmapParameters& parameters) {
    GenmapResult result;
    result.index_directory = work_directory / "genmap_index";
    const auto output_directory = work_directory / "genmap_output";
    result.output_prefix = output_directory / "genout";
    result.frequency_raw = result.output_prefix.string() + ".freq16";
    result.log = work_directory / "genmap.log";

    std::filesystem::remove(work_directory / "genmap_input.fa");
    std::filesystem::remove_all(result.index_directory);
    std::filesystem::remove_all(output_directory);
    std::filesystem::remove(result.log);
    std::filesystem::create_directories(output_directory);

    if (prepared.mode == ComparisonMode::multi_self) {
        result.input_fasta = prepared.target.fasta;
        result.catalog = prepared.catalog;
    } else {
        result.input_fasta = work_directory / "genmap_input.fa";
        const std::vector<std::string> roles =
            input_paths.size() == 1
                ? std::vector<std::string>{"self"}
                : std::vector<std::string>{"target", "query"};
        result.catalog = combine_fastas(input_paths, roles, result.input_fasta);
    }

    const std::string executable = parameters.executable.string();
    result.version = capture_process({executable, "--version"});
    run_process({executable, "index", "-F", result.input_fasta.string(),
                 "-I", result.index_directory.string()},
                result.log);

    const std::vector<std::string> map_command{
        executable,
        "map",
        "-K",
        std::to_string(parameters.kmer_length),
        "-E",
        std::to_string(parameters.errors),
        "-I",
        result.index_directory.string(),
        "-O",
        result.output_prefix.string(),
        "-r",
        "-fl"};
    result.command = format_command(map_command);
    run_process(map_command, result.log);

    if (!std::filesystem::is_regular_file(result.frequency_raw)) {
        throw std::runtime_error("GenMap did not create raw frequency output: " +
                                 result.frequency_raw.string());
    }
    return result;
}

}  // namespace vallescope2
