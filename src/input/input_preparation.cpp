#include "vallescope2/input_preparation.hpp"

#include "vallescope2/fasta.hpp"
#include "vallescope2/fasta_index.hpp"

#include <stdexcept>
#include <unordered_set>

namespace vallescope2 {
namespace {

void set_index_paths(PreparedFasta& prepared,
                     const std::filesystem::path& output_directory,
                     const std::string& basename) {
    prepared.index = output_directory / (basename + ".fai");
    prepared.gzi = output_directory / (basename + ".gzi");
}

void add_indexed_sequences(SequenceCatalog& catalog,
                           const PreparedFasta& prepared,
                           const std::string& sample,
                           const std::filesystem::path& source,
                           const std::string& role) {
    for (const auto& sequence :
         read_fasta_index(prepared.fasta, prepared.index, prepared.gzi)) {
        catalog.add({sequence.name,
                     sample,
                     "1",
                     sequence.name,
                     source,
                     sequence.length,
                     role});
    }
}

}  // namespace

PreparationResult prepare_inputs(
    const std::vector<std::filesystem::path>& input_paths,
    const std::filesystem::path& output_directory,
    const bool combine) {
    if (input_paths.empty()) {
        throw std::runtime_error("at least one FASTA file is required");
    }

    std::vector<std::string> samples;
    std::unordered_set<std::string> unique_samples;
    for (const auto& path : input_paths) {
        if (!std::filesystem::is_regular_file(path)) {
            throw std::runtime_error("input is not a regular file: " +
                                     path.string());
        }
        const std::string sample = sample_name_from_path(path);
        if (!unique_samples.insert(sample).second) {
            throw std::runtime_error("duplicate sample name: " + sample);
        }
        samples.push_back(sample);
    }

    std::filesystem::create_directories(output_directory);
    for (const std::string filename : {"target.fa", "target.fa.fai",
                                       "target.fa.gzi", "query.fa",
                                       "query.fa.fai", "query.fa.gzi",
                                       "combined.fa", "combined.fa.fai",
                                       "combined.fa.gzi", "sequences.tsv"}) {
        std::error_code error;
        std::filesystem::remove(output_directory / filename, error);
    }

    PreparationResult result;
    result.sequence_table = output_directory / "sequences.tsv";

    if (combine) {
        result.mode = ComparisonMode::multi_self;
        result.target.fasta = output_directory / "combined.fa";
        set_index_paths(result.target, output_directory, "combined.fa");
        result.query = result.target;
        result.catalog = combine_fastas(
            input_paths, std::vector<std::string>(input_paths.size(), "self"),
            result.target.fasta);
        build_fasta_index(result.target.fasta, result.target.index);
    } else if (input_paths.size() == 1) {
        result.mode = ComparisonMode::self;
        result.target.fasta = input_paths[0];
        set_index_paths(result.target, output_directory, "target.fa");
        result.query = result.target;
        build_fasta_index(result.target.fasta, result.target.index,
                          result.target.gzi);
        add_indexed_sequences(result.catalog, result.target, samples[0],
                              input_paths[0], "self");
    } else if (input_paths.size() == 2) {
        result.mode = ComparisonMode::target_query;
        result.target.fasta = input_paths[0];
        result.query.fasta = input_paths[1];
        set_index_paths(result.target, output_directory, "target.fa");
        set_index_paths(result.query, output_directory, "query.fa");
        build_fasta_index(result.target.fasta, result.target.index,
                          result.target.gzi);
        build_fasta_index(result.query.fasta, result.query.index,
                          result.query.gzi);
        add_indexed_sequences(result.catalog, result.target, samples[0],
                              input_paths[0], "target");
        add_indexed_sequences(result.catalog, result.query, samples[1],
                              input_paths[1], "query");
    } else {
        throw std::runtime_error(
            "more than two FASTA files require the --combine option");
    }

    result.sequence_count = result.catalog.records().size();
    result.total_bases = result.catalog.total_bases();
    result.catalog.write_tsv(result.sequence_table);
    return result;
}

}  // namespace vallescope2
