#include "vallescope2/anchor/anchor_detection.hpp"

#include "vallescope2/anchor/anchor_output.hpp"
#include "vallescope2/anchor/anchor_selection.hpp"
#include "vallescope2/anchor/window_scoring.hpp"
#include "vallescope2/common/workspace.hpp"
#include "vallescope2/context/anchor_grouping.hpp"
#include "vallescope2/context/structural_tokens.hpp"
#include "vallescope2/context/structural_context.hpp"
#include "vallescope2/correspondence/assignment.hpp"
#include "vallescope2/fasta_index.hpp"
#include "vallescope2/input_preparation.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace vallescope2 {
namespace {

void report_preparation(const PreparationResult& result,
                        const ProgramOptions& options) {
    std::cerr << "Prepared " << result.sequence_count << " sequence(s), "
              << result.total_bases << " bp from " << options.input_paths.size()
              << " FASTA file(s).\n";
    if (result.mode == ComparisonMode::self) std::cerr << "Mode: all-vs-all self comparison.\n";
    else if (result.mode == ComparisonMode::target_query) std::cerr << "Mode: target-query comparison.\n";
    else std::cerr << "Mode: combined multi-FASTA all-vs-all self comparison.\n";
    if (!options.debug) return;
    std::cerr << "Target FASTA: " << result.target.fasta << '\n'
              << "Target index: " << result.target.index << '\n';
    if (result.query.fasta != result.target.fasta) {
        std::cerr << "Query FASTA: " << result.query.fasta << '\n'
                  << "Query index: " << result.query.index << '\n';
    }
    std::cerr << "Sequence table: " << result.sequence_table << '\n';
}

}  // namespace

void run_anchor_detection(const ProgramOptions& options) {
    Workspace workspace(options.debug);
    const auto prepared = prepare_inputs(options.input_paths, workspace.path(), options.combine);
    report_preparation(prepared, options);

    const auto genmap = run_genmap(options.input_paths, prepared, workspace.path(), options.genmap);
    std::cerr << "GenMap profile complete (K=" << options.genmap.kmer_length
              << ", E=" << options.genmap.errors << ").\n";
    if (options.debug) {
        std::cerr << "GenMap raw frequencies: " << genmap.frequency_raw << '\n'
                  << "GenMap log: " << genmap.log << '\n';
    }

    std::ofstream window_output;
    const auto window_score_path = std::filesystem::current_path() / "window_scores.tsv";
    if (options.dump_window_scores) {
        window_output.open(window_score_path);
        if (!window_output) throw std::runtime_error("cannot create window score output");
        window_output << "sequence\tstart\tend\tscore\n" << std::setprecision(17);
    }

    CandidateSpool spool(workspace.path(), options.min_center_distance);
    const auto scoring_start = std::chrono::steady_clock::now();
    const auto scoring = score_windows(
        genmap.frequency_raw, genmap.catalog,
        {options.anchor_length, options.genmap.kmer_length},
        [&](const WindowScore& window) {
            spool.add({window.start, window.score, 0, window.sequence_id,
                       static_cast<std::uint32_t>(window.end - window.start)});
            if (window_output) {
                window_output << genmap.catalog.records()[window.sequence_id].sequence_id
                              << '\t' << window.start << '\t' << window.end
                              << '\t' << window.score << '\n';
            }
        });
    const std::chrono::duration<double> scoring_time =
        std::chrono::steady_clock::now() - scoring_start;
    std::cerr << "Scored " << scoring.window_count << " window(s) at "
              << options.anchor_length << " bp in " << scoring_time.count() << " s.\n";
    if (scoring.gap_positions > 0)
        std::cerr << "Skipped " << scoring.gap_positions << " uncovered k-mer position(s).\n";
    if (scoring.saturated_positions > 0)
        std::cerr << "Warning: " << scoring.saturated_positions
                  << " GenMap frequency value(s) reached 65535.\n";
    if (options.dump_window_scores) std::cerr << "Window scores: " << window_score_path << '\n';

    const auto selection_start = std::chrono::steady_clock::now();
    const auto selection = select_anchors_external(
        spool, genmap.catalog.records().size(), genmap.catalog.total_bases(),
        {options.min_center_distance, options.target_density});
    const std::chrono::duration<double> selection_time =
        std::chrono::steady_clock::now() - selection_start;

    const auto current = std::filesystem::current_path();
    const auto anchor_bed = current / "anchors.bed";
    const auto anchor_metadata = current / "anchors.meta.json";
    const auto genmap_log = current / "genmap.log";
    const auto grouped_anchors = current / "grouped_anchors.tsv";
    const auto anchor_groups = current / "anchor_groups.tsv";
    const auto structural_tokens = current / "structural_tokens.tsv";
    const auto structural_token_metadata =
        current / "structural_tokens.meta.json";
    const auto anchor_contexts = current / "anchor_contexts.tsv";
    const auto context_groups = current / "context_groups.tsv";
    const auto tmus = current / "tmus.tsv";
    const auto structural_context_metadata =
        current / "structural_contexts.meta.json";
    const auto assignments = current / "assignments.tsv";
    const auto assignment_metadata = current / "assignments.meta.json";
    write_anchor_bed(anchor_bed, selection.anchors, genmap.catalog);

    const auto context_index = workspace.path() / "context_input.fa.fai";
    build_fasta_index(genmap.input_fasta, context_index);
    const auto grouping = group_anchors(
        anchor_bed, genmap.input_fasta, context_index, options.anchor_length,
        grouped_anchors, anchor_groups);
    const auto tokenization = build_structural_tokens(
        grouped_anchors, options.distance_bin_size, structural_tokens,
        structural_token_metadata);
    const auto context_start = std::chrono::steady_clock::now();
    const auto contexts = build_structural_contexts(
        structural_tokens, prepared.sequence_table,
        options.context_radius_tokens, anchor_contexts, context_groups, tmus,
        structural_context_metadata);
    const std::chrono::duration<double> context_time =
        std::chrono::steady_clock::now() - context_start;
    const auto assignment_start = std::chrono::steady_clock::now();
    const auto assignment = assign_context_correspondences(
        grouped_anchors, anchor_contexts, prepared.sequence_table, assignments,
        assignment_metadata, options.assignment);
    const std::chrono::duration<double> assignment_time =
        std::chrono::steady_clock::now() - assignment_start;
    write_anchor_metadata(
        anchor_metadata,
        {options.anchor_length, options.min_center_distance, options.target_density,
         genmap.catalog.total_bases(), options.genmap.kmer_length,
         options.genmap.errors, options.genmap.threads,
         genmap.command, genmap.version,
         options.input_paths, genmap.input_fasta},
        selection);
    std::filesystem::copy_file(genmap.log, genmap_log,
                               std::filesystem::copy_options::overwrite_existing);

    for (const auto& path : {spool.candidate_path(), spool.group_path()}) {
        std::error_code error;
        std::filesystem::remove(path, error);
    }
    std::cerr << "Selected " << selection.anchors.size() << " anchor(s); actual density "
              << selection.actual_density << " anchors/Mb in "
              << selection_time.count() << " s.\n"
              << "Anchor BED: " << anchor_bed << '\n'
              << "Anchor metadata: " << anchor_metadata << '\n'
              << "Grouped anchors: " << grouped_anchors << '\n'
              << "Anchor groups: " << anchor_groups << '\n'
              << "Constructed " << grouping.group_count << " anchor group(s); "
              << grouping.ungrouped_anchor_count
              << " anchor(s) containing N were not grouped.\n"
              << "Structural tokens: " << structural_tokens << '\n'
              << "Structural token metadata: " << structural_token_metadata
              << '\n'
              << "Constructed " << tokenization.token_count()
              << " structural token(s) using " << options.distance_bin_size
              << " bp distance bins.\n"
              << "Anchor contexts: " << anchor_contexts << '\n'
              << "Context groups: " << context_groups << '\n'
              << "tMUS table: " << tmus << '\n'
              << "Structural context metadata: "
              << structural_context_metadata << '\n'
              << "Constructed " << contexts.anchor_count << " context(s), "
              << contexts.context_group_count << " canonical context group(s), and "
              << contexts.tmus_count << " sample-specific tMUS in "
              << context_time.count() << " s.\n"
              << "Assignments: " << assignments << '\n'
              << "Assignment metadata: " << assignment_metadata << '\n'
              << "Assigned " << assignment.primary_count << " primary, "
              << assignment.ambiguous_query_count << " ambiguous query, and "
              << assignment.unmatched_count << " unmatched query anchor(s) across "
              << assignment.ordered_pair_count << " ordered pair(s) in "
              << assignment_time.count() << " s.\n"
              << "GenMap log: " << genmap_log << '\n';
}

}  // namespace vallescope2
