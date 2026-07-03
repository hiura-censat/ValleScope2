#include "vallescope2/alignment/base_alignment.hpp"
#include "vallescope2/anchor/anchor_detection.hpp"

#include "vallescope2/anchor/anchor_output.hpp"
#include "vallescope2/anchor/anchor_selection.hpp"
#include "vallescope2/anchor/window_scoring.hpp"
#include "vallescope2/bundle/chaining.hpp"
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

void stream_file_to_stdout(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open final PAF output");
    std::cout << input.rdbuf();
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

    const auto output_directory = workspace.path();
    const auto anchor_bed = output_directory / "anchors.bed";
    const auto anchor_metadata = output_directory / "anchors.meta.json";
    const auto genmap_log = output_directory / "genmap.log";
    const auto grouped_anchors = output_directory / "grouped_anchors.tsv";
    const auto anchor_groups = output_directory / "anchor_groups.tsv";
    const auto structural_tokens = output_directory / "structural_tokens.tsv";
    const auto structural_token_metadata =
        output_directory / "structural_tokens.meta.json";
    const auto anchor_contexts = output_directory / "anchor_contexts.tsv";
    const auto context_groups = output_directory / "context_groups.tsv";
    const auto tmus = output_directory / "tmus.tsv";
    const auto structural_context_metadata =
        output_directory / "structural_contexts.meta.json";
    const auto assignments = output_directory / "assignments.tsv";
    const auto assignment_metadata = output_directory / "assignments.meta.json";
    const auto correspondences = output_directory / "anchor_correspondences.tsv";
    const auto correspondence_metadata =
        output_directory / "anchor_correspondences.meta.json";
    const auto chains = output_directory / "chains.tsv";
    const auto chain_anchors = output_directory / "chain_anchors.tsv";
    const auto chain_metadata = output_directory / "chains.meta.json";
    const auto refined_chains = output_directory / "refined_chains.tsv";
    const auto refined_chain_anchors = output_directory / "refined_chain_anchors.tsv";
    const auto refined_chain_metadata = output_directory / "refined_chains.meta.json";
    const auto bundle_alignments = output_directory / "bundle_alignments.paf";
    const auto bundle_alignment_metadata =
        output_directory / "bundle_alignments.meta.json";
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
        correspondences, correspondence_metadata, assignment_metadata,
        options.assignment);
    const std::chrono::duration<double> assignment_time =
        std::chrono::steady_clock::now() - assignment_start;
    const auto chaining_start = std::chrono::steady_clock::now();
    const auto chaining = build_anchor_chains(
        correspondences, assignments, grouped_anchors, chains, chain_anchors,
        chain_metadata, refined_chains, refined_chain_anchors,
        refined_chain_metadata, options.chaining);
    const std::chrono::duration<double> chaining_time =
        std::chrono::steady_clock::now() - chaining_start;
    BaseAlignmentResult base_alignment;
    std::chrono::duration<double> base_alignment_time{0};
    if (options.base_align) {
        const auto base_alignment_start = std::chrono::steady_clock::now();
        BaseAlignmentParameters base_parameters;
        base_parameters.anchor_length = options.anchor_length;
        base_parameters.max_bundle_align_bp = options.max_bundle_align_bp;
        base_parameters.max_patch_gap_bp = options.max_patch_gap_bp;
        base_parameters.patch_window_bp = options.patch_window_bp;
        base_parameters.min_patch_identity = options.min_patch_identity;
        base_parameters.max_wfa_memory_gb = options.max_wfa_memory_gb;
        base_parameters.bundle_trim_overlap = options.chaining.chain_trim_overlap;
        base_parameters.chain_extension = options.chain_extension;
        base_parameters.chain_extension_predecessors =
            options.chaining.predecessor_count;
        base_parameters.max_chain_extension_bp =
            options.max_chain_extension_bp;
        base_parameters.max_chain_gap = options.chaining.max_chain_gap;
        base_parameters.chain_max_gap_ratio =
            options.chaining.chain_max_gap_ratio;
        base_parameters.gap_cost_model = options.chaining.gap_cost_model;
        base_parameters.gap_weight = options.chaining.gap_weight;
        base_parameters.gap_unit = options.chaining.gap_unit;
        base_parameters.min_chain_extension_anchors =
            options.min_chain_extension_anchors;
        base_parameters.min_chain_extension_score =
            options.min_chain_extension_score;
        base_parameters.min_copy_support_anchors =
            options.min_copy_support_anchors;
        base_alignment = align_chain_bundles(
            {chains}, {chain_anchors},
            grouped_anchors, correspondences, assignments, genmap.input_fasta,
            context_index, bundle_alignments, bundle_alignment_metadata,
            base_parameters);
        base_alignment_time =
            std::chrono::steady_clock::now() - base_alignment_start;
    }
    write_anchor_metadata(
        anchor_metadata,
        {options.anchor_length, options.min_center_distance, options.target_density,
         genmap.catalog.total_bases(), options.genmap.kmer_length,
         options.genmap.errors, options.genmap.threads,
         genmap.command, genmap.version,
         options.input_paths, genmap.input_fasta},
        selection);
    if (std::filesystem::absolute(genmap.log) !=
        std::filesystem::absolute(genmap_log)) {
        std::filesystem::copy_file(
            genmap.log, genmap_log,
            std::filesystem::copy_options::overwrite_existing);
    }

    for (const auto& path : {spool.candidate_path(), spool.group_path()}) {
        std::error_code error;
        std::filesystem::remove(path, error);
    }
    std::cerr << "Selected " << selection.anchors.size()
              << " anchor(s); actual density " << selection.actual_density
              << " anchors/Mb in " << selection_time.count() << " s.\n";
    if (options.debug) {
        std::cerr << "Anchor BED: " << anchor_bed << '\n'
                  << "Anchor metadata: " << anchor_metadata << '\n'
                  << "Grouped anchors: " << grouped_anchors << '\n'
                  << "Anchor groups: " << anchor_groups << '\n';
    }
    std::cerr << "Constructed " << grouping.group_count << " anchor group(s); "
              << grouping.ungrouped_anchor_count
              << " anchor(s) containing N were not grouped.\n";
    if (options.debug) {
        std::cerr << "Structural tokens: " << structural_tokens << '\n'
                  << "Structural token metadata: " << structural_token_metadata
                  << '\n';
    }
    std::cerr << "Constructed " << tokenization.token_count()
              << " structural token(s) using " << options.distance_bin_size
              << " bp distance bins.\n";
    if (options.debug) {
        std::cerr << "Anchor contexts: " << anchor_contexts << '\n'
                  << "Context groups: " << context_groups << '\n'
                  << "tMUS table: " << tmus << '\n'
                  << "Structural context metadata: "
                  << structural_context_metadata << '\n';
    }
    std::cerr << "Constructed " << contexts.anchor_count << " context(s), "
              << contexts.context_group_count << " canonical context group(s), and "
              << contexts.tmus_count << " sample-specific tMUS in "
              << context_time.count() << " s.\n";
    if (options.debug) {
        std::cerr << "Assignments: " << assignments << '\n'
                  << "Assignment metadata: " << assignment_metadata << '\n'
                  << "Anchor correspondences: " << correspondences << '\n'
                  << "Anchor correspondence metadata: " << correspondence_metadata
                  << '\n';
    }
    std::cerr << "Assigned " << assignment.primary_count << " primary, "
              << assignment.ambiguous_query_count << " ambiguous query, and "
              << assignment.unmatched_count << " unmatched query anchor(s) across "
              << assignment.ordered_pair_count << " ordered pair(s) in "
              << assignment_time.count() << " s; merged "
              << assignment.correspondence_count
              << " anchor correspondence(s) using "
              << (options.assignment.pair_merge_mode ==
                          PairMergeMode::reciprocal
                      ? "reciprocal"
                      : "union")
              << " mode.\n";
    if (options.debug) {
        std::cerr << "Chains: " << chains << '\n'
                  << "Chain anchors: " << chain_anchors << '\n'
                  << "Chain metadata: " << chain_metadata << '\n'
                  << "Refined chains: " << refined_chains << '\n'
                  << "Refined chain anchors: " << refined_chain_anchors << '\n'
                  << "Refined chain metadata: " << refined_chain_metadata << '\n';
    }
    std::cerr << "Built " << chaining.chain_count << " top chain(s) from "
              << chaining.candidate_count << " correspondence candidate(s), using "
              << chaining.chain_anchor_count << " chained anchor pair(s); trimmed "
              << chaining.trimmed_chain_count << " of "
              << chaining.raw_chain_count << " raw chain(s); refined "
              << chaining.refined_chain_count << " chain(s), using "
              << chaining.refined_chain_anchor_count << " anchor pair(s) in "
              << chaining_time.count() << " s.\n";
    if (options.base_align) {
        if (options.debug) {
            std::cerr << "Bundle alignments: " << bundle_alignments << '\n'
                      << "Bundle alignment metadata: "
                      << bundle_alignment_metadata << '\n';
        }
        std::cerr << "Base-aligned " << base_alignment.aligned_bundle_count
                  << " of " << base_alignment.bundle_count
                  << " bundle(s); skipped "
                  << base_alignment.skipped_bundle_count
                  << "; chain-extended " << base_alignment.extension_partial_count
                  << " partial path(s) from " << base_alignment.extension_count
                  << " adjacent gap(s); patched " << base_alignment.patch_count
                  << " gap(s) into " << base_alignment.patched_bundle_count
                  << " merged bundle(s); final-trimmed "
                  << base_alignment.bundle_trim_count
                  << " bundle(s) before alignment; "
                  << base_alignment.anchor_guided_alignment_count
                  << " anchor-guided alignment(s) in "
                  << base_alignment_time.count() << " s.\n";
    }
    if (options.debug) std::cerr << "GenMap log: " << genmap_log << '\n';
    if (options.base_align) stream_file_to_stdout(bundle_alignments);
}

}  // namespace vallescope2
