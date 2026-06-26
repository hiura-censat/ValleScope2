#include "vallescope2/alignment/base_alignment.hpp"
#include "vallescope2/context/anchor_grouping.hpp"
#include "vallescope2/context/structural_tokens.hpp"
#include "vallescope2/context/structural_context.hpp"
#include "vallescope2/correspondence/assignment.hpp"
#include "vallescope2/bundle/chaining.hpp"
#include "vallescope2/fasta_index.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(const bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

}  // namespace

int main() {
    const auto directory = std::filesystem::temp_directory_path() /
        ("vallescope2-grouping-test-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        std::filesystem::create_directories(directory);
        const auto fasta = directory / "input.fa";
        const auto index = directory / "input.fa.fai";
        const auto bed = directory / "anchors.bed";
        const auto grouped = directory / "grouped.tsv";
        const auto groups = directory / "groups.tsv";
        const auto token_input = directory / "token_input.tsv";
        const auto tokens = directory / "tokens.tsv";
        const auto token_metadata = directory / "tokens.meta.json";
        const auto context_tokens = directory / "context_tokens.tsv";
        const auto sequence_table = directory / "sequences.tsv";
        const auto anchor_contexts = directory / "anchor_contexts.tsv";
        const auto context_groups = directory / "context_groups.tsv";
        const auto tmus = directory / "tmus.tsv";
        const auto context_metadata = directory / "contexts.meta.json";
        const auto assignment_grouped = directory / "assignment_grouped.tsv";
        const auto assignment_contexts = directory / "assignment_contexts.tsv";
        const auto assignments = directory / "assignments.tsv";
        const auto assignment_metadata = directory / "assignments.meta.json";
        const auto correspondences = directory / "anchor_correspondences.tsv";
        const auto correspondence_metadata =
            directory / "anchor_correspondences.meta.json";
        const auto chains = directory / "chains.tsv";
        const auto chain_anchors = directory / "chain_anchors.tsv";
        const auto chain_metadata = directory / "chains.meta.json";
        const auto refined_chains = directory / "refined_chains.tsv";
        const auto refined_chain_anchors = directory / "refined_chain_anchors.tsv";
        const auto refined_chain_metadata = directory / "refined_chains.meta.json";
        const auto bundle_paf = directory / "bundle_alignments.paf";
        const auto bundle_alignment_metadata =
            directory / "bundle_alignments.meta.json";
        const auto base_chains = directory / "base_chains.tsv";
        {
            std::ofstream output(fasta);
            output << ">chr1\nACGTACGTNNACGT\n";
        }
        {
            std::ofstream output(bed);
            output << "chr1\t0\t6\ta1\n"
                   << "chr1\t2\t8\ta2\n"
                   << "chr1\t6\t12\ta3\n";
        }

        vallescope2::build_fasta_index(fasta, index);
        const auto result = vallescope2::group_anchors(
            bed, fasta, index, 6, grouped, groups);
        require(result.anchor_count == 3, "unexpected anchor count");
        require(result.grouped_anchor_count == 2, "unexpected grouped count");
        require(result.ungrouped_anchor_count == 1, "N anchor was grouped");
        require(result.group_count == 1, "reverse complements formed different groups");
        require(vallescope2::reverse_complement("ARYN") == "NRYT",
                "IUPAC reverse complement failed");

        std::ifstream input(grouped);
        const std::string contents((std::istreambuf_iterator<char>(input)),
                                   std::istreambuf_iterator<char>());
        require(contents.find("a1\tchr1\t0\t6\tACGTAC\tACGTAC") != std::string::npos,
                "forward canonical sequence missing");
        require(contents.find("a2\tchr1\t2\t8\tGTACGT\tACGTAC") != std::string::npos,
                "reverse canonical sequence missing");
        require(contents.find("a3\tchr1\t6\t12\tGTNNAC\tGTNNAC\t.\t.\tfalse") !=
                    std::string::npos,
                "N anchor was assigned a group ID");

        {
            std::ofstream output(token_input);
            output << "anchor_id\tsequence_id\tstart\tend\tanchor_seq"
                      "\tcanonical_seq\tanchor_group_id\torientation\tgroupable\n"
                   << "a1\ts1\t0\t6\tAAAAAA\tAAAAAA\tg1\t+\ttrue\n"
                   << "n1\ts1\t20\t26\tAANNAA\tAANNAA\t.\t.\tfalse\n"
                   << "a2\ts1\t49\t55\tCCCCCC\tCCCCCC\tg2\t+\ttrue\n"
                   << "a3\ts1\t99\t105\tGGGGGG\tCCCCCC\tg2\t-\ttrue\n"
                   << "a4\ts2\t0\t6\tTTTTTT\tAAAAAA\tg1\t-\ttrue\n";
        }
        const auto token_result = vallescope2::build_structural_tokens(
            token_input, 50, tokens, token_metadata);
        require(token_result.anchor_token_count == 4,
                "unexpected anchor token count");
        require(token_result.distance_token_count == 2,
                "unexpected distance token count");
        require(token_result.sequence_count == 2,
                "distance crossed a sequence boundary");
        require(token_result.skipped_anchor_count == 1,
                "ungroupable anchor was not skipped");
        std::ifstream token_stream(tokens);
        const std::string token_contents(
            (std::istreambuf_iterator<char>(token_stream)),
            std::istreambuf_iterator<char>());
        require(token_contents.find("s1\t1\tD\tD:0\t.\t49\t0") !=
                    std::string::npos,
                "49 bp distance bin failed");
        require(token_contents.find("s1\t3\tD\tD:1\t.\t50\t1") !=
                    std::string::npos,
                "50 bp distance bin failed");
        require(token_contents.find("s2\t0\tA\tA:g1\ta4\t.\t.") !=
                    std::string::npos,
                "sequence token index was not reset");

        {
            std::ofstream output(sequence_table);
            output << "sequence_id\tsample\thaplotype\toriginal_id\tsource_file"
                      "\tlength\trole\n"
                   << "s1a\tsample1\t1\tc1\ta.fa\t10\tself\n"
                   << "s1b\tsample1\t1\tc2\ta.fa\t10\tself\n"
                   << "s2seq\tsample2\t1\tc3\tb.fa\t10\tself\n";
        }
        {
            std::ofstream output(context_tokens);
            output << "sequence_id\ttoken_index\ttoken_type\ttoken\tanchor_id"
                      "\traw_distance\tdistance_bin\n"
                   << "s1a\t0\tA\tA:x\ta1\t.\t.\n"
                   << "s1a\t1\tD\tD:1\t.\t50\t1\n"
                   << "s1a\t2\tA\tA:y\ta2\t.\t.\n"
                   << "s1b\t0\tA\tA:y\ta3\t.\t.\n"
                   << "s1b\t1\tD\tD:1\t.\t50\t1\n"
                   << "s1b\t2\tA\tA:x\ta4\t.\t.\n"
                   << "s2seq\t0\tA\tA:z\tb1\t.\t.\n"
                   << "s2seq\t1\tD\tD:2\t.\t100\t2\n"
                   << "s2seq\t2\tA\tA:q\tb2\t.\t.\n";
        }
        const auto context_result = vallescope2::build_structural_contexts(
            context_tokens, sequence_table, 2, anchor_contexts, context_groups,
            tmus, context_metadata);
        require(context_result.anchor_count == 6, "unexpected context count");
        require(context_result.context_group_count == 2,
                "reverse contexts did not canonicalize together");
        require(context_result.tmus_count == 2,
                "reverse-equivalent substrings were treated as unique");
        std::ifstream anchor_context_stream(anchor_contexts);
        const std::string anchor_context_contents(
            (std::istreambuf_iterator<char>(anchor_context_stream)),
            std::istreambuf_iterator<char>());
        require(anchor_context_contents.find(
                    "canonical_context_key\tforward_context_tokens\t"
                    "tmus_keys\tsample1_alpha_prime") != std::string::npos,
                "forward_context_tokens header is missing");
        const auto b1 = anchor_context_contents.find("b1\ts2seq\tsample2\t0\t0\t3\t3\t");
        require(b1 != std::string::npos, "sample2 context missing");
        const auto b1_end = anchor_context_contents.find('\n', b1);
        const auto b1_line = anchor_context_contents.substr(b1, b1_end - b1);
        require(b1_line.find("\tA:z|D:2|A:q\t") != std::string::npos,
                "forward context tokens are incorrect");
        require(b1_line.rfind("\t0\t2") != std::string::npos,
                "sample-specific alpha prime is incorrect");

        {
            std::ofstream output(assignment_grouped);
            output << "anchor_id\tsequence_id\tstart\tend\tanchor_seq"
                      "\tcanonical_seq\tanchor_group_id\torientation\tgroupable\n"
                   << "t1\tseqA\t0\t50\tA\tA\tg1\t+\ttrue\n"
                   << "t2\tseqA\t100\t150\tA\tA\tg2\t+\ttrue\n"
                   << "t3\tseqA\t200\t250\tA\tA\tg3\t+\ttrue\n"
                   << "t4\tseqA\t300\t350\tA\tA\tg3\t+\ttrue\n"
                   << "q1\tseqB\t0\t50\tA\tA\tg1\t+\ttrue\n"
                   << "q2\tseqB\t100\t150\tA\tA\tg2\t+\ttrue\n"
                   << "q3\tseqB\t200\t250\tA\tA\tg3\t+\ttrue\n"
                   << "q4\tseqB\t300\t350\tA\tA\tg4\t+\ttrue\n";
        }
        {
            std::ofstream output(assignment_contexts);
            output << "anchor_id\tsequence_id\tsample\tcenter_token_index"
                      "\tcontext_start_token\tcontext_end_token\tcontext_length"
                      "\tcanonical_context_key\tforward_context_tokens"
                      "\tsampleA_alpha_prime\tsampleB_alpha_prime\n"
                   << "t1\tseqA\tsampleA\t0\t0\t3\t3\tk1\tA:x|D:1|A:y\t0\t3\n"
                   << "t2\tseqA\tsampleA\t0\t0\t3\t3\tkt2\tA:b|D:1|A:a\t0\t5\n"
                   << "t3\tseqA\tsampleA\t0\t0\t1\t1\tkt3\tA:m\t0\t6\n"
                   << "t4\tseqA\tsampleA\t0\t0\t1\t1\tkt4\tA:m\t0\t6\n"
                   << "q1\tseqB\tsampleB\t0\t0\t3\t3\tk1\tA:x|D:1|A:y\t2\t0\n"
                   << "q2\tseqB\tsampleB\t0\t0\t3\t3\tkq2\tA:a|D:1|A:b\t5\t0\n"
                   << "q3\tseqB\tsampleB\t0\t0\t1\t1\tkq3\tA:m\t6\t0\n"
                   << "q4\tseqB\tsampleB\t0\t0\t1\t1\tkq4\tA:z\t7\t0\n";
        }
        {
            std::ofstream output(sequence_table);
            output << "sequence_id\tsample\thaplotype\toriginal_id\tsource_file"
                      "\tlength\trole\n"
                   << "seqA\tsampleA\t1\tseqA\ta.fa\t1000\tself\n"
                   << "seqB\tsampleB\t1\tseqB\tb.fa\t1000\tself\n";
        }
        vallescope2::AssignmentParameters assignment_parameters;
        assignment_parameters.beta_tolerance = 2;
        assignment_parameters.primary_margin = 5;
        const auto assignment_result =
            vallescope2::assign_context_correspondences(
                assignment_grouped, assignment_contexts, sequence_table,
                assignments, correspondences, correspondence_metadata,
                assignment_metadata, assignment_parameters);
        require(assignment_result.ordered_pair_count == 2,
                "unexpected ordered pair count");
        require(assignment_result.exact_context_primary_count >= 1,
                "exact context primary assignment missing");
        require(assignment_result.edit_distance_primary_count >= 1,
                "edit-distance primary assignment missing");
        require(assignment_result.ambiguous_query_count >= 1,
                "ambiguous assignment missing");
        require(assignment_result.unmatched_count >= 1,
                "unmatched assignment missing");
        std::ifstream assignment_stream(assignments);
        const std::string assignment_contents(
            (std::istreambuf_iterator<char>(assignment_stream)),
            std::istreambuf_iterator<char>());
        require(assignment_contents.find(
                    "sampleA\tsampleB\tq1\tt1\tprimary\texact_context\t."
                ) != std::string::npos,
                "exact primary row is incorrect");
        require(assignment_contents.find(
                    "sampleA\tsampleB\tq2\tt2\tprimary\tedit_distance\t-"
                ) != std::string::npos,
                "reverse-strand edit-distance row is incorrect");
        require(assignment_contents.find(
                    "sampleA\tsampleB\tq3\tt3\tambiguous\tedit_distance"
                ) != std::string::npos &&
                assignment_contents.find(
                    "sampleA\tsampleB\tq3\tt4\tambiguous\tedit_distance"
                ) != std::string::npos,
                "ambiguous rows are incorrect");
        require(assignment_contents.find(
                    "sampleA\tsampleB\tq4\t.\tunmatched\tunmatched"
                ) != std::string::npos,
                "unmatched row is incorrect");
        std::ifstream correspondence_stream(correspondences);
        const std::string correspondence_contents(
            (std::istreambuf_iterator<char>(correspondence_stream)),
            std::istreambuf_iterator<char>());
        require(correspondence_contents.find(
                    "sampleA\tsampleB\tt1\tq1\tboth"
                ) != std::string::npos,
                "reciprocal correspondence is missing");
        const auto chaining_result = vallescope2::build_anchor_chains(
            correspondences, assignments, assignment_grouped, chains, chain_anchors,
            chain_metadata, refined_chains, refined_chain_anchors,
            refined_chain_metadata, {2, 50000, 1.0, 50, 1, 0.0, 50000, 1});
        require(chaining_result.candidate_count >= 1,
                "chain candidates were not loaded");
        require(chaining_result.chain_count >= 1,
                "chain was not constructed");
        std::ifstream chain_anchor_stream(chain_anchors);
        const std::string chain_anchor_contents(
            (std::istreambuf_iterator<char>(chain_anchor_stream)),
            std::istreambuf_iterator<char>());
        require(chain_anchor_contents.find("t1\tq1") != std::string::npos,
                "expected chained anchor pair is missing");
        {
            std::ofstream output(base_chains);
            output << "chain_id\tsample_a\tsample_b\tsequence_a\tsequence_b"
                      "\tassign_strand\tn_candidates\tchain_score"
                      "\tref_start\tref_end\tquery_start\tquery_end\n"
                   << "0\tsampleA\tsampleB\tchr1\tchr1\t+\t1\t1\t0\t6\t0\t6\n";
        }
        const auto base_alignment_result = vallescope2::align_chain_bundles(
            base_chains, fasta, index, bundle_paf, bundle_alignment_metadata,
            {6, 1000, 1000000});
        require(base_alignment_result.aligned_bundle_count >= 1,
                "bundle base alignment was not produced");
        std::ifstream bundle_paf_stream(bundle_paf);
        const std::string bundle_paf_contents(
            (std::istreambuf_iterator<char>(bundle_paf_stream)),
            std::istreambuf_iterator<char>());
        require(bundle_paf_contents.find("cg:Z:") != std::string::npos,
                "bundle PAF is missing cg:Z CIGAR");

        std::filesystem::remove_all(directory);
        return 0;
    } catch (const std::exception& error) {
        std::filesystem::remove_all(directory);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
