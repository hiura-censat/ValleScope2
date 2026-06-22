#include "vallescope2/context/anchor_grouping.hpp"
#include "vallescope2/context/structural_tokens.hpp"
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
        require(contents.find("a3\tchr1\t6\t12\tGTNNAC\tGTNNAC\t.\t.") !=
                    std::string::npos,
                "N anchor was assigned a group ID");

        {
            std::ofstream output(token_input);
            output << "anchor_id\tsequence_id\tstart\tend\tanchor_seq"
                      "\tcanonical_seq\tanchor_group_id\tnumeric_group_id"
                      "\torientation\tgroupable\n"
                   << "a1\ts1\t0\t6\tAAAAAA\tAAAAAA\tg1\t0\t+\ttrue\n"
                   << "n1\ts1\t20\t26\tAANNAA\tAANNAA\t.\t.\t.\tfalse\n"
                   << "a2\ts1\t49\t55\tCCCCCC\tCCCCCC\tg2\t1\t+\ttrue\n"
                   << "a3\ts1\t99\t105\tGGGGGG\tCCCCCC\tg2\t1\t-\ttrue\n"
                   << "a4\ts2\t0\t6\tTTTTTT\tAAAAAA\tg1\t0\t-\ttrue\n";
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

        std::filesystem::remove_all(directory);
        return 0;
    } catch (const std::exception& error) {
        std::filesystem::remove_all(directory);
        std::cerr << error.what() << '\n';
        return 1;
    }
}
