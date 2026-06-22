#include "vallescope2/anchor/anchor_output.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace vallescope2 {
namespace {

std::string json_escape(const std::string& value) {
    std::string escaped;
    for (const char character : value) {
        switch (character) {
            case '\\': escaped += "\\\\"; break;
            case '"': escaped += "\\\""; break;
            case '\n': escaped += "\\n"; break;
            case '\r': escaped += "\\r"; break;
            case '\t': escaped += "\\t"; break;
            default: escaped += character; break;
        }
    }
    return escaped;
}

}  // namespace

void write_anchor_bed(const std::filesystem::path& path,
                      const std::vector<AnchorCandidate>& anchors,
                      const SequenceCatalog& catalog) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot create anchor BED: " + path.string());
    }
    output << std::setprecision(17);
    for (std::size_t i = 0; i < anchors.size(); ++i) {
        const auto& anchor = anchors[i];
        if (anchor.sequence_id >= catalog.records().size()) {
            throw std::runtime_error("invalid anchor sequence ID");
        }
        output << catalog.records()[anchor.sequence_id].sequence_id << '\t'
               << anchor.start << '\t'
               << anchor.start + anchor.length << "\tanchor_" << std::setw(8)
               << std::setfill('0') << i + 1 << std::setfill(' ') << '\t'
               << anchor.score << "\t.\n";
    }
}

void write_anchor_metadata(const std::filesystem::path& path,
                           const AnchorMetadata& metadata,
                           const AnchorSelectionResult& selection) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("cannot create anchor metadata: " +
                                 path.string());
    }
    output << std::setprecision(17)
           << "{\n"
           << "  \"anchor_length\": " << metadata.anchor_length << ",\n"
           << "  \"min_center_distance\": "
           << metadata.min_center_distance << ",\n"
           << "  \"target_density\": " << metadata.target_density << ",\n"
           << "  \"actual_density\": " << selection.actual_density << ",\n"
           << "  \"score_formula\": \"mean(log1p(kmer_frequency))\",\n"
           << "  \"score_cutoff\": ";
    if (selection.has_score_cutoff) {
        output << selection.score_cutoff;
    } else {
        output << "null";
    }
    output << ",\n"
           << "  \"n_anchors\": " << selection.anchors.size() << ",\n"
           << "  \"sequence_length\": " << metadata.sequence_length << ",\n"
           << "  \"genmap_k\": " << metadata.genmap_k << ",\n"
           << "  \"genmap_e\": " << metadata.genmap_e << ",\n"
           << "  \"genmap_threads\": " << metadata.genmap_threads << ",\n"
           << "  \"genmap_command\": \""
           << json_escape(metadata.genmap_command) << "\",\n"
           << "  \"genmap_version\": \""
           << json_escape(metadata.genmap_version) << "\",\n"
           << "  \"ref_fasta\": \""
           << json_escape(metadata.input_fastas.front().string()) << "\",\n"
           << "  \"query_fasta\": \""
           << json_escape((metadata.input_fastas.size() > 1
                               ? metadata.input_fastas[1]
                               : metadata.input_fastas.front()).string())
           << "\",\n"
           << "  \"input_fastas\": [";
    for (std::size_t i = 0; i < metadata.input_fastas.size(); ++i) {
        if (i > 0) {
            output << ", ";
        }
        output << '"' << json_escape(metadata.input_fastas[i].string()) << '"';
    }
    output << "],\n"
           << "  \"combined_fasta_sha256\": \""
           << sha256_file(metadata.genmap_input) << "\"\n"
           << "}\n";
}

}  // namespace vallescope2
