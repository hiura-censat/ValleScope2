#include "vallescope2/fasta.hpp"

#include <htslib/hts.h>
#include <htslib/kstring.h>

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <unordered_set>

namespace vallescope2 {
namespace {

struct HtsFileDeleter {
    void operator()(htsFile* file) const {
        if (file != nullptr) hts_close(file);
    }
};

struct KString {
    kstring_t value{0, 0, nullptr};
    ~KString() { std::free(value.s); }
};

bool ends_with(const std::string& value, const std::string& suffix) {
    return value.size() >= suffix.size() &&
           value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string sequence_name(const std::string_view header) {
    const auto end = header.find_first_of(" \t");
    const std::string name(header.substr(0, end));
    if (name.empty()) throw std::runtime_error("FASTA header has no sequence name");
    if (name.find('#') != std::string::npos) {
        throw std::runtime_error(
            "sequence name contains reserved character '#': " + name);
    }
    return name;
}

void append_fasta(const std::filesystem::path& input_path,
                  const std::string& sample,
                  const std::string& role,
                  std::ostream& output,
                  SequenceCatalog& catalog) {
    std::unique_ptr<htsFile, HtsFileDeleter> input(hts_open(input_path.c_str(), "r"));
    if (!input) {
        throw std::runtime_error("cannot open FASTA file: " + input_path.string());
    }

    SequenceRecord current;
    std::string output_line;
    output_line.reserve(80);
    KString line;

    const auto finish = [&]() {
        if (current.sequence_id.empty()) return;
        if (current.length == 0) {
            throw std::runtime_error("sequence '" + current.original_id +
                                     "' in " + input_path.string() + " is empty");
        }
        if (!output_line.empty()) {
            output << output_line << '\n';
            output_line.clear();
        }
        catalog.add(current);
    };

    int status = 0;
    while ((status = hts_getline(input.get(), '\n', &line.value)) >= 0) {
        const std::string_view text(line.value.s, line.value.l);
        if (text.empty()) continue;
        if (text.front() == '>') {
            finish();
            current = SequenceRecord{};
            current.sample = sample;
            current.original_id = sequence_name(text.substr(1));
            current.sequence_id = sample + "#1#" + current.original_id;
            current.source_fasta = input_path;
            current.role = role;
            output << '>' << current.sequence_id << '\n';
            continue;
        }
        if (current.sequence_id.empty()) {
            throw std::runtime_error("expected FASTA header in " +
                                     input_path.string());
        }
        for (const unsigned char character : text) {
            if (std::isspace(character)) continue;
            const char base = static_cast<char>(std::toupper(character));
            constexpr std::string_view valid = "ACGTRYSWKMBDHVN";
            if (valid.find(base) == std::string_view::npos) {
                throw std::runtime_error("invalid nucleotide in " +
                                         input_path.string());
            }
            output_line.push_back(base);
            ++current.length;
            if (output_line.size() == 80) {
                output << output_line << '\n';
                output_line.clear();
            }
        }
    }
    if (status < -1) {
        throw std::runtime_error("error while reading FASTA: " +
                                 input_path.string());
    }
    finish();
}

}  // namespace

std::string sample_name_from_path(const std::filesystem::path& path) {
    std::string name = path.filename().string();
    if (ends_with(name, ".gz")) name.resize(name.size() - 3);
    for (const std::string extension : {".fasta", ".fna", ".fas", ".fa"}) {
        if (ends_with(name, extension)) {
            name.resize(name.size() - extension.size());
            break;
        }
    }
    if (name.empty() || name.find('#') != std::string::npos) {
        throw std::runtime_error("invalid sample name: " + name);
    }
    return name;
}

SequenceCatalog combine_fastas(
    const std::vector<std::filesystem::path>& input_paths,
    const std::vector<std::string>& roles,
    const std::filesystem::path& output_path) {
    if (input_paths.empty() || input_paths.size() != roles.size()) {
        throw std::runtime_error("invalid FASTA combination request");
    }
    std::ofstream output(output_path);
    if (!output) {
        throw std::runtime_error("cannot create combined FASTA: " +
                                 output_path.string());
    }
    SequenceCatalog catalog;
    std::unordered_set<std::string> samples;
    for (std::size_t i = 0; i < input_paths.size(); ++i) {
        const std::string sample = sample_name_from_path(input_paths[i]);
        if (!samples.insert(sample).second) {
            throw std::runtime_error("duplicate sample name: " + sample);
        }
        append_fasta(input_paths[i], sample, roles[i], output, catalog);
    }
    if (!output) throw std::runtime_error("failed to write combined FASTA");
    return catalog;
}

}  // namespace vallescope2
