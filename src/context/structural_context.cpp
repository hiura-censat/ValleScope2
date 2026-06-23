#include "vallescope2/context/structural_context.hpp"

#include "vallescope2/common/sha256.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace vallescope2 {
namespace {

using TokenId = std::uint32_t;

struct AnchorPosition {
    std::string anchor_id;
    std::uint32_t token_position = 0;
    std::vector<std::uint64_t> alpha_prime;
};

struct ContigTokens {
    std::string sequence_id;
    std::uint32_t sample_id = 0;
    std::vector<TokenId> tokens;
    std::vector<std::uint32_t> anchor_prefix{0};
    std::vector<AnchorPosition> anchors;
};

struct Corpus {
    std::vector<std::string> samples;
    std::unordered_map<std::string, std::uint32_t> sample_ids;
    std::unordered_map<std::string, std::uint32_t> sequence_samples;
    std::vector<std::string> token_labels;
    std::unordered_map<std::string, TokenId> token_ids;
    std::vector<bool> token_is_anchor;
    std::vector<ContigTokens> contigs;
};

struct Occurrence {
    std::uint32_t count = 0;
    std::uint32_t contig = 0;
    std::uint32_t start = 0;
    std::uint32_t length = 0;
};

struct Fingerprint {
    std::uint64_t first = 0;
    std::uint64_t second = 0;
    std::uint32_t length = 0;

    bool operator==(const Fingerprint& other) const {
        return first == other.first && second == other.second &&
               length == other.length;
    }
};

struct FingerprintHash {
    std::size_t operator()(const Fingerprint& value) const {
        return static_cast<std::size_t>(
            value.first ^ (value.second + 0x9e3779b97f4a7c15ULL +
                           (value.first << 6) + (value.first >> 2)) ^
            value.length);
    }
};

struct FrequencyBucket {
    Occurrence first;
    std::vector<Occurrence> collisions;
};

struct TMus {
    std::uint32_t sample = 0;
    std::uint32_t contig = 0;
    std::uint32_t start = 0;
    std::uint32_t length = 0;
};

struct PatternTag {
    std::uint32_t sample = 0;
    std::uint32_t length = 0;
};

struct TrieNode {
    std::unordered_map<TokenId, std::uint32_t> next;
    std::uint32_t failure = 0;
    std::vector<PatternTag> output;
};

bool contains_anchor(const ContigTokens& contig,
                     const std::uint32_t start,
                     const std::uint32_t length) {
    return contig.anchor_prefix[start + length] > contig.anchor_prefix[start];
}

int compare_forward_reverse(const std::vector<TokenId>& tokens,
                            const std::uint32_t start,
                            const std::uint32_t length) {
    for (std::uint32_t i = 0; i < length; ++i) {
        const TokenId forward = tokens[start + i];
        const TokenId reverse = tokens[start + length - i - 1];
        if (forward < reverse) return -1;
        if (forward > reverse) return 1;
    }
    return 0;
}

bool canonical_is_reverse(const std::vector<TokenId>& tokens,
                          const std::uint32_t start,
                          const std::uint32_t length) {
    return compare_forward_reverse(tokens, start, length) > 0;
}

Fingerprint fingerprint(const std::vector<TokenId>& tokens,
                        const std::uint32_t start,
                        const std::uint32_t length) {
    const bool reverse = compare_forward_reverse(tokens, start, length) > 0;
    std::uint64_t first = 1469598103934665603ULL;
    std::uint64_t second = 0x9e3779b97f4a7c15ULL;
    for (std::uint32_t i = 0; i < length; ++i) {
        const TokenId value = reverse ? tokens[start + length - i - 1]
                                      : tokens[start + i];
        first ^= static_cast<std::uint64_t>(value) + 1;
        first *= 1099511628211ULL;
        second ^= static_cast<std::uint64_t>(value) + 0x9e3779b9U +
                  (second << 6) + (second >> 2);
    }
    return {first, second, length};
}

std::vector<TokenId> canonical_tokens(const ContigTokens& contig,
                                      const std::uint32_t start,
                                      const std::uint32_t length) {
    std::vector<TokenId> result(contig.tokens.begin() + start,
                                contig.tokens.begin() + start + length);
    if (canonical_is_reverse(contig.tokens, start, length)) {
        std::reverse(result.begin(), result.end());
    }
    return result;
}

std::string join_tokens(const std::vector<TokenId>& tokens,
                        const Corpus& corpus) {
    std::string text;
    for (const TokenId token : tokens) {
        if (!text.empty()) text += ',';
        text += corpus.token_labels.at(token);
    }
    return text;
}

std::string join_token_range(const ContigTokens& contig,
                             const Corpus& corpus,
                             const std::uint32_t begin,
                             const std::uint32_t end) {
    std::string text;
    for (std::uint32_t index = begin; index < end; ++index) {
        if (!text.empty()) text += '|';
        text += corpus.token_labels.at(contig.tokens.at(index));
    }
    return text;
}

std::string stable_serialization(const std::vector<TokenId>& tokens,
                                 const Corpus& corpus) {
    std::string result;
    for (const TokenId token : tokens) {
        const auto& label = corpus.token_labels.at(token);
        const std::uint32_t size = static_cast<std::uint32_t>(label.size());
        result.append(reinterpret_cast<const char*>(&size), sizeof(size));
        result += label;
    }
    return result;
}

void load_sequence_table(const std::filesystem::path& path, Corpus& corpus) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open sequence table: " + path.string());
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::istringstream fields(line);
        std::string sequence_id, sample, haplotype, original_id, source, length, role;
        if (!(fields >> sequence_id >> sample >> haplotype >> original_id >> source >> length >> role)) {
            throw std::runtime_error("invalid sequence table row");
        }
        auto found = corpus.sample_ids.find(sample);
        std::uint32_t sample_id = 0;
        if (found == corpus.sample_ids.end()) {
            sample_id = static_cast<std::uint32_t>(corpus.samples.size());
            corpus.sample_ids.emplace(sample, sample_id);
            corpus.samples.push_back(sample);
        } else {
            sample_id = found->second;
        }
        corpus.sequence_samples.emplace(sequence_id, sample_id);
    }
}

TokenId intern_token(Corpus& corpus, const std::string& label, const bool is_anchor) {
    const auto found = corpus.token_ids.find(label);
    if (found != corpus.token_ids.end()) return found->second;
    const TokenId id = static_cast<TokenId>(corpus.token_labels.size());
    corpus.token_ids.emplace(label, id);
    corpus.token_labels.push_back(label);
    corpus.token_is_anchor.push_back(is_anchor);
    return id;
}

void load_tokens(const std::filesystem::path& path, Corpus& corpus) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open structural tokens: " + path.string());
    std::string line;
    std::getline(input, line);
    std::string current_sequence;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::istringstream fields(line);
        std::string sequence_id, index, type, token, anchor_id, raw_distance, bin;
        if (!(fields >> sequence_id >> index >> type >> token >> anchor_id >> raw_distance >> bin)) {
            throw std::runtime_error("invalid structural token row");
        }
        if (sequence_id != current_sequence) {
            const auto sample = corpus.sequence_samples.find(sequence_id);
            if (sample == corpus.sequence_samples.end()) {
                throw std::runtime_error("structural token sequence has no sample: " + sequence_id);
            }
            ContigTokens contig;
            contig.sequence_id = sequence_id;
            contig.sample_id = sample->second;
            corpus.contigs.push_back(std::move(contig));
            current_sequence = sequence_id;
        }
        auto& contig = corpus.contigs.back();
        const bool is_anchor = type == "A";
        if (!is_anchor && type != "D") throw std::runtime_error("invalid token type");
        const TokenId id = intern_token(corpus, token, is_anchor);
        const auto position = static_cast<std::uint32_t>(contig.tokens.size());
        contig.tokens.push_back(id);
        contig.anchor_prefix.push_back(contig.anchor_prefix.back() + (is_anchor ? 1 : 0));
        if (is_anchor) contig.anchors.push_back({anchor_id, position, {}});
    }
}

using FrequencyMap =
    std::unordered_map<Fingerprint, FrequencyBucket, FingerprintHash>;

TokenId canonical_token_at(const ContigTokens& contig,
                           const std::uint32_t start,
                           const std::uint32_t length,
                           const std::uint32_t offset) {
    return canonical_is_reverse(contig.tokens, start, length)
               ? contig.tokens[start + length - offset - 1]
               : contig.tokens[start + offset];
}

bool same_substring(const Corpus& corpus,
                    const Occurrence& representative,
                    const std::uint32_t contig_id,
                    const std::uint32_t start,
                    const std::uint32_t length) {
    if (representative.length != length) return false;
    const auto& left = corpus.contigs[representative.contig];
    const auto& right = corpus.contigs[contig_id];
    for (std::uint32_t i = 0; i < length; ++i) {
        if (canonical_token_at(left, representative.start, length, i) !=
            canonical_token_at(right, start, length, i)) return false;
    }
    return true;
}

Occurrence* find_occurrence(FrequencyBucket& bucket,
                            const Corpus& corpus,
                            const std::uint32_t contig,
                            const std::uint32_t start,
                            const std::uint32_t length) {
    if (same_substring(corpus, bucket.first, contig, start, length)) {
        return &bucket.first;
    }
    for (auto& occurrence : bucket.collisions) {
        if (same_substring(corpus, occurrence, contig, start, length)) {
            return &occurrence;
        }
    }
    return nullptr;
}

const Occurrence* lookup_occurrence(const FrequencyMap& frequency,
                                    const Corpus& corpus,
                                    const std::uint32_t contig,
                                    const std::uint32_t start,
                                    const std::uint32_t length) {
    const auto found = frequency.find(
        fingerprint(corpus.contigs[contig].tokens, start, length));
    if (found == frequency.end()) return nullptr;
    if (same_substring(corpus, found->second.first, contig, start, length)) {
        return &found->second.first;
    }
    for (const auto& occurrence : found->second.collisions) {
        if (same_substring(corpus, occurrence, contig, start, length)) {
            return &occurrence;
        }
    }
    return nullptr;
}

std::vector<TMus> find_tmus(const Corpus& corpus, const std::uint32_t maximum_length) {
    std::vector<FrequencyMap> frequencies(corpus.samples.size());
    for (std::uint32_t c = 0; c < corpus.contigs.size(); ++c) {
        const auto& contig = corpus.contigs[c];
        auto& frequency = frequencies[contig.sample_id];
        for (std::uint32_t start = 0; start < contig.tokens.size(); ++start) {
            const auto limit = std::min<std::uint32_t>(maximum_length,
                static_cast<std::uint32_t>(contig.tokens.size() - start));
            for (std::uint32_t length = 1; length <= limit; ++length) {
                if (!contains_anchor(contig, start, length)) continue;
                const auto key = fingerprint(contig.tokens, start, length);
                const auto found = frequency.find(key);
                if (found == frequency.end()) {
                    FrequencyBucket bucket;
                    bucket.first = {1, c, start, length};
                    frequency.emplace(key, std::move(bucket));
                } else {
                    auto* occurrence = find_occurrence(
                        found->second, corpus, c, start, length);
                    if (occurrence == nullptr) {
                        found->second.collisions.push_back({1, c, start, length});
                    } else if (occurrence->count == 1) {
                        occurrence->count = 2;
                    }
                }
            }
        }
    }

    std::vector<TMus> result;
    for (std::uint32_t sample = 0; sample < frequencies.size(); ++sample) {
        const auto& frequency = frequencies[sample];
        const auto consider = [&](const Occurrence& occurrence) {
            if (occurrence.count != 1) return;
            const auto& contig = corpus.contigs[occurrence.contig];
            bool minimal = true;
            if (occurrence.length > 1) {
                for (const std::uint32_t trim_left : {0U, 1U}) {
                    const std::uint32_t start = occurrence.start + trim_left;
                    const std::uint32_t length = occurrence.length - 1;
                    if (!contains_anchor(contig, start, length)) continue;
                    const auto* shorter = lookup_occurrence(
                        frequency, corpus, occurrence.contig, start, length);
                    if (shorter != nullptr && shorter->count == 1) {
                        minimal = false;
                        break;
                    }
                }
            }
            if (minimal) result.push_back(
                {sample, occurrence.contig, occurrence.start, occurrence.length});
        };
        for (const auto& item : frequency) {
            consider(item.second.first);
            for (const auto& collision : item.second.collisions) consider(collision);
        }
    }
    std::sort(result.begin(), result.end(), [&](const TMus& left, const TMus& right) {
        if (left.sample != right.sample) return left.sample < right.sample;
        const auto left_tokens = canonical_tokens(
            corpus.contigs[left.contig], left.start, left.length);
        const auto right_tokens = canonical_tokens(
            corpus.contigs[right.contig], right.start, right.length);
        return left_tokens < right_tokens;
    });
    return result;
}

void add_pattern(std::vector<TrieNode>& trie,
                 const std::vector<TokenId>& pattern,
                 const PatternTag tag) {
    std::uint32_t node = 0;
    for (const TokenId token : pattern) {
        auto found = trie[node].next.find(token);
        if (found == trie[node].next.end()) {
            const auto next = static_cast<std::uint32_t>(trie.size());
            trie[node].next.emplace(token, next);
            trie.emplace_back();
            node = next;
        } else node = found->second;
    }
    const auto duplicate = std::find_if(trie[node].output.begin(), trie[node].output.end(),
        [&](const PatternTag value) {
            return value.sample == tag.sample && value.length == tag.length;
        });
    if (duplicate == trie[node].output.end()) trie[node].output.push_back(tag);
}

std::vector<TrieNode> make_automaton(const std::vector<TMus>& tmus,
                                    const Corpus& corpus) {
    std::vector<TrieNode> trie(1);
    for (const auto& item : tmus) {
        auto pattern = canonical_tokens(
            corpus.contigs[item.contig], item.start, item.length);
        add_pattern(trie, pattern, {item.sample, item.length});
        std::reverse(pattern.begin(), pattern.end());
        add_pattern(trie, pattern, {item.sample, item.length});
    }
    std::queue<std::uint32_t> pending;
    for (const auto edge : trie[0].next) pending.push(edge.second);
    while (!pending.empty()) {
        const auto node = pending.front();
        pending.pop();
        for (const auto edge : trie[node].next) {
            const TokenId token = edge.first;
            const std::uint32_t child = edge.second;
            std::uint32_t failure = trie[node].failure;
            while (failure != 0 && trie[failure].next.count(token) == 0) {
                failure = trie[failure].failure;
            }
            const auto found = trie[failure].next.find(token);
            if (found != trie[failure].next.end() && found->second != child) {
                trie[child].failure = found->second;
            }
            const auto& inherited = trie[trie[child].failure].output;
            trie[child].output.insert(trie[child].output.end(), inherited.begin(), inherited.end());
            pending.push(child);
        }
    }
    return trie;
}

void assign_alpha(Corpus& corpus,
                  const std::vector<TrieNode>& trie,
                  const std::uint32_t radius) {
    for (auto& contig : corpus.contigs) {
        const auto sample_count = corpus.samples.size();
        std::vector<std::vector<std::int64_t>> difference(
            sample_count, std::vector<std::int64_t>(contig.anchors.size() + 1));
        std::vector<std::uint32_t> positions;
        positions.reserve(contig.anchors.size());
        for (const auto& anchor : contig.anchors) positions.push_back(anchor.token_position);
        std::uint32_t state = 0;
        for (std::uint32_t end = 0; end < contig.tokens.size(); ++end) {
            const TokenId token = contig.tokens[end];
            while (state != 0 && trie[state].next.count(token) == 0) {
                state = trie[state].failure;
            }
            const auto transition = trie[state].next.find(token);
            state = transition == trie[state].next.end() ? 0 : transition->second;
            for (const auto tag : trie[state].output) {
                const std::uint32_t start = end + 1 - tag.length;
                const std::uint32_t minimum_center = end > radius ? end - radius : 0;
                const std::uint64_t maximum_center =
                    static_cast<std::uint64_t>(start) + radius;
                const auto first = std::lower_bound(positions.begin(), positions.end(), minimum_center);
                const auto last = std::upper_bound(positions.begin(), positions.end(), maximum_center);
                if (first == last) continue;
                const auto begin_index = static_cast<std::size_t>(first - positions.begin());
                const auto end_index = static_cast<std::size_t>(last - positions.begin());
                ++difference[tag.sample][begin_index];
                --difference[tag.sample][end_index];
            }
        }
        for (std::uint32_t sample = 0; sample < sample_count; ++sample) {
            std::int64_t value = 0;
            for (std::size_t i = 0; i < contig.anchors.size(); ++i) {
                value += difference[sample][i];
                if (contig.anchors[i].alpha_prime.empty()) {
                    contig.anchors[i].alpha_prime.resize(sample_count);
                }
                contig.anchors[i].alpha_prime[sample] =
                    static_cast<std::uint64_t>(value);
            }
        }
    }
}

std::vector<TokenId> canonical_context(const ContigTokens& contig,
                                       const Corpus& corpus,
                                       const std::uint32_t center,
                                       const std::uint32_t radius,
                                       std::uint32_t& begin,
                                       std::uint32_t& end) {
    begin = center > radius ? center - radius : 0;
    end = std::min<std::uint32_t>(contig.tokens.size(), center + radius + 1);
    std::vector<TokenId> context(contig.tokens.begin() + begin,
                                 contig.tokens.begin() + end);
    auto reverse = context;
    std::reverse(reverse.begin(), reverse.end());
    const auto labels_less = [&](const std::vector<TokenId>& left,
                                 const std::vector<TokenId>& right) {
        return std::lexicographical_compare(
            left.begin(), left.end(), right.begin(), right.end(),
            [&](const TokenId a, const TokenId b) {
                return corpus.token_labels[a] < corpus.token_labels[b];
            });
    };
    if (labels_less(reverse, context)) return reverse;
    return context;
}

void write_metadata(const std::filesystem::path& path,
                    const Corpus& corpus,
                    const std::uint32_t radius,
                    const StructuralContextResult& result,
                    const std::vector<std::uint64_t>& tmus_per_sample) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create structural context metadata");
    output << "{\n  \"context_radius_tokens\": " << radius
           << ",\n  \"maximum_context_length\": " << 2 * radius + 1
           << ",\n  \"uniqueness_scope\": \"sample\",\n"
           << "  \"strand_invariant\": true,\n"
           << "  \"tmus_requires_anchor_token\": true,\n"
           << "  \"tmus_must_be_fully_contained\": true,\n"
           << "  \"alpha_prime_counting\": \"occurrences\",\n"
           << "  \"n_anchors\": " << result.anchor_count << ",\n"
           << "  \"n_context_groups\": " << result.context_group_count << ",\n"
           << "  \"n_tmus\": " << result.tmus_count << ",\n"
           << "  \"samples\": [";
    for (std::size_t i = 0; i < corpus.samples.size(); ++i) {
        if (i) output << ", ";
        output << "{\"name\": \"" << corpus.samples[i]
               << "\", \"n_tmus\": " << tmus_per_sample[i] << '}';
    }
    output << "]\n}\n";
}

}  // namespace

StructuralContextResult build_structural_contexts(
    const std::filesystem::path& structural_tokens,
    const std::filesystem::path& sequence_table,
    const std::uint32_t context_radius_tokens,
    const std::filesystem::path& anchor_context_output,
    const std::filesystem::path& context_group_output,
    const std::filesystem::path& tmus_output,
    const std::filesystem::path& metadata_output) {
    if (context_radius_tokens >
        (std::numeric_limits<std::uint32_t>::max() - 1) / 2) {
        throw std::runtime_error("context radius is too large");
    }
    Corpus corpus;
    load_sequence_table(sequence_table, corpus);
    load_tokens(structural_tokens, corpus);
    const std::uint32_t maximum_length = 2 * context_radius_tokens + 1;
    const auto tmus = find_tmus(corpus, maximum_length);
    auto automaton = make_automaton(tmus, corpus);
    assign_alpha(corpus, automaton, context_radius_tokens);

    std::vector<std::uint64_t> tmus_per_sample(corpus.samples.size());
    std::ofstream tmus_stream(tmus_output);
    if (!tmus_stream) throw std::runtime_error("cannot create tMUS output");
    tmus_stream << "sample\ttmus_id\ttoken_length\tcanonical_tokens"
                   "\tsource_sequence_id\tsource_start_token\tsource_end_token\n";
    std::vector<std::uint64_t> next_tmus_id(corpus.samples.size());
    for (const auto& item : tmus) {
        const auto tokens = canonical_tokens(
            corpus.contigs[item.contig], item.start, item.length);
        tmus_stream << corpus.samples[item.sample] << '\t'
                    << next_tmus_id[item.sample]++ << '\t' << item.length << '\t'
                    << join_tokens(tokens, corpus) << '\t'
                    << corpus.contigs[item.contig].sequence_id << '\t'
                    << item.start << '\t' << item.start + item.length << '\n';
        ++tmus_per_sample[item.sample];
    }

    struct ContextGroup { std::string exact; std::string tokens; std::uint64_t count = 0; };
    std::unordered_map<std::string, ContextGroup> context_groups;
    std::ofstream anchors(anchor_context_output);
    if (!anchors) throw std::runtime_error("cannot create anchor context output");
    anchors << "anchor_id\tsequence_id\tsample\tcenter_token_index"
               "\tcontext_start_token\tcontext_end_token\tcontext_length"
               "\tcanonical_context_key\tforward_context_tokens";
    for (const auto& sample : corpus.samples) anchors << '\t' << sample << "_alpha_prime";
    anchors << '\n';

    StructuralContextResult result;
    for (const auto& contig : corpus.contigs) {
        for (const auto& anchor : contig.anchors) {
            std::uint32_t begin = 0, end = 0;
            const auto context = canonical_context(
                contig, corpus, anchor.token_position, context_radius_tokens,
                begin, end);
            const auto exact = stable_serialization(context, corpus);
            const auto key = sha256_text(exact);
            auto found = context_groups.find(key);
            if (found == context_groups.end()) {
                context_groups.emplace(key, ContextGroup{exact, join_tokens(context, corpus), 1});
            } else {
                if (found->second.exact != exact) {
                    throw std::runtime_error("SHA-256 collision between structural contexts");
                }
                ++found->second.count;
            }
            anchors << anchor.anchor_id << '\t' << contig.sequence_id << '\t'
                    << corpus.samples[contig.sample_id] << '\t'
                    << anchor.token_position << '\t' << begin << '\t' << end
                    << '\t' << end - begin << '\t' << key << '\t'
                    << join_token_range(contig, corpus, begin, end);
            for (const auto value : anchor.alpha_prime) anchors << '\t' << value;
            anchors << '\n';
            ++result.anchor_count;
        }
    }

    std::vector<std::pair<std::string, ContextGroup>> ordered_groups(
        context_groups.begin(), context_groups.end());
    std::sort(ordered_groups.begin(), ordered_groups.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    std::ofstream groups(context_group_output);
    if (!groups) throw std::runtime_error("cannot create context group output");
    groups << "canonical_context_key\tcanonical_tokens\tanchor_count\n";
    for (const auto& item : ordered_groups) {
        groups << item.first << '\t' << item.second.tokens << '\t'
               << item.second.count << '\n';
    }
    result.context_group_count = ordered_groups.size();
    result.tmus_count = tmus.size();
    write_metadata(metadata_output, corpus, context_radius_tokens, result,
                   tmus_per_sample);
    return result;
}

}  // namespace vallescope2
