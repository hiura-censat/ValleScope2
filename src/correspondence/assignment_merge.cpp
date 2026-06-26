#include "assignment_internal.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

namespace vallescope2 {
namespace {

std::string merge_key(const std::uint32_t sample_a,
                      const std::uint32_t sample_b,
                      const std::string& anchor_a,
                      const std::string& anchor_b) {
    constexpr char separator = '\x1f';
    return std::to_string(sample_a) + separator + std::to_string(sample_b) +
           separator + anchor_a + separator + anchor_b;
}

std::string anchor_key(const std::uint32_t sample_a,
                       const std::uint32_t sample_b,
                       const std::string& anchor) {
    constexpr char separator = '\x1f';
    return std::to_string(sample_a) + separator + std::to_string(sample_b) +
           separator + anchor;
}

CanonicalPrimary canonical_primary(const DirectedPrimary& primary) {
    CanonicalPrimary item;
    item.record = &primary;
    if (primary.target_sample < primary.query_sample) {
        item.anchor_a = primary.target_anchor_id;
        item.anchor_b = primary.query_anchor_id;
        item.forward = true;
    } else {
        item.anchor_a = primary.query_anchor_id;
        item.anchor_b = primary.target_anchor_id;
        item.forward = false;
    }
    return item;
}

void write_correspondence_metadata(const std::filesystem::path& path,
                                   const AssignmentParameters& parameters,
                                   const AssignmentResult& result) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot create correspondence metadata");
    output << "{\n"
           << "  \"pair_merge_mode\": \""
           << pair_merge_mode_name(parameters.pair_merge_mode) << "\",\n"
           << "  \"correspondence_count\": " << result.correspondence_count
           << ",\n"
           << "  \"reciprocal_correspondence_count\": "
           << result.reciprocal_correspondence_count << ",\n"
           << "  \"forward_only_correspondence_count\": "
           << result.forward_only_correspondence_count << ",\n"
           << "  \"reverse_only_correspondence_count\": "
           << result.reverse_only_correspondence_count << ",\n"
           << "  \"discordant_correspondence_count\": "
           << result.discordant_correspondence_count << "\n"
           << "}\n";
}

}  // namespace

void write_correspondences(const std::filesystem::path& output_path,
                           const std::filesystem::path& metadata_path,
                           const DataSet& data,
                           const std::vector<DirectedPrimary>& primaries,
                           const AssignmentParameters& parameters,
                           AssignmentResult& result) {
    std::unordered_map<std::string, const DirectedPrimary*> forward;
    std::unordered_map<std::string, const DirectedPrimary*> reverse;
    std::unordered_map<std::string, std::vector<std::string>> forward_by_anchor_a;
    std::unordered_map<std::string, std::vector<std::string>> reverse_by_anchor_b;

    for (const auto& primary : primaries) {
        if (primary.target_sample == primary.query_sample) continue;
        const auto sample_a =
            std::min(primary.target_sample, primary.query_sample);
        const auto sample_b =
            std::max(primary.target_sample, primary.query_sample);
        const auto canonical = canonical_primary(primary);
        const auto key =
            merge_key(sample_a, sample_b, canonical.anchor_a, canonical.anchor_b);
        if (canonical.forward) {
            forward.emplace(key, &primary);
            forward_by_anchor_a[anchor_key(sample_a, sample_b, canonical.anchor_a)]
                .push_back(canonical.anchor_b);
        } else {
            reverse.emplace(key, &primary);
            reverse_by_anchor_b[anchor_key(sample_a, sample_b, canonical.anchor_b)]
                .push_back(canonical.anchor_a);
        }
    }

    std::vector<std::string> keys;
    keys.reserve(forward.size() + reverse.size());
    for (const auto& item : forward) keys.push_back(item.first);
    for (const auto& item : reverse) {
        if (forward.find(item.first) == forward.end()) keys.push_back(item.first);
    }
    std::sort(keys.begin(), keys.end());

    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("cannot create correspondence output");
    output << "sample_a\tsample_b\tanchor_a\tanchor_b\tsupport_direction"
              "\tforward_pair_id\treverse_pair_id\tforward_score\treverse_score"
              "\tforward_edit_distance\treverse_edit_distance"
              "\tforward_strand\treverse_strand"
              "\tforward_method\treverse_method\tdiscordant_anchor\n";

    for (const auto& key : keys) {
        const auto forward_found = forward.find(key);
        const auto reverse_found = reverse.find(key);
        const auto* f =
            forward_found == forward.end() ? nullptr : forward_found->second;
        const auto* r =
            reverse_found == reverse.end() ? nullptr : reverse_found->second;
        const auto& primary = f ? *f : *r;
        const auto sample_a =
            std::min(primary.target_sample, primary.query_sample);
        const auto sample_b =
            std::max(primary.target_sample, primary.query_sample);
        const auto canonical = canonical_primary(primary);

        std::string support_direction;
        std::string discordant_anchor = ".";
        if (f && r) {
            support_direction = "both";
            ++result.reciprocal_correspondence_count;
        } else if (f) {
            const auto opposite = reverse_by_anchor_b.find(
                anchor_key(sample_a, sample_b, canonical.anchor_b));
            if (opposite != reverse_by_anchor_b.end() &&
                std::find(opposite->second.begin(), opposite->second.end(),
                          canonical.anchor_a) == opposite->second.end()) {
                support_direction = "discordant";
                discordant_anchor = opposite->second.front();
                ++result.discordant_correspondence_count;
            } else {
                support_direction = "forward_only";
                ++result.forward_only_correspondence_count;
            }
        } else {
            const auto opposite = forward_by_anchor_a.find(
                anchor_key(sample_a, sample_b, canonical.anchor_a));
            if (opposite != forward_by_anchor_a.end() &&
                std::find(opposite->second.begin(), opposite->second.end(),
                          canonical.anchor_b) == opposite->second.end()) {
                support_direction = "discordant";
                discordant_anchor = opposite->second.front();
                ++result.discordant_correspondence_count;
            } else {
                support_direction = "reverse_only";
                ++result.reverse_only_correspondence_count;
            }
        }

        if (parameters.pair_merge_mode == PairMergeMode::reciprocal &&
            support_direction != "both") {
            if (support_direction == "forward_only")
                --result.forward_only_correspondence_count;
            else if (support_direction == "reverse_only")
                --result.reverse_only_correspondence_count;
            else if (support_direction == "discordant")
                --result.discordant_correspondence_count;
            continue;
        }

        output << data.samples[sample_a] << '\t' << data.samples[sample_b]
               << '\t' << canonical.anchor_a << '\t' << canonical.anchor_b
               << '\t' << support_direction << '\t'
               << (f ? std::to_string(f->pair_id) : ".") << '\t'
               << (r ? std::to_string(r->pair_id) : ".") << '\t'
               << (f ? std::to_string(f->score) : ".") << '\t'
               << (r ? std::to_string(r->score) : ".") << '\t'
               << (f ? std::to_string(f->edit_distance) : ".") << '\t'
               << (r ? std::to_string(r->edit_distance) : ".") << '\t'
               << (f ? std::string(1, f->strand) : ".") << '\t'
               << (r ? std::string(1, r->strand) : ".") << '\t'
               << (f ? f->method : ".") << '\t'
               << (r ? r->method : ".") << '\t' << discordant_anchor << '\n';
        ++result.correspondence_count;
    }

    write_correspondence_metadata(metadata_path, parameters, result);
}

}  // namespace vallescope2
