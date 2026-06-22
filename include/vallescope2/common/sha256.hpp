#pragma once

#include <filesystem>
#include <string>
#include <string_view>

namespace vallescope2 {

std::string sha256_file(const std::filesystem::path& path);
std::string sha256_text(std::string_view text);

}  // namespace vallescope2
