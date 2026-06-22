#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace vallescope2 {

std::string format_command(const std::vector<std::string>& arguments);

void run_process(const std::vector<std::string>& arguments,
                 const std::filesystem::path& log_path);

std::string capture_process(const std::vector<std::string>& arguments);

}  // namespace vallescope2
