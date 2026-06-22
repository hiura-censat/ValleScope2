#include "vallescope2/common/workspace.hpp"

#include <chrono>

namespace vallescope2 {

Workspace::Workspace(const bool keep) : keep_(keep) {
    if (keep_) {
        path_ = std::filesystem::current_path() / "vallescope2-debug";
    } else {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path() /
                ("vallescope2-" + std::to_string(stamp));
    }
    std::filesystem::create_directories(path_);
}

Workspace::~Workspace() {
    if (!keep_) {
        std::error_code error;
        std::filesystem::remove_all(path_, error);
    }
}

}  // namespace vallescope2
