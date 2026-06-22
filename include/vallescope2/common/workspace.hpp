#pragma once

#include <filesystem>

namespace vallescope2 {

class Workspace {
public:
    explicit Workspace(bool keep);
    ~Workspace();
    Workspace(const Workspace&) = delete;
    Workspace& operator=(const Workspace&) = delete;

    const std::filesystem::path& path() const { return path_; }

private:
    bool keep_;
    std::filesystem::path path_;
};

}  // namespace vallescope2
