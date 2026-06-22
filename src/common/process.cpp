#include "vallescope2/common/process.hpp"

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include <array>
#include <cerrno>
#include <fstream>
#include <stdexcept>

namespace vallescope2 {
namespace {

std::vector<char*> make_argv(const std::vector<std::string>& arguments) {
    std::vector<char*> argv;
    argv.reserve(arguments.size() + 1);
    for (const auto& argument : arguments) {
        argv.push_back(const_cast<char*>(argument.c_str()));
    }
    argv.push_back(nullptr);
    return argv;
}

void require_success(const pid_t process, const std::string& message) {
    int status = 0;
    if (waitpid(process, &status, 0) < 0 || !WIFEXITED(status) ||
        WEXITSTATUS(status) != 0) {
        throw std::runtime_error(message);
    }
}

}  // namespace

std::string format_command(const std::vector<std::string>& arguments) {
    std::string text;
    for (const auto& argument : arguments) {
        if (!text.empty()) text += ' ';
        text += argument;
    }
    return text;
}

void run_process(const std::vector<std::string>& arguments,
                 const std::filesystem::path& log_path) {
    std::ofstream log(log_path, std::ios::app);
    if (!log) throw std::runtime_error("cannot open process log");
    log << "$ " << format_command(arguments) << '\n';
    log.close();

    const pid_t process = fork();
    if (process < 0) throw std::runtime_error("failed to start process");
    if (process == 0) {
        const int fd = open(log_path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644);
        if (fd < 0 || dup2(fd, STDOUT_FILENO) < 0 ||
            dup2(fd, STDERR_FILENO) < 0) {
            _exit(126);
        }
        close(fd);
        auto argv = make_argv(arguments);
        execvp(argv.front(), argv.data());
        _exit(errno == ENOENT ? 127 : 126);
    }
    require_success(process, "process failed; see log: " + log_path.string());
}

std::string capture_process(const std::vector<std::string>& arguments) {
    int descriptors[2];
    if (pipe(descriptors) != 0) {
        throw std::runtime_error("failed to create process pipe");
    }
    const pid_t process = fork();
    if (process < 0) {
        close(descriptors[0]);
        close(descriptors[1]);
        throw std::runtime_error("failed to start process");
    }
    if (process == 0) {
        close(descriptors[0]);
        if (dup2(descriptors[1], STDOUT_FILENO) < 0 ||
            dup2(descriptors[1], STDERR_FILENO) < 0) {
            _exit(126);
        }
        close(descriptors[1]);
        auto argv = make_argv(arguments);
        execvp(argv.front(), argv.data());
        _exit(127);
    }

    close(descriptors[1]);
    std::string output;
    std::array<char, 4096> buffer{};
    ssize_t count = 0;
    while ((count = read(descriptors[0], buffer.data(), buffer.size())) > 0) {
        output.append(buffer.data(), static_cast<std::size_t>(count));
    }
    close(descriptors[0]);
    require_success(process, "failed to capture process output");
    const auto newline = output.find('\n');
    return output.substr(0, newline);
}

}  // namespace vallescope2
