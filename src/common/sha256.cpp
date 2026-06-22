#include "vallescope2/common/sha256.hpp"

#include <openssl/evp.h>

#include <array>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace vallescope2 {
namespace {

using Context = std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)>;

Context make_context() {
    Context context(EVP_MD_CTX_new(), EVP_MD_CTX_free);
    if (!context || EVP_DigestInit_ex(context.get(), EVP_sha256(), nullptr) != 1) {
        throw std::runtime_error("failed to initialize SHA-256");
    }
    return context;
}

std::string finish(Context& context) {
    std::array<unsigned char, EVP_MAX_MD_SIZE> digest{};
    unsigned int digest_size = 0;
    if (EVP_DigestFinal_ex(context.get(), digest.data(), &digest_size) != 1) {
        throw std::runtime_error("failed to finalize SHA-256");
    }
    std::ostringstream value;
    value << std::hex << std::setfill('0');
    for (unsigned int i = 0; i < digest_size; ++i) {
        value << std::setw(2) << static_cast<unsigned int>(digest[i]);
    }
    return value.str();
}

}  // namespace

std::string sha256_file(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot open file for SHA-256: " + path.string());
    }
    auto context = make_context();
    std::array<char, 1 << 16> buffer{};
    while (input) {
        input.read(buffer.data(), buffer.size());
        const auto count = input.gcount();
        if (count > 0 &&
            EVP_DigestUpdate(context.get(), buffer.data(), count) != 1) {
            throw std::runtime_error("failed to update SHA-256");
        }
    }
    return finish(context);
}

std::string sha256_text(const std::string_view text) {
    auto context = make_context();
    if (EVP_DigestUpdate(context.get(), text.data(), text.size()) != 1) {
        throw std::runtime_error("failed to update SHA-256");
    }
    return finish(context);
}

}  // namespace vallescope2
