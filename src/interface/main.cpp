#include "vallescope2/anchor/anchor_detection.hpp"
#include "vallescope2/interface/options.hpp"

#include <exception>
#include <iostream>

int main(const int argc, char* argv[]) {
    try {
        vallescope2::run_anchor_detection(vallescope2::parse_arguments(argc, argv));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "vallescope2: error: " << error.what() << '\n';
        return 1;
    }
}
