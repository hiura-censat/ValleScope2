cmake_minimum_required(VERSION 3.16)

project(ValleScope2 VERSION 0.0.0 LANGUAGES CXX)

find_package(OpenSSL REQUIRED)

set(HTSLIB_ROOT "" CACHE PATH "Root directory of the HTSlib installation")
find_path(HTSLIB_INCLUDE_DIR
    NAMES htslib/faidx.h
    HINTS ${HTSLIB_ROOT} $ENV{CONDA_PREFIX}
    PATH_SUFFIXES include
    REQUIRED
)
find_library(HTSLIB_LIBRARY
    NAMES hts
    HINTS ${HTSLIB_ROOT} $ENV{CONDA_PREFIX}
    PATH_SUFFIXES lib
    REQUIRED
)

add_library(HTSlib::HTSlib UNKNOWN IMPORTED)
set_target_properties(HTSlib::HTSlib PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}"
)

add_executable(vallescope2
    src/anchor/anchor_detection.cpp
    src/anchor/anchor_output.cpp
    src/anchor/anchor_selection.cpp
    src/anchor/genmap_runner.cpp
    src/anchor/window_scoring.cpp
    src/common/process.cpp
    src/common/workspace.cpp
    src/interface/main.cpp
    src/interface/options.cpp
    src/input/fasta.cpp
    src/input/fasta_index.cpp
    src/input/input_preparation.cpp
    src/input/sequence_catalog.cpp
)

target_include_directories(vallescope2 PRIVATE include)
target_compile_features(vallescope2 PRIVATE cxx_std_17)
target_compile_options(vallescope2 PRIVATE -Wall -Wextra -Wpedantic)
target_compile_definitions(vallescope2 PRIVATE
    VALLESCOPE2_VERSION="${PROJECT_VERSION}"
)
target_link_libraries(vallescope2 PRIVATE HTSlib::HTSlib OpenSSL::Crypto)

install(TARGETS vallescope2 RUNTIME DESTINATION bin)
