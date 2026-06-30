cmake_minimum_required(VERSION 3.16)

project(ValleScope2 VERSION 0.0.0 LANGUAGES CXX)

include(CTest)

find_package(OpenSSL REQUIRED)

add_subdirectory(external/WFA2-lib EXCLUDE_FROM_ALL)

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
    src/alignment/alignment_engine.cpp
    src/alignment/anchor_guided_alignment.cpp
    src/alignment/base_alignment.cpp
    src/alignment/base_alignment_io.cpp
    src/alignment/base_alignment_output.cpp
    src/alignment/bundle_postprocessing.cpp
    src/anchor/anchor_detection.cpp
    src/anchor/anchor_output.cpp
    src/anchor/anchor_selection.cpp
    src/anchor/genmap_runner.cpp
    src/anchor/window_scoring.cpp
    src/bundle/chaining.cpp
    src/bundle/chaining_dp.cpp
    src/bundle/chaining_io.cpp
    src/bundle/chaining_refinement.cpp
    src/bundle/chaining_trim.cpp
    src/common/process.cpp
    src/common/sha256.cpp
    src/common/workspace.cpp
    src/context/anchor_grouping.cpp
    src/context/structural_tokens.cpp
    src/context/structural_context.cpp
    src/correspondence/assignment.cpp
    src/correspondence/assignment_common.cpp
    src/correspondence/assignment_loader.cpp
    src/correspondence/assignment_matching.cpp
    src/correspondence/assignment_merge.cpp
    src/correspondence/assignment_output.cpp
    src/interface/main.cpp
    src/interface/options.cpp
    src/input/fasta.cpp
    src/input/fasta_index.cpp
    src/input/input_preparation.cpp
    src/input/sequence_catalog.cpp
)

target_include_directories(vallescope2 PRIVATE include external/WFA2-lib)
target_compile_features(vallescope2 PRIVATE cxx_std_17)
target_compile_options(vallescope2 PRIVATE -Wall -Wextra -Wpedantic)
target_compile_definitions(vallescope2 PRIVATE
    VALLESCOPE2_VERSION="${PROJECT_VERSION}"
)
target_link_libraries(vallescope2 PRIVATE
    HTSlib::HTSlib OpenSSL::Crypto wfa2::wfa2cpp_static
)

install(TARGETS vallescope2 RUNTIME DESTINATION bin)

if(BUILD_TESTING)
    add_executable(anchor_grouping_test
        test/anchor_grouping_test.cpp
        src/alignment/alignment_engine.cpp
        src/alignment/anchor_guided_alignment.cpp
        src/alignment/base_alignment.cpp
        src/alignment/base_alignment_io.cpp
        src/alignment/base_alignment_output.cpp
        src/alignment/bundle_postprocessing.cpp
        src/bundle/chaining.cpp
        src/bundle/chaining_dp.cpp
        src/bundle/chaining_io.cpp
        src/bundle/chaining_refinement.cpp
        src/bundle/chaining_trim.cpp
        src/common/sha256.cpp
        src/context/anchor_grouping.cpp
        src/context/structural_tokens.cpp
        src/context/structural_context.cpp
        src/correspondence/assignment.cpp
        src/correspondence/assignment_common.cpp
        src/correspondence/assignment_loader.cpp
        src/correspondence/assignment_matching.cpp
        src/correspondence/assignment_merge.cpp
        src/correspondence/assignment_output.cpp
        src/input/fasta_index.cpp
    )
    target_include_directories(anchor_grouping_test PRIVATE include external/WFA2-lib)
    target_compile_features(anchor_grouping_test PRIVATE cxx_std_17)
    target_compile_options(anchor_grouping_test PRIVATE -Wall -Wextra -Wpedantic)
    target_link_libraries(anchor_grouping_test PRIVATE
        HTSlib::HTSlib OpenSSL::Crypto wfa2::wfa2cpp_static
    )
    add_test(NAME anchor-grouping COMMAND anchor_grouping_test)

    add_executable(edit_distance_test
        test/edit_distance_test.cpp
    )
    target_include_directories(edit_distance_test PRIVATE include)
    target_compile_features(edit_distance_test PRIVATE cxx_std_17)
    target_compile_options(edit_distance_test PRIVATE -Wall -Wextra -Wpedantic)
    add_test(NAME edit-distance COMMAND edit_distance_test)
endif()
