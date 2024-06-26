# ---------------------------------------------------------------------------
# GUIDED RESEARCH
# ---------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.7)
project(guided-research)

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -funroll-loops")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -funroll-loops")
add_compile_options(-Wall -Wextra -Wvla)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

find_package(Threads REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)

include("${CMAKE_SOURCE_DIR}/cmake/clang-tidy.cmake")
include("${CMAKE_SOURCE_DIR}/vendor/googletest.cmake")

# ---------------------------------------------------------------------------
# Includes
# ---------------------------------------------------------------------------

include("${CMAKE_SOURCE_DIR}/include/local.cmake")
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/btree
    ${CMAKE_SOURCE_DIR}/include/buffer_manager
    ${CMAKE_SOURCE_DIR}/include/utils
)

include_directories(SYSTEM
    ${GTEST_INCLUDE_DIR}
)

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

include("${CMAKE_SOURCE_DIR}/src/local.cmake")

set(SRC_CC
    "${CMAKE_SOURCE_DIR}/src/buffer_manager.cc")

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
add_library(grlib STATIC ${SRC_CC})
target_link_libraries(grlib Threads::Threads)

add_executable(bench "${CMAKE_SOURCE_DIR}/src/bench.cc" ${SRC_CC})
target_link_libraries(bench grlib)
target_compile_definitions(bench PRIVATE SINGLE_THREADED)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

include("${CMAKE_SOURCE_DIR}/test/local.cmake")

# ---------------------------------------------------------------------------
# Linting
# ---------------------------------------------------------------------------

add_clang_tidy_target(src_linting "${SRC_CC}")
add_clang_tidy_target(include_linting "${INCLUDE_H}")
add_clang_tidy_target(test_linting "${TEST_CC}")

add_custom_target(lint)
list(APPEND lint_targets include_linting src_linting test_linting)
add_dependencies(lint ${lint_targets})