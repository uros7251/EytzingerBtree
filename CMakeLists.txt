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
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w")
add_compile_options(-Wall -Wextra -Wvla)

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------

find_package(Threads REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)

include("${CMAKE_SOURCE_DIR}/vendor/googletest.cmake")

# ---------------------------------------------------------------------------
# Includes
# ---------------------------------------------------------------------------

include("${CMAKE_SOURCE_DIR}/include/local.cmake")
include_directories(
    ${CMAKE_SOURCE_DIR}/include
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

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

include("${CMAKE_SOURCE_DIR}/test/local.cmake")