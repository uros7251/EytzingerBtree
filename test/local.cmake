# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

set(TEST_CC test/bm_test.cc)

# ---------------------------------------------------------------------------
# Tester
# ---------------------------------------------------------------------------

add_executable(tester test/tester.cc ${TEST_CC})
target_link_libraries(tester grlib gtest gmock Threads::Threads)

enable_testing()
add_test(moderndbs tester)