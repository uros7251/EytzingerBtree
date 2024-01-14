# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

set(TEST_CC test/bm_test.cc test/btree_test.cc test/eytzinger_test.cc)

# ---------------------------------------------------------------------------
# Tester
# ---------------------------------------------------------------------------

add_executable(tester test/tester.cc ${TEST_CC})
target_link_libraries(tester grlib gtest gmock Threads::Threads)

enable_testing()
add_test(guided-research tester)