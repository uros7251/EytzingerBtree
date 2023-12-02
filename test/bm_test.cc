#include "buffer_manager.h"
#include "Swip.h"
#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstring>
#include <gtest/gtest.h>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

namespace {

using BufferManager = guidedresearch::BufferManager;
using BufferFrame = guidedresearch::BufferFrame;
using buffer_full_error = guidedresearch::buffer_full_error;
using Swip = guidedresearch::Swip;

// NOLINTNEXTLINE
TEST(BufferManagerTest, FixSingle) {
    BufferManager buffer_manager{1024, 10};
    auto swip = Swip::fromPID(1ull);
    std::vector<uint64_t> expected_values(1024 / sizeof(uint64_t), 123);
    {
        auto& page = buffer_manager.fix_page(swip, true);
        ASSERT_EQ(&swip.asBufferFrame(), &page); // check if the swip is correctly swizzled
        ASSERT_TRUE(page.get_data());
        std::memcpy(page.get_data(), expected_values.data(), 1024);
        buffer_manager.unfix_page(page, true);
    }
    {
        std::vector<uint64_t> values(1024 / sizeof(uint64_t));
        auto& page = buffer_manager.fix_page(swip, false);
        std::memcpy(values.data(), page.get_data(), 1024);
        buffer_manager.unfix_page(page, true);
        ASSERT_EQ(expected_values, values);
    }
}

// NOLINTNEXTLINE
TEST(BufferManagerTest, MultithreadParallelFix) {
    BufferManager buffer_manager{1024, 10};
    std::vector<std::thread> threads;
    std::vector<Swip> swips(8);
    for (size_t i = 0; i < 8; ++i) {
        swips[i] = Swip::fromPID(i);
    }
    for (size_t i = 0; i < 4; ++i) {
        threads.emplace_back([i, &buffer_manager, &swips] {
            // NOLINTNEXTLINE
            ASSERT_NO_THROW(
                auto& page1 = buffer_manager.fix_page(swips[i], false);
                auto& page2 = buffer_manager.fix_page(swips[i+4], false);
                buffer_manager.unfix_page(page1, false);
                buffer_manager.unfix_page(page2, false);
            );
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    auto pages = buffer_manager.get_all_pages();
    std::sort(pages.begin(), pages.end());
    std::vector<uint64_t> expected_pages{0, 1, 2, 3, 4, 5, 6, 7};
    EXPECT_EQ(expected_pages, pages);
}


// NOLINTNEXTLINE
TEST(BufferManagerTest, MultithreadExclusiveAccess) {
    BufferManager buffer_manager{1024, 10};
    auto swip = Swip::fromPID(0ull);
    {
        auto& page = buffer_manager.fix_page(swip, true);
        ASSERT_TRUE(page.get_data());
        std::memset(page.get_data(), 0, 1024);
        buffer_manager.unfix_page(page, true);
    }
    std::vector<std::thread> threads;
    for (size_t i = 0; i < 4; ++i) {
        threads.emplace_back([&buffer_manager, &swip] {
            for (size_t j = 0; j < 1000; ++j) {
                auto& page = buffer_manager.fix_page(swip, true);
                ASSERT_TRUE(page.get_data());
                uint64_t& value = *reinterpret_cast<uint64_t*>(page.get_data());
                ++value;
                buffer_manager.unfix_page(page, true);
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    EXPECT_EQ(std::vector<uint64_t>{0}, buffer_manager.get_all_pages());
    auto& page = buffer_manager.fix_page(swip, false);
    ASSERT_TRUE(page.get_data());
    uint64_t value = *reinterpret_cast<uint64_t*>(page.get_data());
    buffer_manager.unfix_page(page, false);
    EXPECT_EQ(4000, value);
}

// NOLINTNEXTLINE
TEST(BufferManagerTest, BlockedThreadsHoldsNoLocks) {
  BufferManager buffer_manager{1024, 10};
  auto swip_0 = Swip::fromPID(0ull), swip_1 = Swip::fromPID(1ull);
  for (size_t i = 0; i < 2; ++i ){
    auto& page = buffer_manager.fix_page(i==0 ? swip_0 : swip_1, true);
    ASSERT_TRUE(page.get_data());
    std::memset(page.get_data(), 0, 1024);
    buffer_manager.unfix_page(page, true);
  }
  auto blockedCv = std::condition_variable();
  auto blockingCv = std::condition_variable();
  auto cvMutex = std::mutex();
  std::atomic<bool> isBlocked = false;
  auto blocking = std::thread([&] {
    auto lock = std::unique_lock(cvMutex);
    auto& page = buffer_manager.fix_page(swip_0, true);
    isBlocked = true;
    blockedCv.notify_all();
    uint64_t& value = *reinterpret_cast<uint64_t*>(page.get_data());
    ++value;
    blockingCv.wait(lock);
    ASSERT_EQ(value, 1);
    buffer_manager.unfix_page(page, true);
  });
  auto blocked = std::thread([&] {
    auto lock = std::unique_lock(cvMutex);
    while (!isBlocked) {
      blockedCv.wait(lock);
    }
    lock.unlock();

    auto& page = buffer_manager.fix_page(swip_0, true);
    uint64_t& value = *reinterpret_cast<uint64_t*>(page.get_data());
    ASSERT_EQ(value, 1);
    buffer_manager.unfix_page(page, false);
  });

  {
    auto lock = std::unique_lock(cvMutex);
    while (!isBlocked) {
      blockedCv.wait(lock);
    }
  }

  std::vector<std::thread> threads;
  for (size_t i = 0; i < 4; ++i) {
    threads.emplace_back([&buffer_manager, &swip_1] {
      for (size_t j = 0; j < 1000; ++j) {
        // Threads operating on page 1 should not be blocked by the blocking thread opn page 0
        auto& page = buffer_manager.fix_page(swip_1, false);
        ASSERT_TRUE(page.get_data());
        uint64_t& value = *reinterpret_cast<uint64_t*>(page.get_data());
        ASSERT_EQ(value, 0);
        buffer_manager.unfix_page(page, false);
      }
    });
  }
  for (auto& thread : threads) {
    thread.join();
  }

  blockingCv.notify_all();
  blocking.join();
  blocked.join();
}

// NOLINTNEXTLINE
TEST(BufferManagerTest, MultithreadBufferFull) {
    BufferManager buffer_manager{1024, 10};
    std::atomic<uint64_t> num_buffer_full = 0;
    std::atomic<uint64_t> finished_threads = 0;
    std::vector<Swip> swips(16);
    for (size_t i = 0; i < 16; ++i) {
        swips[i] = Swip::fromPID(i);
    }
    std::vector<std::thread> threads;
    for (size_t i = 0; i < 4; ++i) {
        threads.emplace_back([i, &buffer_manager, &num_buffer_full, &finished_threads, &swips] {
            std::vector<BufferFrame*> pages;
            pages.reserve(4);
            for (size_t j = 0; j < 4; ++j) {
                try {
                    pages.push_back(&buffer_manager.fix_page(swips[i + j * 4], false));
                } catch (const buffer_full_error&) {
                    ++num_buffer_full;
                }
            }
            ++finished_threads;
            // Busy wait until all threads have finished.
            while (finished_threads.load() < 4) {}
            for (auto* page : pages) {
                buffer_manager.unfix_page(*page, false);
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    EXPECT_EQ(10, buffer_manager.get_all_pages().size());
    EXPECT_EQ(6, num_buffer_full.load());
}


// NOLINTNEXTLINE
/*TEST(BufferManagerTest, MultithreadManyPages) {
    BufferManager buffer_manager{1024, 10};
    std::vector<std::thread> threads;
    for (size_t i = 0; i < 4; ++i) {
        threads.emplace_back([i, &buffer_manager] {
            std::mt19937_64 engine{i};
            std::geometric_distribution<uint64_t> distr{0.1};
            for (size_t j = 0; j < 10000; ++j) {
                ASSERT_NO_THROW(
                    auto& page = buffer_manager.fix_page(distr(engine), false);
                    buffer_manager.unfix_page(page, false);
                );
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
}


// NOLINTNEXTLINE
TEST(BufferManagerTest, MultithreadReaderWriter) {
    {
        // Zero out all pages first
        BufferManager buffer_manager{1024, 10};
        for (uint16_t segment = 0; segment <= 3; ++segment) {
            for (uint64_t segment_page = 0; segment_page <= 100; ++segment_page) {
                uint64_t page_id = (static_cast<uint64_t>(segment) << 48) | segment_page;
                auto& page = buffer_manager.fix_page(page_id, true);
                ASSERT_TRUE(page.get_data());
                std::memset(page.get_data(), 0, 1024);
                buffer_manager.unfix_page(page, true);
            }
        }
        // Let the buffer manager be destroyed here so that the caches are
        // empty before running the actual test.
    }
    BufferManager buffer_manager{1024, 10};
    std::atomic<size_t> aborts = 0;
    std::vector<std::thread> threads;
    for (size_t i = 0; i < 4; ++i) {
        threads.emplace_back([i, &buffer_manager, &aborts] {
            std::mt19937_64 engine{i};
            // 5% of queries are scans.
            std::bernoulli_distribution scan_distr{0.05};
            // Number of pages accessed by a point query is geometrically
            // distributed.
            std::geometric_distribution<size_t> num_pages_distr{0.5};
            // 60% of point queries are reads.
            std::bernoulli_distribution reads_distr{0.6};
            // Out of 20 accesses, 12 are from segment 0, 5 from segment 1,
            // 2 from segment 2, and 1 from segment 3.
            std::discrete_distribution<uint16_t> segment_distr{12.0, 5.0, 2.0, 1.0};
            // Page accesses for point queries are uniformly distributed in
            // [0, 100].
            std::uniform_int_distribution<uint64_t> page_distr{0, 100};
            std::vector<uint64_t> scan_sums(4);
            for (size_t j = 0; j < 100; ++j) {
                uint16_t segment = segment_distr(engine);
                uint64_t segment_shift = static_cast<uint64_t>(segment) << 48;
                if (scan_distr(engine)) {
                    // scan
                    uint64_t scan_sum = 0;
                    for (uint64_t segment_page = 0; segment_page <= 100; ++segment_page) {
                        uint64_t page_id = segment_shift | segment_page;
                        BufferFrame* page = nullptr;
                        while (true) {
                            try {
                                page = &buffer_manager.fix_page(page_id, false);
                                break;
                            } catch (const buffer_full_error&) {
                                // Don't abort scan when the buffer is full, retry
                                // the current page.
                            }
                        }
                        ASSERT_TRUE(page->get_data());
                        uint64_t value = *reinterpret_cast<uint64_t*>(page->get_data());
                        scan_sum += value;
                        buffer_manager.unfix_page(*page, false);
                    }
                    EXPECT_GE(scan_sum, scan_sums[segment]);
                    scan_sums[segment] = scan_sum;
                } else {
                    // point query
                    auto num_pages = num_pages_distr(engine) + 1;
                    // For point queries all accesses but the last are always
                    // reads. Only the last is potentially a write. Also,
                    // all pages but the last are held for the entire duration
                    // of the query.
                    std::vector<BufferFrame*> pages;
                    auto unfix_pages = [&] {
                        for (auto it = pages.rbegin(); it != pages.rend(); ++it) {
                            auto& page = **it;
                            buffer_manager.unfix_page(page, false);
                        }
                        pages.clear();
                    };
                    for (size_t page_number = 0; page_number < num_pages - 1; ++page_number) {
                        uint64_t segment_page = page_distr(engine);
                        uint64_t page_id = segment_shift | segment_page;
                        BufferFrame* page = nullptr;
                        try {
                            page = &buffer_manager.fix_page(page_id, false);
                        } catch (const buffer_full_error&) {
                            // Abort query when buffer is full.
                            ++aborts;
                            goto abort;
                        }
                        pages.push_back(page);
                    }
                    // Unfix all pages before accessing the last one
                    // (potentially exclusively) to avoid deadlocks.
                    unfix_pages();
                    {
                        uint64_t segment_page = page_distr(engine);
                        uint64_t page_id = segment_shift | segment_page;
                        if (reads_distr(engine)) {
                            // read
                            BufferFrame* page = nullptr;
                            try {
                                page = &buffer_manager.fix_page(page_id, false);
                            } catch (const buffer_full_error&) {
                                ++aborts;
                                goto abort;
                            }
                            buffer_manager.unfix_page(*page, false);
                        } else {
                            // write
                            BufferFrame* page = nullptr;
                            try {
                                page = &buffer_manager.fix_page(page_id, true);
                            } catch (const buffer_full_error&) {
                                ++aborts;
                                goto abort;
                            }
                            ASSERT_TRUE(page->get_data());
                            auto& value = *reinterpret_cast<uint64_t*>(page->get_data());
                            ++value;
                            buffer_manager.unfix_page(*page, true);
                        }
                    }
                    abort:
                    unfix_pages();
                }
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    EXPECT_LT(aborts.load(), 20);
}
*/
}  // namespace 