#include "buffer_manager.h"
#include "swip.h"
#include "btree.h"
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
}  // namespace 