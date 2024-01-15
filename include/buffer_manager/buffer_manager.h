#ifndef INCLUDE_BUFFER_MANAGER_H
#define INCLUDE_BUFFER_MANAGER_H

#include <cstddef>
#include <cstdint>
#include <exception>
#include <unordered_map>
#include <vector>
#include <list>
#include <shared_mutex>
#include <mutex>
#include "buffer_manager/swip.h"

namespace guidedresearch {

class BufferFrame {
   private:
   friend class BufferManager;

   u64 pid;
   std::shared_mutex latch;
   bool exclusive = false;
   bool dirty = false;

   std::vector<char> data;

   public:
   /// Constructor.
   BufferFrame(u64 pid) noexcept;
   // Constructor.
   BufferFrame(BufferFrame&& o) noexcept;
   // Assignment
   BufferFrame& operator=(BufferFrame&& o) noexcept;

   /// Returns a pointer to this page's data.
   char* get_data() { return data.data(); }
};


class BufferManager {
   private:
   std::mutex bm_lock;
   size_t page_size, page_count;
   std::vector<BufferFrame> frames;

   // cold path
   void fix_page_slow(Swip& swip);

   public:
   BufferManager(const BufferManager&) = delete;
   BufferManager(BufferManager&&) = delete;
   BufferManager& operator=(const BufferManager&) = delete;
   BufferManager& operator=(BufferManager&&) = delete;

   /// Constructor.
   /// @param[in] page_size  Size in bytes that all pages will have.
   /// @param[in] page_count Maximum number of pages that should reside in
   //                        memory at the same time.
   BufferManager(size_t page_size, size_t page_count);

   /// Destructor. Writes all dirty pages to disk.
   ~BufferManager();

   /// Returns size of a page
   [[nodiscard]] size_t get_page_size() const { return page_size; }

   /// Returns a reference to a `BufferFrame` object for a given page id. When
   /// the page is not loaded into memory, it is read from disk. Otherwise the
   /// loaded page is used.
   /// When the page cannot be loaded because the buffer is full, throws the
   /// exception `buffer_full_error`.
   /// Is thread-safe w.r.t. other concurrent calls to `fix_page()` and
   /// `unfix_page()`.
   /// @param[in] swip  Swip that holds a reference to the page
   /// @param[in] exclusive If `exclusive` is true, the page is locked
   ///                      exclusively. Otherwise it is locked
   ///                      non-exclusively (shared).
   inline BufferFrame& fix_page(Swip& swip, bool exclusive) {
      if (swip.isUNSWIZZLED()) {
         fix_page_slow(swip);
      }
      auto &bf = swip.asBufferFrame();
      if (exclusive) {
         bf.latch.lock();
         bf.exclusive = true;
      }
      else {
         bf.latch.lock_shared();
      }
      return bf;
   }

   /// Takes a `BufferFrame` reference that was returned by an earlier call to
   /// `fix_page()` and unfixes it. When `is_dirty` is / true, the page is
   /// written back to disk eventually.
   void unfix_page(BufferFrame& page, bool is_dirty);

   /// Returns a list of all page ids that are currently loaded in memory.
   /// For debugging purposes, not thread-safe.
   std::vector<u64> get_all_pages() const;

   /// Returns the segment id for a given page id which is contained in the 16
   /// most significant bits of the page id.
   static constexpr uint16_t get_segment_id(u64 page_id) {
      return page_id >> 48;
   }

   /// Returns the page id within its segment for a given page id. This
   /// corresponds to the 48 least significant bits of the page id.
   static constexpr uint64_t get_segment_page_id(u64 page_id) {
      return page_id & ((1ull << 48) - 1);
   }

   static constexpr uint64_t get_page_id(uint16_t segment_id, uint64_t segment_page_id) {
      return (static_cast<uint64_t>(segment_id) << 48) ^ segment_page_id;
   }
};

class buffer_full_error
: public std::exception {
public:
    [[nodiscard]] const char* what() const noexcept override {
        return "buffer is full";
    }
};
} // namespace guidedresearch

#endif
