#include "buffer_manager.h"
#include <cassert>

/*
This is only a dummy implementation of a buffer manager. It does not do any
disk I/O and does not respect the page_count.
*/


namespace guidedresearch {
BufferFrame::BufferFrame(u64 pid) noexcept : pid(pid) {

}

BufferFrame::BufferFrame(BufferFrame&& o) noexcept
   : pid(o.pid), exclusive(o.exclusive), dirty(o.dirty), data(move(o.data)) {

}

BufferFrame& BufferFrame::operator=(BufferFrame&& o) noexcept {
   pid = o.pid;
   exclusive = o.exclusive;
   dirty = o.dirty;
   data = move(o.data);
   return *this;
};

BufferManager::BufferManager(size_t page_size, size_t page_count) 
   : page_size(page_size), page_count(page_count)
{
   frames.reserve(page_count);
}

BufferManager::~BufferManager() = default;


BufferFrame& BufferManager::fix_page(Swip& swip, bool exclusive) {
   if (swip.isUNSWIZZLED()) {
      // cold access, take the global lock
      std::lock_guard<std::mutex> lock {bm_lock};
      // check again if the page was not swizzled in the meantime by another thread
      if (swip.isUNSWIZZLED()) {
         if (frames.size() >= page_count) {
            throw buffer_full_error();
         }
         auto &bf = frames.emplace_back(swip.asPageID());
         bf.data.resize(page_size, 0);
         swip.warm(&bf);
      }
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


void BufferManager::unfix_page(BufferFrame& page, bool is_dirty) {
   if (page.exclusive) {
      page.exclusive = false;
      page.dirty = is_dirty;
      page.latch.unlock();
   }
   else {
      page.latch.unlock_shared();
   }
}

std::vector<u64> BufferManager::get_all_pages() const
{
   std::vector<u64> result;
   for (auto &bf : frames) {
      result.push_back(bf.pid);
   }
   return result;
}

}  // namespace guidedresearch
