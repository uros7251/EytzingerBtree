#include "buffer_manager.h"
#include <cassert>

/*
This is only a dummy implementation of a buffer manager. It does not do any
disk I/O and does not respect the page_count.
*/


namespace guidedresearch {
BufferFrame::BufferFrame(u64 pid, u8* data) noexcept : pid(pid), dirty(false), data(data) {

}

BufferFrame::BufferFrame(BufferFrame&& o) noexcept : pid(o.pid), dirty(o.dirty), data(o.data) {

}

BufferFrame& BufferFrame::operator=(BufferFrame&& o) noexcept {
   pid = o.pid;
   dirty = o.dirty;
   data = o.data;
   return *this;
};


void BufferManager::unswizzle_random() {                                        
   auto* frame = &frames[0]; // choose a random hot page instead
   bool found;
   while (true) { // loop until above constraint is satisfied
      found = true;
      for (auto &swip : frame->get_swips()) {
         if (swip.isUNSWIZZLED()) {
            frame = &swip.asBufferFrame();
            found = false;
            break;
         }
      }
      if (found) break;
   }
   // unswizzle
}

BufferManager::BufferManager(size_t page_size, size_t page_count) 
   : page_size(page_size), page_count(page_count)
{
   buffer = new u8[page_size*page_count];
}

BufferManager::~BufferManager() {
   delete[] buffer;
}


BufferFrame& BufferManager::fix_page(Swip& swip, bool exclusive) {
   if (swip.isSWIZZLED()) {
      // if the page is swizzled, the swip holds a pointer
      return swip.asBufferFrame();
   }
   // cold access, take the global lock
   std::lock_guard lock {bm_lock};
   // check if the page is in the cooling stage
   auto pid = swip.asPageID();
   BufferFrame* bf { nullptr };
   if (auto hm_it = hashmap.find(pid); hm_it != hashmap.end()) {
      // page is in memory, swizzle it back
      auto f_it = hm_it->second; // get an iterator to the element of the fifo queue
      bf = *f_it;
      swip.warm(bf);
      hashmap.erase(hm_it);
      fifo.erase(f_it);
   }
   else {
      // page is not in memory
      // for now, just allocate a new frame in memory
      bf = &frames.emplace_back(pid, &buffer[page_size * frames.size()]);
   }
   if (too_few_cooling()) {
      unswizzle_random();
   }
   if (exclusive) {
      bf->latch.lock();
   }
   return *bf;
}


void BufferManager::unfix_page(BufferFrame& page, bool is_dirty) {
   if (is_dirty) {
      page.dirty = true;
      page.latch.unlock();
   }
}

}  // namespace guidedresearch
