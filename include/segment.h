#ifndef INCLUDE_GUIDED_RESEARCH_SEGMENT_H_
#define INCLUDE_GUIDED_RESEARCH_SEGMENT_H_

#include "buffer_manager.h"
#include <array>
#include <optional>

namespace guidedresearch {

class Segment {
   public:
   /// Constructor.
   /// @param[in] segment_id       Id of the segment.
   /// @param[in] buffer_manager   The buffer manager that should be used by the segment.
   Segment(uint16_t segment_id, BufferManager& buffer_manager)
      : segment_id(segment_id), used_pages_count(0), buffer_manager(buffer_manager) {}

   /// The segment id
   uint16_t segment_id;
   /// The next free page;
   std::atomic_uint64_t used_pages_count;

   protected:
   /// The buffer manager
   BufferManager& buffer_manager;
   /// Pages to be reused
   std::forward_list<BufferFrame*> free_pages;
   /// Mutex for free_pages
   std::mutex free_pages_mutex;

   Swip allocate_page() {
      std::lock_guard lock(free_pages_mutex);
      if (free_pages.empty()) {
         Swip swip (buffer_manager.get_page_id(segment_id, used_pages_count++));
         return swip;
      }
      else {
         Swip swip = Swip(free_pages.front());
         free_pages.pop_front();
         return swip;
      }
   }

   void deallocate_page(BufferFrame *page) {
      std::lock_guard lock(free_pages_mutex);
      free_pages.push_front(page);
   }
};

} // namespace guidedresearch

#endif // INCLUDE_GUIDED_RESEARCH_SEGMENT_H_