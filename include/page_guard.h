#ifndef INCLUDE_PAGE_GUARD
#define INCLUDE_PAGE_GUARD

#include "buffer_manager.h"

namespace guidedresearch {

class PageGuard {
protected:
    BufferManager& bm;
    BufferFrame *page;

    PageGuard(BufferManager &bm, BufferFrame *page) : bm(bm), page(page) {}
public:
    PageGuard() = delete;
    PageGuard(const PageGuard& other) = delete;
    BufferFrame* operator->() { assert(page); return page; }
    BufferFrame* bf_ptr() const { return page; }
};
}
#endif