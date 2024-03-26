#ifndef INCLUDE_PAGE_GUARD
#define INCLUDE_PAGE_GUARD

#include "buffer_manager/buffer_manager.h"

namespace guidedresearch {

/// A RAII wrapper for a BufferFrame.
template<bool Exclusive = false>
class PageGuard {
private:
    bool dirty;
protected:
    BufferManager& bm;
    BufferFrame *page;

    PageGuard(BufferManager &bm, BufferFrame *page, bool dirty = false) : dirty(dirty), bm(bm), page(page) {}
    void reset() { page = nullptr; if constexpr (Exclusive) dirty = false; }
    void unfix(int) { if (page) bm.unfix_page(*page, dirty); }
public:
    PageGuard() = delete;
    PageGuard(const PageGuard& other) = delete;
    PageGuard(PageGuard&& other) : PageGuard(other.bm, other.page, other.dirty) {
        other.reset();
    }
    PageGuard(BufferManager &bm) :PageGuard(bm, nullptr) {}
    PageGuard(BufferManager &bm, Swip &swip) : PageGuard(bm, &bm.fix_page(swip, Exclusive)) {}
    ~PageGuard() { unfix(); }

    PageGuard& operator=(const PageGuard& other) = delete;
    PageGuard& operator=(PageGuard&& other) {
        assert(&bm == &other.bm);
        if (this != &other) {
            unfix(0);
            page = other.page;
            if constexpr (Exclusive) dirty = other.dirty;
            other.reset();
        }
        return *this;
    }

    void mark_dirty() { if constexpr(Exclusive) dirty = true; }
    void fix(Swip& swip) { unfix(0); page = &bm.fix_page(swip, Exclusive); if constexpr (Exclusive) dirty = false; }
    void unfix() { unfix(0); reset(); }

    BufferFrame* operator->() { assert(page); return page; }
    BufferFrame* bf_ptr() const { return page; }
};

using SharedPage = PageGuard<false>;
using UniquePage = PageGuard<true>;

}
#endif