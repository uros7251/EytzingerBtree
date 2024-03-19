#ifndef INCLUDE_PAGE_GUARD
#define INCLUDE_PAGE_GUARD

#include "buffer_manager/buffer_manager.h"

namespace guidedresearch {
    
template<bool Exclusive = false>
class PageGuard {
private:
    bool dirty;
protected:
    BufferManager& bm;
    BufferFrame *page;

    PageGuard(BufferManager &bm, BufferFrame *page, bool dirty = false) : dirty(dirty), bm(bm), page(page) {}
    void reset() { page = nullptr; if constexpr (Exclusive) dirty = false; }
public:
    PageGuard() = delete;
    PageGuard(const PageGuard& other) = delete;
    PageGuard(PageGuard&& other) : PageGuard(other.bm, other.page, other.dirty) {
        other.reset();
    }
    PageGuard(BufferManager &bm) :PageGuard(bm, nullptr) {}
    PageGuard(BufferManager &bm, Swip &swip) : PageGuard(bm, &bm.fix_page(swip, Exclusive)) {}
    ~PageGuard() { unfix(); }

    PageGuard& operator=(PageGuard&& other) {
        assert(&bm == &other.bm);
        if (this != &other) {
            unfix();
            page = other.page;
            if constexpr (Exclusive) dirty = other.dirty;
            other.reset();
        }
        return *this;
    }
    PageGuard& operator=(const PageGuard& other) = delete;

    void mark_dirty() { if constexpr(Exclusive) dirty = true; }
    void fix(Swip& swip) { unfix(); page = &bm.fix_page(swip, true); }
    void unfix() { if (page) bm.unfix_page(*page, dirty); reset(); }

    BufferFrame* operator->() { assert(page); return page; }
    BufferFrame* bf_ptr() const { return page; }
};
}
#endif