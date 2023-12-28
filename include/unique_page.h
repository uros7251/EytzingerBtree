#ifndef INCLUDE_UNIQUE_PAGE
#define INCLUDE_UNIQUE_PAGE

#include "buffer_manager.h"

namespace guidedresearch {

class UniquePage {
    BufferManager& bm;
    BufferFrame* page;
    bool dirty;

    void reset() { page = nullptr; dirty = false; }
public:
    UniquePage(BufferManager& bm) : bm(bm), page(nullptr) {}
    UniquePage(UniquePage&& other) : bm(other.bm), page(other.page) {
        other.reset();
    }
    UniquePage(const UniquePage& other) = delete;
    UniquePage(BufferManager& bm, Swip& swip) : bm(bm), page(&bm.fix_page(swip, true)) {}
    ~UniquePage() { unfix(); }

    UniquePage& operator=(UniquePage&& other) {
        assert(&bm == &other.bm);
        if (this != &other) {
            unfix();
            page = other.page;
            dirty = other.dirty;
            other.reset();
        }
        return *this;
    }
    UniquePage& operator=(const UniquePage& other) = delete;
    BufferFrame* operator->() { assert(page); return page; }
    BufferFrame* bf_ptr() const { return page; }

    void mark_dirty() { dirty = true; }
    void fix(Swip& swip) { unfix(); page = &bm.fix_page(swip, true); }
    void unfix() { if (page) bm.unfix_page(*page, dirty); reset(); }
};

} // namespace guidedresearch
#endif