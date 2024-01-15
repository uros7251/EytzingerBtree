#ifndef INCLUDE_UNIQUE_PAGE
#define INCLUDE_UNIQUE_PAGE

#include "buffer_manager/page_guard.h"

namespace guidedresearch {

class UniquePage : public PageGuard {
    bool dirty;

    void reset() { page = nullptr; dirty = false; }
public:
    UniquePage(BufferManager& bm) : PageGuard(bm, nullptr) {}
    UniquePage(UniquePage&& other) : PageGuard(other.bm, other.page) {
        other.reset();
    }
    UniquePage(const UniquePage& other) = delete;
    UniquePage(BufferManager& bm, Swip& swip) : PageGuard(bm, &bm.fix_page(swip, true)) {}
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

    void mark_dirty() { dirty = true; }
    void fix(Swip& swip) { unfix(); page = &bm.fix_page(swip, true); }
    void unfix() { if (page) bm.unfix_page(*page, dirty); reset(); }
};

} // namespace guidedresearch
#endif