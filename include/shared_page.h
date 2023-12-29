#ifndef INCLUDE_SHARED_PAGE
#define INCLUDE_SHARED_PAGE

#include "page_guard.h"

namespace guidedresearch {

class SharedPage : public PageGuard {
    /// @brief the number of objects holding a read reference to 'page'
    uint16_t *count = nullptr;
    
    void reset() { page = nullptr; count = nullptr; }
public:
    SharedPage(BufferManager& bm) : PageGuard(bm, nullptr), count(nullptr) {}
    SharedPage(BufferManager& bm, Swip& swip) : PageGuard(bm, &bm.fix_page(swip, false)), count(new uint16_t(1u)) {}
    SharedPage(const SharedPage& other) : PageGuard(other.bm, other.page), count(other.count) {
        ++(*count);
    }
    ~SharedPage() { unfix(); }

    SharedPage& operator=(const SharedPage& other) {
        if (this != &other) {
            assert(&bm == &other.bm);
            unfix();
            page = other.page;
            count = other.count;
            ++(*count);
        }
        return *this;
    }

    SharedPage& operator=(SharedPage&& other) {
        if (this != &other) {
            page = other.page;
            count = other.count;
            other.reset();
        }
        return *this;
    }

    void fix(Swip& swip) { 
        unfix(); 
        page = &bm.fix_page(swip, false); 
        count = new uint16_t(1u);
    }

    /// This method does not neccessarily unfix the page
    void unfix() {
        if (count && !(--(*count))) {
            assert(page);
            delete count;
            bm.unfix_page(*page, false);
        }
        reset();
    }
};
}
#endif