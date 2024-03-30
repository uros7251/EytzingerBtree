#ifndef INCLUDE_ITERATORS_H
#define INCLUDE_ITERATORS_H

#include <cstdint>
#include <bit>
#include <cassert>

namespace guidedresearch {

template<uint16_t Capacity, uint16_t B>
struct Precalc {
    uint16_t sizes[10];
    uint16_t array[Capacity+1];
    uint16_t index[Capacity+1];

    constexpr Precalc() :array{} {
        array[0] = 0;
        index[0] = 0;
        for (uint16_t k=1, i=0; k<Capacity; ++k) {
            i = next(i/B, i%B);
            array[k] = i;
            index[i] = k;
        }
        array[Capacity] = 0;
        index[Capacity] = Capacity;
        sizes[0]=Capacity;
        for (auto i=1; i<10; ++i) {
            sizes[i]=Capacity/(1<<i);
        }
        
    }

    private:
    static constexpr auto child_block(uint16_t i, uint16_t j) {
        return i*(B+1)+j;
    }
    static constexpr auto parent_elem_index(uint16_t block_index) {
        return block_index%(B+1);
    }
    static constexpr uint16_t parent_block_index(uint16_t block_index) {
        return block_index/(B+1);
    }
    constexpr auto next(uint16_t i, uint16_t j) {
        if (B*(1+child_block(i,j)) < Capacity) {
            // go to right child
            i=child_block(i,j)+1;
            j=0;
            while (B*child_block(i,j) < Capacity) i=child_block(i,j);
            return i*B+j;
        }
        ++j; // increment j
        while (j>=B || i*B+j>=Capacity) {
            // climb up until you're not rightmost child and you're within bounds
            std::tie(i,j) = std::pair{parent_block_index(i), parent_elem_index(i)};
        }
        return i*B+j;
    }
};

template<uint16_t Capacity,uint16_t B>
struct StaticBlockIterator {
    static constexpr Precalc<Capacity, B> precalc{};
    uint16_t i;
    uint16_t n;

    StaticBlockIterator() = delete;
    StaticBlockIterator(uint16_t i, uint16_t n) :i(precalc.index[i]), n(n) {}

    bool operator==(const StaticBlockIterator &other) const = default;
    uint16_t operator*() { return precalc.array[i]; }

    StaticBlockIterator& operator++() {
        assert(i<Capacity);
        do ++i; while (precalc.array[i]>=n);
        if (i==Capacity) i=0;
        return *this;
    }

    StaticBlockIterator& operator--() {
        if (i==0) i=Capacity;
        do --i; while (precalc.array[i]>=n);
        return *this;
    }

    static StaticBlockIterator begin(uint16_t size) { assert(size); return ++StaticBlockIterator(0u, size); }
    static StaticBlockIterator rbegin(uint16_t size) { return --StaticBlockIterator(Capacity, size); }
    static StaticBlockIterator end(uint16_t size)  { return StaticBlockIterator(0u, size); }
    static StaticBlockIterator rend(uint16_t size)  { return StaticBlockIterator(0u, size); }
};

/// @brief Iterator for in-order traversal of the keys stored as static btree
/// @tparam Block size in units of keys
template<uint32_t B>
struct EytzingerIterator {
    // /// @brief The current block index
    // uint32_t i;
    // /// @brief The current index within the block
    // uint32_t j;
    /// @brief The current index
    uint32_t k;
    /// @brief The node count
    uint32_t n;

    EytzingerIterator() = delete;
    EytzingerIterator(uint32_t i, uint32_t size) : k(i), n(size) {}
    EytzingerIterator(const EytzingerIterator&) = default;

    bool operator==(const EytzingerIterator &other) const { assert(n==other.n); return k == other.k; }
    EytzingerIterator& operator=(const EytzingerIterator &other) = default;

    uint32_t operator*() { return k; }

    /// @brief Calculates the index of the beginning of child block of entry k
    /// @return Returns the index of the beginning of child block of entry k
    uint32_t child() { 
        // i*(B+1)+j = i*B+j+i = k+i
        return B*(k+k/B);
    } 
    /// @brief Calculates the (block, offset) pair such that its child block is i
    /// @return Returns the entry (block, offset) such that its child block is i
    std::pair<uint32_t, uint32_t> parent(uint32_t i) { return {i/(B+1), i%(B+1)}; }

    /// @brief Advance the iterator to the next key in in-order traversal
    EytzingerIterator& operator++() {
        if (B+child() < n) {
            // go to right child
            k = B+child();
            // then left as much as possible
            while (child() < n) k=child();
            return *this;
        }
        uint32_t j = k % B + 1, i = k / B;
        if (k+1u<n && j<B) {
            ++k;
            return *this;
        }
        do {
            std::tie(i,j) = parent(i);
        } while (j>=B);
        k = i*B+j;
        return *this;
    }

    /// @brief Advance the iterator to the next key in in-order traversal
    EytzingerIterator& operator--() {
        if (child() < n) {
            k=child()+B-1; // go to left child's rightmost entry
            // go right until either rightmost entry is out of bounds or 
            while (k < n && child()+B < n) k=child()+2*B-1;
            // at this point either:
            // 1. k is out of bounds but k-B+1 is in bounds or
            // 2. k is in bounds
            k = k >= n ? n-1u : k;
            return *this;
        }
        if (k % B > 0) {
            // go left
            --k;
            return *this;
        }
        auto j = k%B, i = k/B;
        while (j==0 && i!=0) {
            // climb up until you're not leftmost child or you're at root block
            std::tie(i,j) = parent(i);
        }
        if (j) --j;
        k = i*B+j;
        return *this;
    }

    /// @brief Create an iterator pointing to the first key in in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the first key in in-order traversal
    static EytzingerIterator begin(uint32_t size) { assert(size); return ++EytzingerIterator(0u, size); }
    
    /// @brief Create an iterator pointing to the last key in reverse in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the last key in in-order traversal
    static EytzingerIterator rbegin(uint32_t size) {
        if (B >= size) return EytzingerIterator(size-1, size);
        uint32_t i=0;
        while(true) {
            if ((i+1)*B >= size) return EytzingerIterator(size-1, size);
            else if (((i+1)*B+i)*B >= size) return EytzingerIterator((i+1)*B-1, size);
            i = (i+1)*B+i; // = i*(B+1)+B
        }
    }

    static EytzingerIterator end(uint32_t size)  { return EytzingerIterator(0u, size); }
    static EytzingerIterator rend(uint32_t size)  { return EytzingerIterator(0u, size); }
};

/// @brief Iterator for in-order traversal of the keys stored as static binary tree
template<>
struct EytzingerIterator<1u> {
    /// @brief The current index
    uint32_t k;
    /// @brief The node count
    uint32_t n;

    EytzingerIterator() = delete;
    EytzingerIterator(uint32_t k, uint32_t size) :k(k), n(size) {}
    EytzingerIterator(const EytzingerIterator&) = default;

    bool operator==(const EytzingerIterator &other) const { assert(n==other.n); return k == other.k; }
    EytzingerIterator& operator=(const EytzingerIterator &other) = default;

    uint32_t operator*() { return k; }

    /// @brief Advance the iterator to the next key in in-order traversal
    EytzingerIterator& operator++() {
        if (2*k+1<n) {
            // go right
            k = 2*k+1;
            // then left as much as possible
            while (2*k<n) k = 2*k;
        }
        else {
            k = k >> (std::countr_one(k)+1); // climb up until you're the left child and then climb up one more time
        }
        return *this;
    }

    /// @brief Advance the iterator to the previous key in in-order traversal
    EytzingerIterator& operator--() {
        if (2*k<n) {
            // go left
            k = 2*k;
            // then right as much as possible
            while (2*k+1<n) k = 2*k+1;
        }
        else {
            k = k >> (std::countr_zero(k)+1); // climb up until you're the right child and then climb up one more time
        }
        return *this;
    }

    /// @brief Create an iterator pointing to the first key in in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the first key in in-order traversal
    static EytzingerIterator begin(uint32_t size) { assert(size); return EytzingerIterator((1u << (32 - std::countl_zero(size-1u))) >> 1, size); } // greatest power of 2 less than n
    
    /// @brief Create an iterator pointing to the last key in reverse in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the last key in in-order traversal
    static EytzingerIterator rbegin(uint32_t size) { return EytzingerIterator(((1u << (32 - std::countl_zero(size))) >> 1) - 1u, size); } // (greatest power of 2 not bigger than n) - 1

    static EytzingerIterator end(uint32_t size)  { return EytzingerIterator(0u, size); }
    static EytzingerIterator rend(uint32_t size)  { return EytzingerIterator(0u, size); }
};

/// @brief Iterator for in-order traversal of the keys stored as sorted array
struct OrderedIterator {
    /// @brief The current index
    uint32_t k;
    /// @brief The node count
    uint32_t n;

    OrderedIterator() = delete;
    OrderedIterator(uint32_t k, uint32_t size) :k(k), n(size) {}
    OrderedIterator(const OrderedIterator&) = default;

    bool operator==(const OrderedIterator &other) const { assert(n==other.n); return k == other.k; }
    OrderedIterator& operator=(const OrderedIterator &other) = default;

    uint32_t operator*() { return k; }

    /// @brief Advance the iterator to the next key in in-order traversal
    OrderedIterator& operator++() {
        k = k+1==n ? 0 : k+1;
        return *this;
    }

    /// @brief Advance the iterator to the previous key in in-order traversal
    OrderedIterator& operator--() {
        k = k ? k-1 : n-1;
        return *this;
    }

    /// @brief Create an iterator pointing to the first key in in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the first key in in-order traversal
    static OrderedIterator begin(uint32_t size) { assert(size); return OrderedIterator(std::min(size-1u, 1u), size); }
    
    /// @brief Create an iterator pointing to the last key in reverse in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the last key in in-order traversal
    static OrderedIterator rbegin(uint32_t size) { assert(size); return OrderedIterator(size-1u, size); }

    static OrderedIterator end(uint32_t size)  { return OrderedIterator(0u, size); }
    static OrderedIterator rend(uint32_t size)  { return OrderedIterator(0u, size); }
};

} // namespace guidedresearch

#endif