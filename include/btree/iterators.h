#ifndef INCLUDE_ITERATORS_H
#define INCLUDE_ITERATORS_H

#include <cstdint>
#include <bit>
#include <cassert>

namespace guidedresearch {

/// @brief Iterator for in-order traversal of the keys stored as static btree
/// @tparam Block size in units of keys
template<uint16_t B>
struct BlockIterator {
    // /// @brief The current block index
    // uint16_t i;
    // /// @brief The current index within the block
    // uint16_t j;
    /// @brief The current index
    uint16_t k;
    /// @brief The node count
    uint16_t n;

    BlockIterator() = delete;
    BlockIterator(uint32_t i, uint16_t size) : k(i), n(size) {}

    bool operator==(const BlockIterator &other) const { return k == other.k; }

    uint32_t operator*() { return k; }

    /// @brief Calculates the index of the beginning of child block of entry k
    /// @return Returns the index of the beginning of child block of entry k
    uint16_t child() { 
        // i*(B+1)+j = i*B+j+i = k+i
        return B*(k+k/B);
    } 
    /// @brief Calculates the (block, offset) pair such that its child block is i
    /// @return Returns the entry (block, offset) such that its child block is i
    std::pair<uint16_t, uint16_t> parent(uint16_t i) { return {i/(B+1), i%(B+1)}; }

    /// @brief Advance the iterator to the next key in in-order traversal
    BlockIterator& operator++() {
        if (B+child() < n) {
            // go to right child
            k = B+child();
            // then left as much as possible
            while (child() < n) k=child();
            return *this;
        }
        uint16_t j = k % B + 1, i = k / B;
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
    BlockIterator& operator--() {
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
    static BlockIterator begin(uint16_t size) { assert(size); return ++BlockIterator(0u, size); }
    
    /// @brief Create an iterator pointing to the last key in reverse in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the last key in in-order traversal
    static BlockIterator rbegin(uint16_t size) {
        if (B >= size) return BlockIterator(size-1, size);
        uint16_t i=0;
        while(true) {
            if ((i+1)*B >= size) return BlockIterator(size-1, size);
            else if (((i+1)*B+i)*B >= size) return BlockIterator((i+1)*B-1, size);
            i = (i+1)*B+i; // = i*(B+1)+B
        }
    }

    static BlockIterator end(uint16_t size)  { return BlockIterator(0u, size); }
    static BlockIterator rend(uint16_t size)  { return BlockIterator(0u, size); }
};

/// @brief Iterator for in-order traversal of the keys stored as static binary tree
struct EytzingerIterator {
    /// @brief The current index
    uint32_t k;
    /// @brief The node count
    uint16_t n;

    EytzingerIterator() = delete;
    EytzingerIterator(uint32_t k, uint16_t size) :k(k), n(size) {}

    bool operator==(const EytzingerIterator &other) const = default;

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
    static EytzingerIterator begin(uint16_t size) { assert(size); return EytzingerIterator((1u << (32 - std::countl_zero(size-1u))) >> 1, size); } // greatest power of 2 less than n
    
    /// @brief Create an iterator pointing to the last key in reverse in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the last key in in-order traversal
    static EytzingerIterator rbegin(uint16_t size) { return EytzingerIterator(((1u << (16 - std::countl_zero(size))) >> 1) - 1u, size); } // (greatest power of 2 not bigger than n) - 1

    static EytzingerIterator end(uint16_t size)  { return EytzingerIterator(0u, size); }
    static EytzingerIterator rend(uint16_t size)  { return EytzingerIterator(0u, size); }
};

/// @brief Iterator for in-order traversal of the keys stored as sorted array
struct OrderedIterator {
    /// @brief The current index
    uint32_t k;
    /// @brief The node count
    uint16_t n;

    OrderedIterator() = delete;
    OrderedIterator(uint32_t k, uint16_t size) :k(k), n(size) {}

    bool operator==(const OrderedIterator &other) const = default;

    uint32_t operator*() { return k; }

    /// @brief Advance the iterator to the next key in in-order traversal
    OrderedIterator& operator++() {
        ++k;
        if (k==n) k=0;
        return *this;
    }

    /// @brief Advance the iterator to the previous key in in-order traversal
    OrderedIterator& operator--() {
        if (!k) k = n;
        --k;
        return *this;
    }

    /// @brief Create an iterator pointing to the first key in in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the first key in in-order traversal
    static OrderedIterator begin(uint16_t size) { assert(size); return OrderedIterator(1, size); } // greatest power of 2 less than n
    
    /// @brief Create an iterator pointing to the last key in reverse in-order traversal
    /// @param size The node count
    /// @return Iterator pointing to the last key in in-order traversal
    static OrderedIterator rbegin(uint16_t size) { return OrderedIterator(size-1, size); } // (greatest power of 2 not bigger than n) - 1

    static OrderedIterator end(uint16_t size)  { return OrderedIterator(0u, size); }
    static OrderedIterator rend(uint16_t size)  { return OrderedIterator(0u, size); }
};

} // namespace guidedresearch

#endif