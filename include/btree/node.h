#ifndef INCLUDE_NODE_H
#define INCLUDE_NODE_H

#include <cstddef>
#include <cstdint>

namespace guidedresearch {

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct Node {
    /// The level in the tree.
    uint16_t level;
    /// The number of children.
    uint16_t count;

    // Constructor
    Node(uint16_t level, uint16_t count)
        : level(level), count(count) {}

    /// Is the node a leaf node?
    bool is_leaf() const { return level == 0; }

    /// Offset of member 'keys'
    static constexpr size_t keys_offset = ((sizeof(Node)+sizeof(KeyT)-1)/sizeof(KeyT))*sizeof(KeyT); // take into account alignment requirements of KeyT

    /// Get the index of the first key that is not less than than a provided key.
    /// @param[in] key          The key that should be inserted.
    /// @return                 The index of the first key that is not less than the provided key and boolean indicating if the key is equal.
    std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
        const KeyT *keys = reinterpret_cast<KeyT*>(reinterpret_cast<size_t>(this)+keys_offset);
        if ((!is_leaf() && Node::count <= 1) || (is_leaf() && Node::count==0)) return {0u, false}; // no keys
        
        ComparatorT less;
        uint32_t i=0, n = Node::count-(!is_leaf());
        // branchless binary search
        while (n>1) {
            auto half = n/2;
            n -= half; // ceil(n/2)
            // __builtin_prefetch(&keys[i+n/2-1]); // prefetch left
            // __builtin_prefetch(&keys[i+half+n/2-1]); // prefetch right
            i += less(keys[i+half-1], key) * half; // hopefully compiler translates this to cmov
        }

        return ((!is_leaf() && i==Node::count-2u) || (is_leaf() && i==Node::count-1u)) && less(keys[i], key) ? // check if last key is less than key
            std::make_pair(i+1, false) : // return index one after last key
            std::make_pair(i, keys[i] == key); 
    }
};
}
#endif