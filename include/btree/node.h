#ifndef INCLUDE_NODE_H
#define INCLUDE_NODE_H

#include <cstddef>
#include <cstdint>
#include <btree/iterators.h>

#define CACHELINE 64

namespace guidedresearch {

enum class NodeLayout {
    SORTED,
    EYTZINGER,
    EYTZINGER_SIMD
};

/// @brief Picks iterator based on node layout
/// @tparam layout 
/// @return Iterator Type
template <NodeLayout layout, typename KeyT>
consteval auto IteratorPicker() {
    if constexpr (layout == NodeLayout::EYTZINGER) {
        return std::type_identity<EytzingerIterator>{};
    } else if constexpr (layout == NodeLayout::EYTZINGER_SIMD) {
        return std::type_identity<BlockIterator<CACHELINE/sizeof(KeyT)>>{};
    } else {
        return std::type_identity<OrderedIterator>{};
    }
}

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
};
}
#endif