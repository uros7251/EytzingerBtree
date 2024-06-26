#ifndef INCLUDE_NODE_H
#define INCLUDE_NODE_H

#include <cstddef>
#include <cstdint>
#include <btree/iterators.h>

#define kInf std::numeric_limits<KeyT>::max()
#define kNegInf std::numeric_limits<KeyT>::min()

#define CACHELINE 64

namespace guidedresearch {

enum class NodeLayout {
    SORTED,
    EYTZINGER,
    EYTZINGER_SIMD,
    ORDERED     // also sorted, but using iterators for updates
};

/// @brief Picks iterator based on node layout
/// @tparam layout 
/// @return Iterator Type
template <NodeLayout layout, typename KeyT, uint16_t Cap = 0>
consteval auto IteratorPicker() {
    if constexpr (layout == NodeLayout::EYTZINGER) {
        return std::type_identity<EytzingerIterator<1u>>{};
    } else if constexpr (layout == NodeLayout::EYTZINGER_SIMD) {
        // return std::type_identity<StaticBlockIterator<Cap, CACHELINE/sizeof(KeyT)>>{};
        return std::type_identity<EytzingerIterator<CACHELINE/sizeof(KeyT)>>{};
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