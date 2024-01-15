#ifndef INCLUDE_LEAF_NODE_H_
#define INCLUDE_LEAF_NODE_H_

#include "swip.h"
#include "node.h"

namespace guidedresearch {

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct LeafNode: public guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize> {

    using Node = guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize>;
    /// The capacity of a node.
    static constexpr uint32_t kCapacity = (PageSize-sizeof(Node))/(sizeof(KeyT)+sizeof(ValueT));
    
    inline static constexpr uint32_t max_keys() { return kCapacity; }
    inline static constexpr uint32_t max_values() { return kCapacity; }

    /// The keys.
    KeyT keys[LeafNode::max_keys()]; // adjust this
    /// The values.
    ValueT values[LeafNode::max_values()]; // adjust this
    /// Constructor.
    LeafNode() : Node(0, 0) { assert(reinterpret_cast<size_t>(&keys[0])==reinterpret_cast<size_t>(this)+Node::keys_offset); }
    LeafNode(const LeafNode &other) : Node(other) {
        copy();
    }
    
    /// Assignment operator
    LeafNode& operator=(const LeafNode &other) {
        copy(other);
        return *this;
    }

    /// @brief check if node is less than half full
    /// @return 
    bool merge_needed() const { return Node::count < max_values()/2; } 

    /// Insert a key.
    /// @param[in] key          The key that should be inserted.
    /// @param[in] value        The value that should be inserted.
    void insert(const KeyT &key, const ValueT &value) {
        if (Node::count == kCapacity) throw std::runtime_error("Not enough space!");
        auto [index, found] = Node::lower_bound(key);
        if (!found && index < Node::count) { // if key should be inserted in the end, no need to move
            for (auto i=Node::count; i>index; --i) {
                keys[i] = keys[i-1];
                values[i] = values[i-1];
            }
        }
        keys[index] = key;
        values[index] = value;
        if (!found) ++Node::count;
    }

    /// Erase a key.
    bool erase(const KeyT &key) {
        // try to find key
        auto [index, found] = Node::lower_bound(key);
        if (!found) return false; // key not found
        if (index < Node::count-1u) { // if key is last key, no need to move
            for (auto i=index; i<Node::count-1u; ++i) {
                keys[i] = keys[i+1];
                values[i] = values[i+1]; 
            }
        }
        --Node::count;
        return true;
    }

    /// Split the node.
    /// @param[in] buffer       The buffer for the new page.
    /// @return                 The separator key.
    KeyT split(char* buffer) {
        // for example, keys = [1,2,3,4,5,6,7,8,9,10], count = 10
        // old_node->keys = [1,2,3,4,5], old_node->count = 5
        // new_node->keys = [6,7,8,9,10], new_node->count = 5
        LeafNode* new_node = new (buffer) LeafNode();
        // old_node retains ceil(count/2) keys
        // new_node gets floor(count/2) keys
        // move keys and values starting from index ceil(count/2), in total floor(count/2)
        for (auto i=0u, pos=(Node::count+1)/2; i<Node::count/2; ++i) {
            new_node->keys[i] = keys[pos+i];
            new_node->values[i] = values[pos+i];
        }
        new_node->count = Node::count/2; // floor(count/2)
        Node::count -= Node::count/2; // ceil(count/2)
        return keys[Node::count-1]; // return (new) rightmost key of old_node
    }

    void merge(LeafNode &other) {
        for (uint16_t i=0; i<other.count; ++i) {
            keys[i+Node::count] = other.keys[i];
            values[i+Node::count] = other.values[i];
        }
        Node::count += other.count;
    }

    /// Returns the keys.
    std::vector<KeyT> get_key_vector() {
        return std::vector<KeyT>(keys, keys+Node::count);
    }

    /// Returns the values.
    std::vector<ValueT> get_value_vector() {
        return std::vector<ValueT>(values, values+Node::count);
    }

    static void rebalance(LeafNode &left, LeafNode &right, KeyT &separator) {
        if (left.count > right.count) {
            uint16_t to_shift = (left.count-right.count)/2; // shift half the difference
            // make space for keys and values from the left node
            for (uint16_t i=right.count; i>0;) {
                --i;
                right.keys[i+to_shift] = right.keys[i];
                right.values[i+to_shift] = right.values[i];
            }
            for (uint16_t i=0, pos=left.count-to_shift; i<to_shift; ++i) {
                right.keys[i] = left.keys[pos+i];
                right.values[i] = left.values[pos+i];
            }
            separator = left.keys[left.count-to_shift-1];
            left.count -= to_shift;
            right.count += to_shift;
        }
        else {
            uint16_t to_shift = (right.count-left.count)/2; // shift half the difference
            
            for (uint16_t i=0; i<to_shift; ++i) {
                left.keys[i+left.count] = right.keys[i];
                left.values[i+left.count] = right.values[i];
            }
            for (uint16_t i=0; i<right.count-to_shift; ++i) {
                right.keys[i] = right.keys[i+to_shift];
                right.values[i] = right.values[i+to_shift];
            }
            separator = left.keys[left.count+to_shift-1];
            left.count += to_shift;
            right.count -= to_shift;
        }
    }

private:
    void copy(const LeafNode &other) {
        Node::operator=(other);
        for (uint32_t i=0; i<Node::count; ++i) {
            keys[i] = other.keys[i];
            values[i] = other.values[i];
        }
    }
};
}

#endif