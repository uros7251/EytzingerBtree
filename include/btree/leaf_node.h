#ifndef INCLUDE_LEAF_NODE_H_
#define INCLUDE_LEAF_NODE_H_

#include "btree/node.h"

namespace guidedresearch {

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize, bool FastInsertion>
struct LeafNode;

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct LeafNode<KeyT, ValueT, ComparatorT, PageSize, false> : public guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize> {

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
        for (auto i=index; i<Node::count-1u; ++i) {
            keys[i] = keys[i+1];
            values[i] = values[i+1]; 
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
        for (auto i=0, pos=(Node::count+1)/2; i<Node::count/2; ++i) {
            new_node->keys[i] = keys[pos+i];
            new_node->values[i] = values[pos+i];
        }
        new_node->count = Node::count/2; // floor(count/2)
        Node::count -= Node::count/2; // ceil(count/2)
        return keys[Node::count-1]; // return (new) rightmost key of old_node
    }

    void merge(LeafNode &other) {
        assert(Node::count+other.count <= max_keys());
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

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct LeafNode<KeyT, ValueT, ComparatorT, PageSize, true> : public guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize> {

    using Node = guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize>;
    /// The capacity of a node.
    static constexpr uint32_t kCapacity = (PageSize-sizeof(Node))/(sizeof(KeyT)+sizeof(ValueT)+sizeof(uint16_t));
    
    inline static constexpr uint32_t max_keys() { return kCapacity; }
    inline static constexpr uint32_t max_values() { return kCapacity; }

    /// The keys.
    KeyT keys[LeafNode::max_keys()]; // adjust this
    /// The values.
    ValueT values[LeafNode::max_values()]; // adjust this
    /// The indices.
    uint16_t indices[LeafNode::max_keys()];
    /// Constructor.
    LeafNode() : Node(0, 0) { }
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

    /// Get the index of the first key that is not less than than a provided key.
    /// @param[in] key          The key that should be inserted.
    /// @param[in] indirect     If false (default), return the index of the key in the keys array, otherwise return the index i such that keys[indices[i]] == key.
    /// @return                 The index of the first key that is not less than the provided key and boolean indicating if the key is equal.
    std::pair<uint32_t, bool> lower_bound(const KeyT &key, bool indirect=false) {
        if (Node::is_leaf() && Node::count==0) return {0u, false}; // no keys
        
        ComparatorT less;
        uint32_t i=0, n = Node::count;
        // branchless binary search
        while (n>1) {
            auto half = n/2;
            n -= half; // ceil(n/2)
            __builtin_prefetch(&indices[i+n/2-1]); // prefetch left
            __builtin_prefetch(&indices[i+half+n/2-1]); // prefetch right
            i += less(keys[indices[i+half-1]], key) * half; // hopefully compiler translates this to cmov
        }

        return i==Node::count-1u && less(keys[indices[i]], key) ? // check if last key is less than key
            std::make_pair(i+1, false) : // return index one after last key
            std::make_pair(indirect ? i : indices[i], keys[indices[i]] == key); 
    }

    /// Insert a key.
    /// @param[in] key          The key that should be inserted.
    /// @param[in] value        The value that should be inserted.
    void insert(const KeyT &key, const ValueT &value) {
        if (Node::count == kCapacity) throw std::runtime_error("Not enough space!");
        auto [index, found] = lower_bound(key, true);
        if (found) {
            values[indices[index]] = value;
            return;
        }
        keys[Node::count] = key;
        values[Node::count] = value;
        if (index < Node::count) { // if key should be inserted in the end, no need to move
            for (auto i=Node::count; i>index; --i) {
                indices[i] = indices[i-1];
            }
        }
        indices[index] = Node::count;
        ++Node::count;
    }

    /// Erase a key.
    bool erase(const KeyT &key) {
        // try to find key
        auto [index, found] = lower_bound(key, true);
        if (!found) return false; // key not found
        // find i such that indices[i] == Node::count-1
        auto [index_last, found_last] = lower_bound(keys[Node::count-1], true);
        assert(found_last);
        // relink indices 
        keys[indices[index]] = keys[Node::count-1];
        values[indices[index]] = values[Node::count-1];
        indices[index_last] = indices[index];
        // erase last key
        for (auto i=index; i<Node::count-1u; ++i) {
            indices[i] = indices[i+1];
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
        for (auto i=0, pos=(Node::count+1)/2; i<Node::count/2; ++i) {
            new_node->indices[i] =  i;
            new_node->keys[i] = keys[indices[pos+i]];
            new_node->values[i] = values[indices[pos+i]];
        }
        new_node->count = Node::count/2; // floor(count/2)
        // rearrange first half of keys and values in original node (using space in new node) 
        for (auto i=0, pos=Node::count/2; i<(Node::count+1)/2; ++i) {
            new_node->keys[i+pos] = keys[indices[i]];
            new_node->values[i+pos] = values[indices[i]];
        }
        for (auto i=0, pos=Node::count/2; i<(Node::count+1)/2; ++i) {
            indices[i] = i;
            keys[i] = new_node->keys[i+pos];
            values[i] = new_node->values[i+pos];
        }
        Node::count -= Node::count/2; // ceil(count/2)
        // return (new) rightmost key of old_node
        return keys[Node::count-1]; // note that here indices[Node::count-1] == Node::count-1
    }

    void merge(LeafNode &other) {
        assert(Node::count+other.count <= max_keys());
        for (uint16_t i=0; i<other.count; ++i) {
            indices[i+Node::count] = i+Node::count;
            keys[i+Node::count] = other.keys[other.indices[i]];
            values[i+Node::count] = other.values[other.indices[i]];
        }
        Node::count += other.count;
    }

    /// Returns the keys.
    std::vector<KeyT> get_key_vector() {
        std::vector<KeyT> sorted(Node::count);
        for (uint16_t i=0; i<Node::count; ++i) {
            sorted[i] = keys[indices[i]];
        }
        return sorted;
    }

    /// Returns the values.
    std::vector<ValueT> get_value_vector() {
        std::vector<ValueT> sorted(Node::count);
        for (uint16_t i=0; i<Node::count; ++i) {
            sorted[i] = values[indices[i]];
        }
        return sorted;
    }

    static void rebalance(LeafNode &left, LeafNode &right, KeyT &separator) {
        if (left.count > right.count) {
            // shift right
            uint16_t to_shift = (left.count-right.count)/2; // shift half the difference
            // make space for keys and values from the left node
            for (uint16_t i=right.count; i>0;) {
                --i;
                right.indices[i+to_shift] = right.indices[i]+to_shift; // since we are moving everything right, we need to shift the indices as well
                right.keys[i+to_shift] = right.keys[i];
                right.values[i+to_shift] = right.values[i];
            }
            // move keys and values from left node
            for (uint16_t i=0, pos=left.count-to_shift; i<to_shift; ++i) {
                right.indices[i] = i;
                right.keys[i] = left.keys[left.indices[pos+i]];
                right.values[i] = left.values[left.indices[pos+i]];
            }
            // rearrange keys and values in left node
            for (int i=0, j=left.count-to_shift, n=left.count-to_shift; i<n; ++i) {
                if (left.indices[i] >= n) {
                    while (left.indices[j] >= n) ++j;
                    assert(j < left.count);
                    left.keys[left.indices[j]] = left.keys[left.indices[i]];
                    left.values[left.indices[j]] = left.values[left.indices[i]];
                    left.indices[i] = left.indices[j];
                    ++j;
                }
            }
            separator = left.keys[left.indices[left.count-to_shift-1]];
            left.count -= to_shift;
            right.count += to_shift;
        }
        else {
            // shift left
            uint16_t to_shift = (right.count-left.count)/2; // shift half the difference
            // move keys and values from right node
            for (uint16_t i=0, pos=left.count; i<to_shift; ++i) {
                left.indices[i+pos] = i+pos;
                left.keys[i+pos] = right.keys[right.indices[i]];
                left.values[i+pos] = right.values[right.indices[i]];
            }
            // rearrange keys and values in right node
            for (int i=to_shift, j=0, n=right.count-to_shift; i<right.count; ++i) {
                if (right.indices[i] >= n) {
                    while (right.indices[j] >= n) ++j;
                    assert(j < to_shift);
                    right.keys[right.indices[j]] = right.keys[right.indices[i]];
                    right.values[right.indices[j]] = right.values[right.indices[i]];
                    right.indices[i] = right.indices[j];
                    ++j;
                }
            }
            // move keys and values to the beginning of the right node
            for (int i=0; i<right.count-to_shift; ++i) {
                right.indices[i] = right.indices[i+to_shift];
            }
            separator = left.keys[left.indices[left.count+to_shift-1]];
            left.count += to_shift;
            right.count -= to_shift;
        }
    }

private:
    void copy(const LeafNode &other) {
        Node::operator=(other);
        for (uint32_t i=0; i<Node::count; ++i) {
            indices[i] = other.indices[i];
            keys[i] = other.keys[i];
            values[i] = other.values[i];
        }
    }
};
}

#endif