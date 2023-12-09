#ifndef INCLUDE_BTREE_H
#define INCLUDE_BTREE_H

#include <type_traits>
#include <cstring>

#include "buffer_manager.h"
#include "segment.h"
#include "swip.h"

namespace guidedresearch {

//------------------------------------------------------------------------------------------------
// The algorithmic part of this BTree implementation is heavily inspired
// by prof. Neumann's article https://databasearchitects.blogspot.com/2022/06/btreeoperations.html
//------------------------------------------------------------------------------------------------

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct BTree : public Segment {
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

    struct InnerNode: public Node {
        /// Remark: it might be better to store the rightmost key to make nodes more uniform and algorithms more concise

        /// The capacity of a node.
        static constexpr uint32_t kCapacity = (PageSize-sizeof(Node))/(sizeof(KeyT)+sizeof(Swip));

        inline static constexpr uint32_t max_keys() { return kCapacity; }
        inline static constexpr uint32_t max_children() { return kCapacity+1; }

        /// The keys.
        KeyT keys[InnerNode::max_keys()]; 
        /// The values.
        Swip children[InnerNode::max_children()]; 

        /// Constructor.
        InnerNode() : Node(0, 0) {}
        InnerNode(uint16_t level) : Node(level, 0) {}

        /// this method can be virtual, but we avoid that on purpose
        /// @brief check if node is less than half full
        /// @return 
        bool merge_needed() const { return Node::count < max_children()/2; } 

        /// Get the index of the first key that is not less than than a provided key.
        /// @param[in] key          The key that should be inserted.
        /// @return                 The index of the first key that is not less than the provided key and boolean indicating if the key is equal.
        std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
            
            if (Node::count <= 1) return {0u, false}; // no keys
            
            ComparatorT less;
            uint32_t i=0, n = Node::count-1;
            // branchless binary search
            while (n>1) {
                auto half = n/2;
                i += less(keys[i+half-1], key) * half; // hopefully compiler translates this to cmov
                n -= half; // ceil(n/2)
            }

            return i==Node::count-2u && less(keys[i], key) ? // check if last key is less than key
                std::make_pair(i+1, false) : // return index one after last key
                std::make_pair(i, keys[i] == key); 
        }

        /// Insert a key.
        /// @param[in] key          The key that should be inserted.
        /// @param[in] child        The child that should be inserted.
        void insert_split(const KeyT &key, Swip child) {
            // for example, insert_split(3, 40)
            // BEFORE: keys = [1, 4, 9], pages = [8, 3, 13, 7]
            // AFTER: keys = [1, 3, 4, 9], pages = [8, 3, 40, 13, 7]

            if (Node::count-1 == InnerNode::kCapacity) throw std::runtime_error("insert_split(): Not enough space!");
            auto [index, _] = lower_bound(key); // find index to insert key
            std::memmove(&keys[index+1], &keys[index], sizeof(KeyT)*(Node::count-1-index)); // move keys
            std::memmove(&children[index+2], &children[index+1], sizeof(uint64_t)*(Node::count-1-index)); // move children
            keys[index] = key; // insert new separator
            children[index+1] = child; // insert child
            Node::count++;
        }

        /// @brief Erases a key and assosiacted child
        /// @param index Index of the key and child to be erased
        void erase(uint16_t index) {
            if (index < Node::count-1u) {
                std::memmove(&keys[index], &keys[index+1], sizeof(KeyT)*(Node::count-(index+2)));
                std::memmove(&children[index], &children[index+1], sizeof(Swip)*(Node::count-(index+1)));
            }
            --Node::count;
        }

        /// Split the node.
        /// @param[in] buffer       The buffer for the new page.
        /// @return                 The separator key.
        KeyT split(char* buffer) {
            // for example, keys = [1,2,3,4,5,6,7,8,9,10], count = 11 (number of children, not keys)
            // old_node->keys = [1,2,3,4,5], old_node->count = 6
            // new_node->keys = [7,8,9,10], new_node->count = 5
            // separator = 6
            auto* new_node = new (buffer) InnerNode(Node::level);
            // old_node retains floor((count-1)/2) keys
            // new_node gets floor(count/2)-1 keys
            // move keys starting from index floor((count-1)/2)+1=floor((count+1)/2) (skip separator), in total ceil((count-1)/2)-1=floor(count/2)-1
            std::memcpy(&new_node->keys[0], &keys[(Node::count+1)/2], sizeof(KeyT)*(Node::count/2-1));
            // move children starting from index ceil(count/2), in total floor(count/2)
            std::memcpy(&new_node->children[0], &children[(Node::count+1)/2], sizeof(ValueT)*(Node::count/2));
            new_node->count = Node::count/2; // floor(count/2)
            Node::count -= Node::count/2; // ceil(count/2)
            return keys[Node::count-1]; // Node::count-1 = new number of keys in old_node
        }

        /// Returns the keys.
        std::vector<KeyT> get_key_vector() {
            return std::vector<KeyT>(&keys[0], &keys[0]+Node::count);
        }
    };

    struct LeafNode: public Node {
        /// The capacity of a node.
        static constexpr uint32_t kCapacity = (PageSize-sizeof(Node)-sizeof(uint64_t))/(sizeof(KeyT)+sizeof(ValueT));
        
        inline static constexpr uint32_t max_keys() { return kCapacity; }
        inline static constexpr uint32_t max_values() { return kCapacity; }

        /// The keys.
        KeyT keys[LeafNode::max_keys()]; // adjust this
        /// The values.
        ValueT values[LeafNode::max_values()]; // adjust this
        /// Constructor.
        LeafNode() : Node(0, 0) {}

        /// @brief check if node is less than half full
        /// @return 
        bool merge_needed() const { return Node::count < max_values()/2; } 

        /// Get the index of the first key that is not less than than a provided key.
        std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
            if (Node::count==0) return std::make_pair(0u, false); // no keys
            
            ComparatorT less;
            uint32_t i = 0, n = Node::count;
            // branchless binary search
            while (n>1) {
                auto half = n/2;
                i += less(keys[i+half-1], key) * half; // hopefully compiler translates this to cmov
                n -= half; // ceil(n/2)
            }

            return i==Node::count-1u && less(keys[i], key) ? // check if last key is less than key
                std::make_pair(i+1, false) : // return index one after last key
                std::make_pair(i, keys[i] == key);
        }

        /// Insert a key.
        /// @param[in] key          The key that should be inserted.
        /// @param[in] value        The value that should be inserted.
        void insert(const KeyT &key, const ValueT &value) {
            if (Node::count == kCapacity) throw std::runtime_error("Not enough space!");
            auto [index, found] = lower_bound(key);
            if (!found && index < Node::count) { // if key should be inserted in the end, no need to move
                std::memmove(&keys[index+1], &keys[index], sizeof(KeyT)*(Node::count-index));
                std::memmove(&values[index+1], &values[index], sizeof(ValueT)*(Node::count-index));
            }
            keys[index] = key;
            values[index] = value;
            if (!found) ++Node::count;
        }

        /// Erase a key.
        void erase(const KeyT &key) {
            // try to find key
            auto [index, found] = lower_bound(key);
            if (!found) return; // key not found
            erase(index, 0);
        }

        /// Split the node.
        /// @param[in] buffer       The buffer for the new page.
        /// @return                 The separator key.
        KeyT split(char* buffer) {
            // for example, keys = [1,2,3,4,5,6,7,8,9,10], count = 10
            // old_node->keys = [1,2,3,4,5], old_node->count = 5
            // new_node->keys = [6,7,8,9,10], new_node->count = 5
            auto* new_node = new (buffer) LeafNode();
            // old_node retains ceil(count/2) keys
            // new_node gets floor(count/2) keys
            // move keys starting from index ceil(count/2), in total floor(count/2)
            std::memcpy(&new_node->keys[0], &keys[(Node::count+1)/2], sizeof(KeyT)*(Node::count/2));
            // move values starting from index ceil(count/2), in total floor(count/2)
            std::memcpy(&new_node->values[0], &values[(Node::count+1)/2], sizeof(ValueT)*(Node::count/2));
            new_node->count = Node::count/2; // floor(count/2)
            Node::count -= Node::count/2; // ceil(count/2)
            return keys[Node::count-1]; // return (new) rightmost key of old_node
        }

        /// Returns the keys.
        std::vector<KeyT> get_key_vector() {
            return std::vector<KeyT>(keys, keys+Node::count);
        }

        /// Returns the values.
        std::vector<ValueT> get_value_vector() {
            return std::vector<ValueT>(values, values+Node::count);
        }

        private:
        inline void erase(uint16_t index, int) {
            if (index < Node::count-1u) { // if key is last key, no need to move
                std::memmove(&keys[index], &keys[index+1], sizeof(KeyT)*(Node::count-(index+1)));
                std::memmove(&values[index], &values[index+1], sizeof(ValueT)*(Node::count-(index+1)));
            }
            --Node::count;
        }
    };

    /// The root.
    Swip root;
    /// The next free page;
    uint64_t next_available_page;

    /// Constructor.
    BTree(uint16_t segment_id, BufferManager &buffer_manager)
        : Segment(segment_id, buffer_manager), root(BufferManager::get_page_id(segment_id, 0)), next_available_page(1) {
        
        auto &page = buffer_manager.fix_page(root, true);
        new (page.get_data()) LeafNode();
        buffer_manager.unfix_page(page, true);
    }

    /// Destructor.
    ~BTree() = default;

    /// Lookup an entry in the tree.
    /// @param[in] key      The key that should be searched.
    /// @return             Whether the key was in the tree.
    std::optional<ValueT> lookup(const KeyT &key) {
        // q: why we need lock coupling?
        // a: because it might happen that some other thread is running insert operation and it 
        // splits a node. So, for example, the key we are looking for ends up on a different node
        BufferFrame *page = &buffer_manager.fix_page(root, false), *parent_page = nullptr;
        auto *node = reinterpret_cast<Node*>(page->get_data());
        while (!node->is_leaf()) {
            // traverse inner nodes
            InnerNode *inner_node = reinterpret_cast<InnerNode*>(node);
            auto [index, match] = inner_node->lower_bound(key);
            Swip swip = inner_node->children[index];
            parent_page = page;
            // lock coupling
            page = &buffer_manager.fix_page(swip, false);
            buffer_manager.unfix_page(*parent_page, false); 
            node = reinterpret_cast<Node*>(page->get_data());
        }
        // leaf level
        auto* leaf_node = reinterpret_cast<LeafNode*>(node);
        auto [index, match] = leaf_node->lower_bound(key);
        std::optional<ValueT> result;
        if (match) result = {leaf_node->values[index]};
        buffer_manager.unfix_page(*page, false);
        return result;
    }

    template<typename T>
    requires std::is_base_of_v<Node, T>
    void handle_root(BufferFrame* *child_page, BufferFrame* *parent_page) {
        // create a new node and copy the contents of the child node to it
        // modify initial child node such that it's now an InnerNode with only one child - newly created node
        Swip swip = next_page_id();
        auto &new_page = buffer_manager.fix_page(swip, true);
        T *new_node = new (new_page.get_data()) T(); 
        *new_node = *reinterpret_cast<T*>((*child_page)->get_data()); // default copy assignment
        InnerNode *old_node = new ((*child_page)->get_data()) InnerNode(new_node->level+1);
        old_node->children[0]=swip;
        old_node->count=1;
        *parent_page = *child_page;
        *child_page = &new_page;
    }

    BufferFrame* handle_root(uint16_t level) {
        Swip swip = next_page_id();
        auto &page = buffer_manager.fix_page(swip, true);
        InnerNode *parent = new (page.get_data()) InnerNode(level+1); // placement new
        parent->count = 1;
        parent->children[0] = root;
        root = swip;
        return &page;
    }

    template<typename T>
    requires std::is_base_of_v<Node, T>
    BufferFrame* split_node(BufferFrame *child_page, BufferFrame *parent_page, KeyT target) {
        // split
        T *node = reinterpret_cast<T*>(child_page->get_data());
        InnerNode *parent_node = reinterpret_cast<InnerNode*>(parent_page->get_data());
        Swip swip = next_page_id();
        auto &page = buffer_manager.fix_page(swip, true);
        KeyT new_separator = node->split(page.get_data());
        // insert new separator into parent node
        parent_node->insert_split(new_separator, swip);
        buffer_manager.unfix_page(*parent_page, true);
        if (target > new_separator) {
            buffer_manager.unfix_page(*child_page, true);
            return &page;
        }
        else {
            buffer_manager.unfix_page(page, true);
            return child_page;
        }
    }

    /// Inserts a new entry into the tree.
    /// @param[in] key      The key that should be inserted.
    /// @param[in] value    The value that should be inserted.
    void insert(const KeyT &key, const ValueT &value) {
        BufferFrame *child_page = &buffer_manager.fix_page(root, true),
            *parent_page = nullptr;
        Node *child = reinterpret_cast<Node*>(child_page->get_data());
        if (child->is_leaf()) {
            LeafNode *child_as_leaf = reinterpret_cast<LeafNode*>(child);
            if (child_as_leaf->count == LeafNode::max_values()) {
                handle_root<LeafNode>(&child_page, &parent_page);
                child = reinterpret_cast<Node*>(child_page->get_data());
            }
        }
        else {
            InnerNode *child_as_inner = reinterpret_cast<InnerNode*>(child);
            if (child_as_inner->count == InnerNode::max_children()) {
                handle_root<InnerNode>(&child_page, &parent_page);
                child = reinterpret_cast<Node*>(child_page->get_data());
            }
        }
        while (!child->is_leaf()) {
            if (child->count == InnerNode::max_children()) {
                child_page = split_node<InnerNode>(child_page, parent_page, key);
            }
            else if (parent_page) {
                buffer_manager.unfix_page(*parent_page, false); // no writes to parent
            }
            InnerNode *child_as_inner = reinterpret_cast<InnerNode*>(child_page->get_data());
            auto [index, match] = child_as_inner->lower_bound(key);
            Swip &swip = child_as_inner->children[index]; // it's crucial that swip is a reference
            parent_page = child_page;
            child_page = &buffer_manager.fix_page(swip, true);
            child = reinterpret_cast<Node*>(child_page->get_data());
        }
        // at the bottom
        if (child->count == LeafNode::max_values()) {
            child_page = split_node<LeafNode>(child_page, parent_page, key);
        }
        else if (parent_page) buffer_manager.unfix_page(*parent_page, false);
        
        LeafNode *child_as_leaf = reinterpret_cast<LeafNode*>(child_page->get_data());
        child_as_leaf->insert(key, value);
        buffer_manager.unfix_page(*child_page, true);
    }

    /// Erase a key.
    void erase(const KeyT &key) {
        BufferFrame *child_page = &buffer_manager.fix_page(root, true),
            *parent_page = nullptr;
        Node *child = reinterpret_cast<Node*>(child_page->get_data());
        InnerNode *parent = nullptr;
        bool dirty = false;
        while (!child->is_leaf()) {
            parent_page = child_page;
            parent = static_cast<InnerNode*>(child);
            auto [index, match] = parent->lower_bound(key);
            Swip &swip = parent->children[index]; // it's crucial that swip is a reference because fix_page might change it
            child_page = &buffer_manager.fix_page(swip, true);
            child = reinterpret_cast<Node*>(child_page->get_data());
            if (child->is_leaf()) {
                LeafNode *child_as_leaf = static_cast<LeafNode*>(child);
                if (child_as_leaf->merge_needed()) {
                    // this call will unfix the pages that need to be unfixed
                    // note: if we merge into the parent, we should not unfix the parent page
                    child_page = merge(parent_page, child_page, index, key);
                    // we changed the child_page, which is a future parent
                    dirty = true;
                }
                else {
                    buffer_manager.unfix_page(*parent_page, dirty);
                    dirty = false;
                }
            }
            else {
                InnerNode *child_as_inner = static_cast<InnerNode*>(child);
                if (child_as_inner->merge_needed()) {
                    child_page = merge(parent_page, child_page, index, key);
                    dirty = true;
                }
                else {
                    buffer_manager.unfix_page(*parent_page, dirty);
                    dirty = false;
                }
            }
            child = reinterpret_cast<Node*>(child_page->get_data());
        }
        LeafNode *leaf = static_cast<LeafNode*>(child);
        leaf->erase(key);
        buffer_manager.unfix_page(*child_page, true);
    }

    protected:
    inline Swip next_page_id() {
        Swip swip (buffer_manager.get_page_id(segment_id, next_available_page++));
        return swip;
    }

    /// @brief allocates a new node and locks it exclusively
    /// @tparam node type
    /// @return pointer to allocated node
    template<typename T>
    requires std::is_base_of_v<Node, T>
    T* allocate_node() {
        Swip swip ((static_cast<uint64_t>(segment_id) << 48) ^ next_available_page++);
        auto &page = buffer_manager.fix_page(swip, true);
        return new (page.get_data()) T();
    }

    BufferFrame* merge(BufferFrame *parent_page, BufferFrame *child_page, uint16_t child_slot, KeyT target) {
            // hold tight, this is 200+ lines function

            // assumes child and parent pages are already fixed
            BufferFrame *neighboor_page;
            Node *left_node, *right_node;
            InnerNode *parent = reinterpret_cast<InnerNode*>(parent_page->get_data());
            uint16_t left_slot;
            // determine which sibling to merge with
            if (child_slot == 0) {
                // there is only right sibling
                Swip &swip = parent->children[1];
                neighboor_page = &buffer_manager.fix_page(swip, true);
                left_slot = child_slot;
                left_node = reinterpret_cast<Node*>(child_page->get_data());
                right_node = reinterpret_cast<Node*>(neighboor_page->get_data());
            }
            else if (child_slot == parent->count-1) {
                // there is only left sibling
                Swip &swip = parent->children[parent->count-2];
                neighboor_page = &buffer_manager.fix_page(swip, true);
                left_slot = child_slot-1;
                left_node = reinterpret_cast<Node*>(neighboor_page->get_data());
                right_node = reinterpret_cast<Node*>(child_page->get_data());
            }
            else {
                Swip &left_swip = parent->children[child_slot-1],
                    &right_swip = parent->children[child_slot+1];
                BufferFrame *left_neighboor_page = &buffer_manager.fix_page(left_swip, true),
                            *right_neighboor_page = &buffer_manager.fix_page(right_swip, true);
                left_node = reinterpret_cast<Node*>(left_neighboor_page);
                right_node = reinterpret_cast<Node*>(right_neighboor_page);
                if (left_node->count > right_node->count) {
                    left_node = reinterpret_cast<Node*>(child_page->get_data());
                    left_slot = child_slot;
                    neighboor_page = right_neighboor_page;
                    buffer_manager.unfix_page(*left_neighboor_page, false);
                }
                else {
                    right_node = reinterpret_cast<Node*>(child_page->get_data());
                    left_slot = child_slot-1;
                    neighboor_page = left_neighboor_page;
                    buffer_manager.unfix_page(*right_neighboor_page, false);
                }
            }
            uint16_t total_count = left_node->count + right_node->count;
            uint16_t max_count = left_node->is_leaf() ? LeafNode::max_values() : InnerNode::max_children();
            if (total_count <= max_count) {
                // if parent has only these two children, we can merge them into parent.
                // This is possible only if the parent is the root node.
                if (parent->count == 2) {
                    if (left_node->is_leaf()) {
                        LeafNode *left_node_as_leaf = static_cast<LeafNode*>(left_node),
                                 *right_node_as_leaf = static_cast<LeafNode*>(right_node);
                        // turn parent into leaf node
                        parent->level = 0;
                        LeafNode *parent_as_leaf = reinterpret_cast<LeafNode*>(parent);
                        parent_as_leaf->count = left_node->count + right_node->count;
                        for (uint16_t i=0; i<left_node->count; ++i) {
                            parent_as_leaf->keys[i] = left_node_as_leaf->keys[i];
                            parent_as_leaf->values[i] = left_node_as_leaf->values[i];
                        }
                        for (uint16_t i=0; i<right_node->count; ++i) {
                            parent_as_leaf->keys[i+left_node->count] = right_node_as_leaf->keys[i];
                            parent_as_leaf->values[i+left_node->count] = right_node_as_leaf->values[i];
                        }
                    }
                    else {
                        InnerNode *left_node_as_inner = static_cast<InnerNode*>(left_node),
                                  *right_node_as_inner = static_cast<InnerNode*>(right_node);
                        --parent->level; // decrement parent's level
                        parent->count = left_node->count + right_node->count;
                        parent->keys[left_node->count-1] = parent->keys[0]; // move the old separator to the middle
                        for (uint16_t i=0; i<left_node->count-1; ++i) {
                            parent->keys[i] = left_node_as_inner->keys[i];
                            parent->children[i] = left_node_as_inner->children[i];
                        }
                        parent->children[left_node->count-1] = left_node_as_inner->children[left_node->count-1]; // move last child from the left node
                        for (uint16_t i=0; i<right_node->count-1; ++i) {
                            parent->keys[i+left_node->count] = right_node_as_inner->keys[i];
                            parent->children[i+left_node->count] = right_node_as_inner->children[i];
                        }
                        parent->children[parent->count-1] = right_node_as_inner->children[right_node->count-1]; // move last child from the right node
                    }
                    // TODO: push these pages to list of free pages
                    buffer_manager.unfix_page(*child_page, false);
                    buffer_manager.unfix_page(*neighboor_page, false);
                    return parent_page;
                }
                else {
                    // regular merge
                    if (left_node->is_leaf()) {
                        LeafNode *left_node_as_leaf = static_cast<LeafNode*>(left_node),
                                 *right_node_as_leaf = static_cast<LeafNode*>(right_node);
                        for (uint16_t i=0; i<right_node->count; ++i) {
                            left_node_as_leaf->keys[i+left_node->count] = right_node_as_leaf->keys[i];
                            left_node_as_leaf->values[i+left_node->count] = right_node_as_leaf->values[i];
                        }
                    }
                    else {
                        InnerNode *left_node_as_inner = static_cast<InnerNode*>(left_node),
                                 *right_node_as_inner = static_cast<InnerNode*>(right_node);
                        for (uint16_t i=0; i<right_node->count; ++i) {
                            left_node_as_inner->keys[i+left_node->count-1] = right_node_as_inner->keys[i];
                            left_node_as_inner->children[i+left_node->count] = right_node_as_inner->children[i];
                        } 
                    }
                    left_node->count += right_node->count;
                    parent->keys[left_slot] = parent->keys[left_slot+1];
                    // for (uint16_t i=left_slot+1; i+1<parent->count-1; ++i) {
                    //     parent->keys[i] = parent->keys[i+1];
                    //     parent->children[i] = parent->children[i+1];
                    // }
                    // parent->children[parent->count-2] = parent->children[parent->count-1];
                    // --parent->count;
                    parent->erase(left_slot+1);
                    buffer_manager.unfix_page(*parent_page, true);
                    if (left_slot==child_slot) {
                        buffer_manager.unfix_page(*neighboor_page, false);
                        // TODO: push neighboor_page to the list of free page
                        return child_page;
                    }
                    else {
                        // TODO: push child_page to the list of free page
                        buffer_manager.unfix_page(*child_page, false);
                        return neighboor_page;
                    }
                }
            }
            // merge not possible, rebalance
            KeyT &new_separator = parent->keys[left_slot];
            if (left_node->count > right_node->count) {
                // shift from left to right
                uint16_t to_shift = (left_node->count-right_node->count)/2; // shift half the difference
                if (left_node->is_leaf()) {
                    LeafNode *left_node_as_leaf = static_cast<LeafNode*>(left_node),
                             *right_node_as_leaf = static_cast<LeafNode*>(right_node);
                    // make space for keys and values from the left
                    for (uint16_t i=right_node->count; i>0;) {
                        --i;
                        right_node_as_leaf->keys[i+to_shift] = right_node_as_leaf->keys[i];
                        right_node_as_leaf->values[i+to_shift] = right_node_as_leaf->values[i];
                    }
                    for (uint16_t i=0, pos=left_node->count-to_shift; i<to_shift; ++i) {
                        right_node_as_leaf->keys[i] = left_node_as_leaf->keys[pos+i];
                        right_node_as_leaf->values[i] = left_node_as_leaf->values[pos+i];
                    }
                    new_separator = left_node_as_leaf->keys[left_node->count-to_shift-1];
                }
                else {
                    InnerNode *left_node_as_inner = static_cast<InnerNode*>(left_node),
                             *right_node_as_inner = static_cast<InnerNode*>(right_node);
                    // make space for keys and values from the left
                    for (uint16_t i=right_node->count; i>0;) {
                        --i;
                        right_node_as_inner->keys[i+to_shift] = right_node_as_inner->keys[i];
                        right_node_as_inner->children[i+to_shift] = right_node_as_inner->children[i];
                    }
                    right_node_as_inner->keys[to_shift-1] = parent->keys[left_slot]; // move old separator from parent
                    right_node_as_inner->children[to_shift-1] = left_node_as_inner->children[left_node->count-1]; // move rightmost child of left node
                    for (uint16_t i=0, pos=left_node->count-to_shift; i<to_shift-1; ++i) {
                        right_node_as_inner->keys[i] = left_node_as_inner->keys[pos+i];
                        right_node_as_inner->children[i] = left_node_as_inner->children[pos+i];
                    }
                    new_separator = right_node_as_inner->keys[left_node->count-to_shift-1];
                }
                left_node->count-=to_shift;
                right_node->count+=to_shift;
            }
            else {
                // shift from right to left
                uint16_t to_shift = (right_node->count-left_node->count)/2; // shift half the difference
                if (left_node->is_leaf()) {
                    LeafNode *left_node_as_leaf = static_cast<LeafNode*>(left_node),
                             *right_node_as_leaf = static_cast<LeafNode*>(right_node);
                    for (uint16_t i=0; i<to_shift; ++i) {
                        left_node_as_leaf->keys[i+left_node->count] = right_node_as_leaf->keys[i];
                        left_node_as_leaf->values[i+left_node->count] = right_node_as_leaf->values[i];
                    }
                    for (uint16_t i=0; i<right_node->count-to_shift; ++i) {
                        right_node_as_leaf->keys[i] = right_node_as_leaf->keys[i+to_shift];
                        right_node_as_leaf->values[i] = right_node_as_leaf->values[i+to_shift];
                    }
                    new_separator = left_node_as_leaf->keys[left_node->count+to_shift-1];
                }
                else {
                    InnerNode *left_node_as_inner = static_cast<InnerNode*>(left_node),
                             *right_node_as_inner = static_cast<InnerNode*>(right_node);
                    // move key from parent
                    left_node_as_inner->keys[left_node->count-1] = parent->keys[left_slot];
                    // move to_shift-1 keys and children from right to left
                    for (uint16_t i=0; i<to_shift-1; ++i) {
                        left_node_as_inner->keys[i+left_node->count] = right_node_as_inner->keys[i];
                        left_node_as_inner->children[i+left_node->count] = right_node_as_inner->children[i];
                    }
                    // move the last child from left to right
                    left_node_as_inner->children[left_node->count+to_shift-1] = right_node_as_inner->children[to_shift-1];
                    // make the new "first" key from left node the new separator
                    new_separator = right_node_as_inner->keys[to_shift-1];
                    for (uint16_t i=0; i<right_node->count-to_shift; ++i) {
                        right_node_as_inner->keys[i] = right_node_as_inner->keys[i+to_shift];
                        right_node_as_inner->children[i] = right_node_as_inner->children[i+to_shift];
                    }
                }
                left_node->count+=to_shift;
                right_node->count-=to_shift;
            }
            buffer_manager.unfix_page(*parent_page, true);
            if (target <= new_separator) {
                if (left_slot == child_slot) {
                    buffer_manager.unfix_page(*neighboor_page, true);
                    return child_page;
                }
                else {
                    buffer_manager.unfix_page(*child_page, true);
                    return neighboor_page;
                }
            }
            else {
                if (left_slot == child_slot) {
                    buffer_manager.unfix_page(*child_page, true);
                    return neighboor_page;
                }
                else {
                    buffer_manager.unfix_page(*neighboor_page, true);
                    return child_page;
                }
            }
        }
    
};
    
}  // namespace guidedresearch

#endif