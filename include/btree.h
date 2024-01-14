#ifndef INCLUDE_BTREE_H
#define INCLUDE_BTREE_H

#include <type_traits>
#include <forward_list>
#include <bit>
#include "segment.h"
#include "buffer_manager.h"
#include "swip.h"
#include "unique_page.h"
#include "shared_page.h"

namespace guidedresearch {

//------------------------------------------------------------------------------------------------
// The algorithmic part of this BTree implementation is heavily inspired
// by prof. Neumann's article https://databasearchitects.blogspot.com/2022/06/btreeoperations.html
//------------------------------------------------------------------------------------------------

enum class NodeLayout {
    SORTED,
    EYTZINGER
};

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize, NodeLayout layout = NodeLayout::SORTED>
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
                i += less(keys[i+half-1], key) * half; // hopefully compiler translates this to cmov
                n -= half; // ceil(n/2)
            }

            return ((!is_leaf() && i==Node::count-2u) || (is_leaf() && i==Node::count-1u)) && less(keys[i], key) ? // check if last key is less than key
                std::make_pair(i+1, false) : // return index one after last key
                std::make_pair(i, keys[i] == key); 
        }
    };

    struct InnerNode: public Node {
        /// Remark: it might be better to store the rightmost key to make nodes more uniform and algorithms more concise

        /// The capacity of a node.
        static constexpr uint32_t kCapacity = (PageSize-sizeof(Node)-sizeof(Swip))/(sizeof(KeyT)+sizeof(Swip));

        inline static constexpr uint32_t max_keys() { return kCapacity; }
        inline static constexpr uint32_t max_children() { return kCapacity+1; }

        /// The keys.
        KeyT keys[InnerNode::max_keys()]; 
        /// The values.
        Swip children[InnerNode::max_children()]; 

        /// Constructor.
        InnerNode() : Node(0, 0) { assert(reinterpret_cast<size_t>(&keys[0])==reinterpret_cast<size_t>(this)+Node::keys_offset); }
        InnerNode(uint16_t level) : Node(level, 0) {}
        InnerNode(const InnerNode &other) : Node(other) {
            copy();
        }
        
        /// Assignment operator
        InnerNode& operator=(const InnerNode &other) {
            copy(other);
            return *this;
        }

        /// this method can be virtual, but we avoid that on purpose
        /// @brief check if node is less than half full
        /// @return 
        bool merge_needed() const { return Node::count < max_children()/2; } 

        /// Insert a key.
        /// @param[in] key          The key that should be inserted.
        /// @param[in] child        The child that should be inserted.
        void insert_split(const KeyT &key, Swip child) {
            // for example, insert_split(3, 40)
            // BEFORE: keys = [1, 4, 9, ?], children = [8, 3, 13, 7, ?], count = 4
            // AFTER: keys = [1, 3, 4, 9], children = [8, 3, 40, 13, 7], count = 5

            assert(Node::count-1 != InnerNode::kCapacity);
            auto [index, _] = Node::lower_bound(key); // find index to insert key
            for (uint32_t i=Node::count-1; i>index; --i) {
                keys[i] = keys[i-1];
                children[i+1] = children[i];
            }
            keys[index] = key; // insert new separator
            children[index+1] = child; // insert child
            Node::count++;
        }

        /// @brief Erases a key and a child associated with the next key
        /// @param index Index of the key to be erased
        void erase(uint16_t index) {
            // for example, erase(1)
            // BEFORE: keys = [1, 3, 4, 9], children = [8, 3, 40, 13, 7], count = 5
            // AFTER: keys = [1, 4, 9, 9], children = [8, 3, 13, 7, 7], count = 4
            keys[index] = keys[index+1];
            ++index;
            assert(Node::count > 0);
            if (index < Node::count-1u) {
                for (uint32_t i=index; i+2u<Node::count; ++i) {
                    keys[i] = keys[i+1];
                    children[i] = children[i+1];
                }
                children[Node::count-2] = children[Node::count-1];
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
            InnerNode* new_node = new (buffer) InnerNode(Node::level);
            // old_node retains floor((count-1)/2) keys
            // new_node gets floor(count/2)-1 keys
            // move keys starting from index floor((count-1)/2)+1=floor((count+1)/2)=ceil(count/2) (skip separator), in total ceil((count-1)/2)-1=floor(count/2)-1
            // move children starting from index ceil(count/2), in total floor(count/2)
            for (auto i=0u, pos=(Node::count+1)/2; i<Node::count/2-1u; ++i) {
                new_node->keys[i] = keys[i+pos];
                new_node->children[i] = children[i+pos];
            }
            new_node->children[Node::count/2-1u] = children[Node::count-1];
            new_node->count = Node::count/2; // floor(count/2)
            Node::count -= Node::count/2; // ceil(count/2)
            return keys[Node::count-1]; // Node::count-1 = new number of keys in old_node
        }

        void merge(InnerNode &other, const KeyT &old_separator) {
            assert(Node::count > 0);
            keys[Node::count-1] = old_separator;
            for (uint16_t i=0; i<other.count-1u; ++i) {
                keys[i+Node::count] = other.keys[i];
                children[i+Node::count] = other.children[i];
            }
            children[Node::count+other.count-1u] = other.children[other.count-1u];
            Node::count += other.count;
        }

        /// Returns the keys.
        std::vector<KeyT> get_key_vector() {
            return std::vector<KeyT>(&keys[0], &keys[Node::count-1]);
        }

        static void rebalance(InnerNode &left, InnerNode &right, KeyT &separator) {
            if (left.count > right.count) {
                uint16_t to_shift = (left.count-right.count)/2; // shift half the difference
                // make space for keys and values from the left node
                for (uint16_t i=right.count; i>0;) {
                    --i;
                    right.keys[i+to_shift] = right.keys[i];
                    right.children[i+to_shift] = right.children[i];
                }
                right.keys[to_shift-1] = separator; // move old separator from parent
                right.children[to_shift-1] = left.children[left.count-1]; // move rightmost child of left node
                for (uint16_t i=0, pos=left.count-to_shift; i<to_shift-1; ++i) {
                    right.keys[i] = left.keys[pos+i];
                    right.children[i] = left.children[pos+i];
                }
                separator = left.keys[left.count-to_shift-1];
                left.count -= to_shift;
                right.count += to_shift;
            }
            else {
                uint16_t to_shift = (right.count-left.count)/2; // shift half the difference
                // move separator from parent
                left.keys[left.count-1] = separator;
                // move to_shift-1 keys and children from right to left
                for (uint16_t i=0; i<to_shift-1; ++i) {
                    left.keys[i+left.count] = right.keys[i];
                    left.children[i+left.count] = right.children[i];
                }
                // move the last child from right to left
                left.children[left.count+to_shift-1] = right.children[to_shift-1];
                separator = right.keys[to_shift-1];
                for (uint16_t i=0; i<right.count-to_shift; ++i) {
                    right.keys[i] = right.keys[i+to_shift];
                    right.children[i] = right.children[i+to_shift];
                }
                left.count += to_shift;
                right.count -= to_shift;
            }
        }

    private:
        void copy(const InnerNode &other) {
            Node::operator=(other);
            for (uint32_t i=0; i<Node::count-1u; ++i) {
                keys[i] = other.keys[i];
                children[i] = other.children[i];
            }
            children[Node::count-1u] = other.children[Node::count-1u];
        }
    };

    struct LeafNode: public Node {
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

    /// The root.
    Swip root;

    /// Constructor.
    BTree(uint16_t segment_id, BufferManager &buffer_manager)
        : Segment(segment_id, buffer_manager), root(allocate_page()) {
        
        UniquePage page(buffer_manager, root);
        new (page->get_data()) LeafNode();
        page.mark_dirty();
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
        SharedPage page(buffer_manager, root), parent_page(buffer_manager);
        auto *node = reinterpret_cast<Node*>(page->get_data());
        while (!node->is_leaf()) {
            // traverse inner nodes
            InnerNode *inner_node = reinterpret_cast<InnerNode*>(node);
            auto [index, match] = inner_node->lower_bound(key);
            Swip swip = inner_node->children[index];
            parent_page = page;
            // lock coupling
            page.fix(swip);
            node = reinterpret_cast<Node*>(page->get_data());
        }
        // leaf level
        auto* leaf_node = reinterpret_cast<LeafNode*>(node);
        auto [index, match] = leaf_node->lower_bound(key);
        std::optional<ValueT> result;
        if (match) result = {leaf_node->values[index]};
        return result;
    }

    template<typename T>
    requires std::is_base_of_v<Node, T>
    void handle_root(UniquePage& child_page, UniquePage& parent_page) {
        // create a new node and copy the contents of the child node to it
        // modify initial child node such that it's now an InnerNode with only one child - newly created node
        Swip swip = allocate_page();
        UniquePage new_page(buffer_manager, swip);
        T *new_node = new (new_page->get_data()) T(); 
        *new_node = *reinterpret_cast<T*>(child_page->get_data()); // copy
        InnerNode *old_node = new (child_page->get_data()) InnerNode(new_node->level+1);
        old_node->children[0]=swip;
        old_node->count=1;
        parent_page = std::move(child_page);
        child_page = std::move(new_page);
    }

    /// @brief Splits the child page node and updates the child page through reference if neccessary
    /// @tparam T 
    /// @param child_page lvalue reference to the child page
    /// @param parent_page lvalue reference to the parent page
    /// @param target key to be inserted which determines if the child page will be updated 
    template<typename T>
    requires std::is_base_of_v<Node, T>
    void split_node(UniquePage& child_page, UniquePage& parent_page, KeyT target) {
        // split
        T *node = reinterpret_cast<T*>(child_page->get_data());
        InnerNode *parent_node = reinterpret_cast<InnerNode*>(parent_page->get_data());
        Swip swip = allocate_page();
        UniquePage page(buffer_manager, swip);
        KeyT new_separator = node->split(page->get_data());
        // mark pages as dirty
        child_page.mark_dirty();
        page.mark_dirty();
        // insert new separator into parent node
        parent_node->insert_split(new_separator, swip);
        parent_page.mark_dirty();
        if (target > new_separator) {
            child_page = std::move(page);
        }
    }

    /// Inserts a new entry into the tree.
    /// @param[in] key      The key that should be inserted.
    /// @param[in] value    The value that should be inserted.
    void insert(const KeyT &key, const ValueT &value) {
        UniquePage child_page(buffer_manager, root), parent_page(buffer_manager);
        Node *child = reinterpret_cast<Node*>(child_page->get_data());
        if (child->is_leaf()) {
            LeafNode *child_as_leaf = reinterpret_cast<LeafNode*>(child);
            if (child_as_leaf->count == LeafNode::max_values()) {
                handle_root<LeafNode>(child_page, parent_page);
                child = reinterpret_cast<Node*>(child_page->get_data());
            }
        }
        else {
            InnerNode *child_as_inner = reinterpret_cast<InnerNode*>(child);
            if (child_as_inner->count == InnerNode::max_children()) {
                handle_root<InnerNode>(child_page, parent_page);
                child = reinterpret_cast<Node*>(child_page->get_data());
            }
        }
        while (!child->is_leaf()) {
            if (child->count == InnerNode::max_children()) {
                split_node<InnerNode>(child_page, parent_page, key);
            }
            InnerNode *child_as_inner = reinterpret_cast<InnerNode*>(child_page->get_data());
            auto [index, match] = child_as_inner->lower_bound(key);
            Swip &swip = child_as_inner->children[index]; // it's crucial that swip is a reference
            parent_page = std::move(child_page);
            child_page.fix(swip);
            child = reinterpret_cast<Node*>(child_page->get_data());
        }
        // at the bottom
        if (child->count == LeafNode::max_values()) {
            split_node<LeafNode>(child_page, parent_page, key);
        }
        
        LeafNode *child_as_leaf = reinterpret_cast<LeafNode*>(child_page->get_data());
        child_as_leaf->insert(key, value);
        child_page.mark_dirty();
    }

    /// Erase a key.
    void erase(const KeyT &key) {
        UniquePage child_page(buffer_manager, root), parent_page(buffer_manager);
        Node *child = reinterpret_cast<Node*>(child_page->get_data());
        InnerNode *parent = nullptr;
        while (!child->is_leaf()) {
            parent_page = std::move(child_page);
            parent = static_cast<InnerNode*>(child);
            auto [index, match] = parent->lower_bound(key);
            Swip &swip = parent->children[index]; // it's crucial that swip is a reference because fix() might change it
            child_page.fix(swip);
            child = reinterpret_cast<Node*>(child_page->get_data());
            if (child->is_leaf()) {
                LeafNode *child_as_leaf = static_cast<LeafNode*>(child);
                if (child_as_leaf->merge_needed()) {
                    // this call will unfix the pages that need to be unfixed
                    merge(parent_page, child_page, index, key);
                }
            }
            else {
                InnerNode *child_as_inner = static_cast<InnerNode*>(child);
                if (child_as_inner->merge_needed()) {
                    merge(parent_page, child_page, index, key);
                }
            }
            child = reinterpret_cast<Node*>(child_page->get_data());
        }
        LeafNode *leaf = static_cast<LeafNode*>(child);
        leaf->erase(key);
    }

    private:

    void merge(UniquePage& parent_page, UniquePage& child_page, uint16_t child_slot, KeyT target) {
        // hold tight, this is ~200 lines function :)

        // assumes child and parent pages are already fixed
        UniquePage neighboor_page(buffer_manager);
        // BufferFrame *neighboor_page;
        Node *left_node, *right_node;
        InnerNode &parent = *reinterpret_cast<InnerNode*>(parent_page->get_data());
        uint16_t left_slot;
        // determine which sibling to merge with
        if (child_slot == 0) { // NEED FIX FOR EYTZINGER, we need leftmost key in sorted order
            // only right sibling exists
            Swip &swip = parent.children[1];
            neighboor_page.fix(swip);
            left_slot = child_slot;
            left_node = reinterpret_cast<Node*>(child_page->get_data());
            right_node = reinterpret_cast<Node*>(neighboor_page->get_data());
        }
        else if (child_slot == parent.count-1) { // NEED FIX FOR EYTZINGER, we need rightmost key in sorted order
            // only left sibling exists
            Swip &swip = parent.children[parent.count-2];
            neighboor_page.fix(swip);
            left_slot = child_slot-1; // NEED FIX FOR EYTZINGER, we need previous key in sorted order
            left_node = reinterpret_cast<Node*>(neighboor_page->get_data());
            right_node = reinterpret_cast<Node*>(child_page->get_data());
        }
        else {
            Swip &left_swip = parent.children[child_slot-1],
                &right_swip = parent.children[child_slot+1];
            UniquePage left_neighboor_page(buffer_manager, left_swip),
                    right_neighboor_page(buffer_manager, right_swip);
            left_node = reinterpret_cast<Node*>(left_neighboor_page->get_data());
            right_node = reinterpret_cast<Node*>(right_neighboor_page->get_data());
            if (left_node->count > right_node->count) {
                left_node = reinterpret_cast<Node*>(child_page->get_data());
                left_slot = child_slot;
                neighboor_page = std::move(right_neighboor_page);
            }
            else {
                right_node = reinterpret_cast<Node*>(child_page->get_data());
                left_slot = child_slot-1; // NEED FIX FOR EYTZINGER, we need previous key in sorted order
                neighboor_page = std::move(left_neighboor_page);
            }
        }
        uint16_t total_count = left_node->count + right_node->count;
        uint16_t max_count = left_node->is_leaf() ? LeafNode::max_values() : InnerNode::max_children();
        if (total_count <= max_count) {
            // if parent has only these two children, we can merge them into parent.
            // This is possible only if the parent is the root node.
            if (parent.count == 2) {
                if (left_node->is_leaf()) {
                    LeafNode &left_node_as_leaf = *static_cast<LeafNode*>(left_node),
                                &right_node_as_leaf = *static_cast<LeafNode*>(right_node);
                    // turn parent into leaf node
                    LeafNode *parent_as_leaf = new (&parent) LeafNode();
                    parent_as_leaf->merge(left_node_as_leaf);
                    parent_as_leaf->merge(right_node_as_leaf);
                }
                else {
                    InnerNode &left_node_as_inner = *static_cast<InnerNode*>(left_node),
                                &right_node_as_inner = *static_cast<InnerNode*>(right_node);
                    KeyT old_separator = parent.keys[0];
                    parent = left_node_as_inner; // copy
                    parent.merge(right_node_as_inner, old_separator);
                }
                // enable later reuse of these frames
                deallocate_page(child_page.bf_ptr());
                deallocate_page(neighboor_page.bf_ptr());
                child_page = std::move(parent_page);
            }
            else {
                // regular merge
                if (left_node->is_leaf()) {
                    LeafNode &left_node_as_leaf = *static_cast<LeafNode*>(left_node),
                                &right_node_as_leaf = *static_cast<LeafNode*>(right_node);
                    left_node_as_leaf.merge(right_node_as_leaf);
                }
                else {
                    InnerNode *left_node_as_inner = static_cast<InnerNode*>(left_node),
                                *right_node_as_inner = static_cast<InnerNode*>(right_node);
                    left_node_as_inner->merge(*right_node_as_inner, parent.keys[left_slot]);
                }
                parent.erase(left_slot);
                parent_page.mark_dirty();
                if (left_slot==child_slot) {
                    deallocate_page(neighboor_page.bf_ptr());
                }
                else {
                    deallocate_page(child_page.bf_ptr());
                    child_page = std::move(neighboor_page);
                }
            }
            return;
        }
        // merge not possible, rebalance
        KeyT &separator = parent.keys[left_slot];
        if (left_node->is_leaf()) {
            LeafNode &left_node_as_leaf = *static_cast<LeafNode*>(left_node),
                     &right_node_as_leaf = *static_cast<LeafNode*>(right_node);
            LeafNode::rebalance(left_node_as_leaf, right_node_as_leaf, separator);
        }
        else {
            InnerNode &left_node_as_inner = *static_cast<InnerNode*>(left_node),
                      &right_node_as_inner = *static_cast<InnerNode*>(right_node);
            InnerNode::rebalance(left_node_as_inner, right_node_as_inner, separator);
        }
        neighboor_page.mark_dirty();
        child_page.mark_dirty();
        parent_page.mark_dirty();
        if ((target <= separator && left_slot != child_slot) || (target > separator && left_slot == child_slot)) {
            child_page = std::move(neighboor_page);
        }
    }
    
};

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct BTree<KeyT, ValueT, ComparatorT, PageSize, NodeLayout::EYTZINGER> : public Segment {

    using Node = BTree<KeyT,ValueT,ComparatorT,PageSize>::Node;
    
    struct InnerNode: public Node {
        struct Iterator {
            uint32_t k;
            uint16_t n;

            Iterator() = delete;
            Iterator(uint32_t k, uint16_t size) :k(k), n(size) {}

            bool operator==(const Iterator &other) const = default;
            uint32_t operator*() { return k; }
            Iterator& operator++() {
                if (2*k+1<n) {
                    // go right
                    k = 2*k+1;
                    // then left as much as possible
                    while (2*k<n) k = 2*k;
                }
                else {
                    k = k >> (std::countr_one(k)+1); // climb up until you're the right child and then climb up one more time
                }
                return *this;
            }

            Iterator& operator--() {
                if (2*k<n) {
                    // go left
                    k = 2*k;
                    // then right as much as possible
                    while (2*k+1<n) k = 2*k+1;
                }
                else {
                    k = k >> (std::countr_zero(k)+1); // climb up until you're the left child and then climb up one more time
                }
                return *this;
            }

            static Iterator begin(uint16_t size) { assert(size); return Iterator((1u << (32 - std::countl_zero(size-1u))) >> 1, size); } // greatest power of 2 less than n
            static Iterator end(uint16_t size)  { return Iterator(0u, size); }
            static Iterator rbegin(uint16_t size) { return Iterator(((1u << (16 - std::countl_zero(size))) >> 1) - 1u, size); } // (greatest power of 2 not bigger than n) - 1
            static Iterator rend(uint16_t size)  { return Iterator(0u, size); }
        };

        template<bool decreasing>
        struct InOrderIterator {
            /// current index
            uint32_t k;
            /// last index in in-order traversal
            uint32_t last;
            /// smallest invalid index
            uint16_t n;

            InOrderIterator() = delete;
            InOrderIterator(uint16_t size)
                :k((decreasing ? get_last(size) : get_first(size)) << 1),
                last(decreasing ? get_first(size) : get_last(size)),
                n(size) {}
            InOrderIterator(uint16_t size, uint32_t start)
                :k(start),
                last(get_last(size)),
                n(size) {}
            
            bool end() { return k == last; }
            uint32_t next() {
                assert(!end());
                if (decreasing) {
                    if (2*k<n) {
                        // go left
                        k = 2*k;
                        // then right as much as possible
                        while (2*k+1<n) k = 2*k+1;
                        return k;
                    }
                    k = k >> (std::countr_zero(k)+1); // climb up until you're the left child and then climb up one more time
                    return k;
                }
                else {
                    if (2*k+1<n) {
                        // go right
                        k = 2*k+1;
                        // then left as much as possible
                        while (2*k<n) k = 2*k;
                        return k;
                    }
                    k = k >> (std::countr_one(k)+1); // climb up until you're the right child and then climb up one more time
                    return k;
                }
            }

            private:
            static uint32_t get_first(uint16_t n) {
                assert(n);
                return n > 1u ?
                    1u << (31u - std::countl_zero(n-1u)) : // greatest power of 2 less than n
                    0u;
            }
            static uint32_t get_last(uint16_t n) {
                assert(n);
                return n > 1u ?
                    (1u << (15u - std::countl_zero(n)))-1u : // (greatest power of 2 not bigger than n) - 1
                    0u;
            }
        };
        
        /// The capacity of a node.
        static constexpr uint32_t kCapacity = (PageSize-sizeof(Node))/(sizeof(KeyT)+sizeof(Swip))-1;

        inline static constexpr uint32_t max_keys() { return kCapacity; }
        inline static constexpr uint32_t max_children() { return kCapacity+1; }

        /// The keys. Stored beginning with index 1.
        KeyT keys[InnerNode::max_keys()+1]; // we leave first empty 
        /// The children. Child associated with the key at index i is stored at position i. The remaining child is stored at the 0th position.
        Swip children[InnerNode::max_children()]; 

        /// Constructor.
        InnerNode() : Node(0, 0) { }
        InnerNode(uint16_t level) : Node(level, 0) {}

        /// this method can be virtual, but we avoid that on purpose
        /// @brief check if node is less than half full
        /// @return 
        bool merge_needed() const { return Node::count < max_children()/2; } 

        std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
            // explained at https://en.algorithmica.org/hpc/data-structures/binary-search/

            if (Node::count <= 1) return {0u, false}; // no keys
            
            ComparatorT less;
            uint32_t i=1;
            while (i < Node::count) {
                i = 2*i + less(keys[i], key);
            }

            // recover index
            i >>= std::countr_one(i)+1;

            return i ? // check if last key is less than key
                std::make_pair(i, keys[i] == key) :
                std::make_pair(0u, false); 
        }

        /// Insert a key.
        /// @param[in] key          The key that should be inserted.
        /// @param[in] child        The child that should be inserted.
        void insert_split(const KeyT &key, Swip child) {
            // for example, insert_split(3, 40)
            // BEFORE: in_order(keys) = [1, 4, 9, ?], children = [8, 3, 13, 7, ?], count = 4
            // AFTER: in_order(keys) = [1, 3, 4, 9], children = [8, 3, 40, 13, 7], count = 5
            // idea is to create a temporary copy of keys and children, and then recreate a layout using this copy and new pair (key, child) as source
            
            assert(Node::count > 0); // inner node should have at least one child before inserting splits
            assert(Node::count-1 != InnerNode::kCapacity); // more place available

            // trivial case with no keys
            if (Node::count == 1u) {
                // no keys
                keys[1] = key;
                children[1] = children[0];
                children[0] = child;
                ++Node::count;
                return;
            }

            // --------------------------------
            Iterator lead(Node::count, Node::count+1);
            Iterator lag = lead;
            ComparatorT less;
            bool forward;
            if (lead == Iterator::rbegin(Node::count+1)) forward = false;
            else if (lead == Iterator::begin(Node::count+1)) forward = true;
            else {
                ++lead;
                if (less(keys[*lead], key)) {
                    forward = true;
                    keys[*lag] = keys[*lead];
                    children[*lag] = children[*lead];
                    lag = lead;
                }
                else {
                    forward = false;
                    --lead;
                }
            }
            assert(lag == lead);
            if (forward) {
                ++lead;
                for (auto end = Iterator::end(Node::count+1); lead != end && less(keys[*lead], key); lag=lead, ++lead) {
                    keys[*lag] = keys[*lead];
                    children[*lag] = children[*lead];
                }
            }
            else {
                --lead;
                for (auto end = Iterator::rend(Node::count+1); lead != end && less(key, keys[*lead]); lag=lead, --lead) {
                    keys[*lag] = keys[*lead];
                    children[*lag] = children[*lead];
                }
            }
            keys[*lag] = key;
            // update children
            lead = lag, ++lead;
            if (lead == Iterator::end(Node::count+1)) {
                children[*lag] = children[0];
                children[0] = child;
            }
            else {
                children[*lag] = children[*lead];
                children[*lead] = child;
            }
            ++Node::count;

            // --------------------------------
            /* // create copies of keys and swips
            std::vector<KeyT> key_vec(&keys[0], &keys[Node::count]);
            std::vector<Swip> swip_vec(&children[0], &children[Node::count]);

            // create iterators which yield indices over old and new key array size
            auto lag = Iterator::begin(Node::count), lead = Iterator::begin(Node::count+1);

            ComparatorT less;
            uint32_t i;
            bool key_inserted = false;

            // in-order traversal
            for (auto end = Iterator::end(Node::count); lag != end; ++lag, ++lead) {
                i = *lead;
                auto j = *lag;
                children[i-1] = swip_vec[j-1];
                // if key is not already inserted and it's less than the current key
                if (!key_inserted && less(key, keys[j])) {
                    keys[i] = key;
                    ++lead;
                    i = *lead;
                    children[i-1] = child;
                    key_inserted = true;
                }
                keys[i] = key_vec[j];
            }
            if (lead == Iterator::end(Node::count+1)) {
                children[Node::count] = swip_vec[Node::count-1];
            }
            else {
                keys[*lead] = key;
                children[(*lead)-1] = swip_vec[Node::count-1];
                children[Node::count] = child;
            }
            ++Node::count; */
        }

        /// @brief Erases a key and a child associated with the next key
        /// @param index Index of the key to be erased
        void erase(uint16_t index) {
            // for example, erase(1)
            // BEFORE: keys = [1, 3, 4, 9], children = [8, 3, 40, 13, 7], count = 5
            // AFTER: keys = [1, 4, 9, 9], children = [8, 3, 13, 7, 7], count = 4

            // -------------------------------
            Iterator lag(index, Node::count);
            Iterator lead = lag;
            ComparatorT less;
            // if index comes last during in-order traversal
            if (lead == Iterator::rbegin(Node::count)) {
                children[0] = children[index];
            }
            else {
                ++lead;
                children[*lead] = children[*lag];
            }
            if (less(keys[index], keys[Node::count-1])) {
                // go right
                for (; *lag != Node::count-1u; lag=lead, ++lead) {
                    keys[*lag] = keys[*lead];
                    children[*lag] = children[*lead];
                }
            }
            else {
                // go left
                lead = lag; --lead;
                for (; *lag != Node::count-1u; lag=lead, --lead) {
                    keys[*lag] = keys[*lead];
                    children[*lag] = children[*lead];
                }
            }
            --Node::count;
            // -------------------------------
            /* // create copies of keys and swips
            std::vector<KeyT> key_vec(&keys[0], &keys[Node::count]);
            std::vector<Swip> swip_vec(&children[0], &children[Node::count]);

            // create iterators which yield indices over old and new key array size
            InOrderIterator<false> old_iter(Node::count), new_iter(Node::count-1);

            while (!new_iter.end()) {
                auto i = new_iter.next(), j = old_iter.next();
                children[i-1] = swip_vec[j-1];
                if (j == index) {
                    keys[i] = key_vec[old_iter.next()];
                }
                else {
                    keys[i] = key_vec[j];
                }
            }
            if (old_iter.end()) {
                children[Node::count-2] = swip_vec[Node::count-1];
            }
            else {
                children[Node::count-2] = swip_vec[old_iter.next()-1];
            }
            --Node::count; */
        }

        /// Split the node.
        /// @param[in] buffer       The buffer for the new page.
        /// @return                 The separator key.
        KeyT split(char* buffer) {
            InnerNode* new_node = new (buffer) InnerNode(Node::level);
            // new_node->count = Node::count;
            // old_node retains floor((count-1)/2) keys
            // new_node gets floor(count/2)-1 keys
            uint16_t left_node_count = (Node::count+1)/2, right_node_count = Node::count/2;
            assert(left_node_count+right_node_count==Node::count);
            // -------------------------------
            Iterator old_it = Iterator::rbegin(Node::count), left_it = Iterator::rbegin(left_node_count), right_it = Iterator::rbegin(right_node_count);
            new_node->children[0] = children[0];
            for (auto i=0; i<right_node_count-1; ++i, --old_it, --right_it) {
                new_node->keys[*right_it] = keys[*old_it];
                new_node->children[*right_it] = children[*old_it];
            }
            assert(right_it == Iterator::rend(right_node_count));
            
            auto separator = keys[*old_it];
            auto separator_swip = children[*old_it];
            --old_it;
            for (auto i=0; i<left_node_count-1; ++i, --old_it, --left_it) {
                keys[*left_it] = keys[*old_it];
                children[*left_it] = children[*old_it];
            }
            assert(left_it == Iterator::rend(left_node_count));
            assert(old_it == Iterator::rend(Node::count));
            children[0] = separator_swip;
            Node::count = left_node_count;
            new_node->count = right_node_count;
            auto rpids = new_node->get_sorted_pids(), lpids = get_sorted_pids();
            return separator;
            // -------------------------------
            /* InOrderIterator<false> old_node_iter(Node::count), left_node_iter(left_node_count), right_node_iter(right_node_count);
            // write smaller half to the second half of the new node
            uint32_t j, k;
            while (!left_node_iter.end()) {
                assert(!old_node_iter.end());
                j = right_node_count+left_node_iter.next(), k = old_node_iter.next();
                new_node->keys[j] = keys[k];
                new_node->children[j-1] = children[k-1];
            }
            k = old_node_iter.next();
            auto separator = keys[k];
            new_node->children[Node::count-1] = children[k-1];
            // write bigger half to the first half of the new node
            while (!right_node_iter.end()) {
                j = right_node_iter.next(), k = old_node_iter.next();
                new_node->keys[j] = keys[k];
                new_node->children[j-1] = children[k-1]; 
            }
            assert(old_node_iter.end());
            new_node->children[right_node_count-1] = children[Node::count-1];
            // copy smaller half back to the original node
            for (auto i=1u; i<left_node_count; ++i) {
                keys[i] = new_node->keys[i+right_node_count];
                children[i-1] = new_node->children[i+right_node_count-1];
            }
            children[left_node_count-1] = new_node->children[Node::count-1];
            // update counts
            Node::count = left_node_count;
            new_node->count = right_node_count;
            auto rpids = new_node->get_sorted_pids(), lpids = get_sorted_pids();
            return separator; */
        }

        void merge(InnerNode &other, const KeyT &old_separator) {
            // -------------------------------
            Iterator left_it = Iterator::begin(Node::count), right_it = Iterator::begin(other.count), merge_it = Iterator::begin(Node::count+other.count);
            // rearrange children in left node
            for (auto i=0; i<Node::count-1; ++i, ++left_it, ++merge_it) {
                keys[*merge_it] = keys[*left_it];
                children[*merge_it] = children[*left_it];
            }
            assert(left_it == Iterator::end(Node::count));
            // insert old separator
            keys[*merge_it] = old_separator;
            children[*merge_it] = children[0];
            ++merge_it;
            // insert children from right node
            for (auto i=0; i<other.count-1; ++i, ++right_it, ++merge_it) {
                keys[*merge_it] = other.keys[*right_it];
                children[*merge_it] = other.children[*right_it];
            }
            assert(right_it == Iterator::end(other.count));
            children[0] = other.children[0];
            // update counts
            Node::count += other.count;
            other.count = 0;
            // -------------------------------
            /* InOrderIterator<false> merged_iter(Node::count+other.count), left_iter(Node::count), right_iter(other.count);
            uint32_t i, j;
            // make a copy of the last child of left node because it might be overwritten
            Swip last_left = children[Node::count-1];
            // rearrange children in left node
            while (!left_iter.end()) {
                i = merged_iter.next(), j = left_iter.next();
                keys[i] = keys[j];
                children[i-1] = children[j-1];
            }
            // insert old separator
            i = merged_iter.next();
            keys[i] = old_separator;
            children[i-1] = last_left;
            // insert children from right node
            while (!right_iter.end()) {
                i = merged_iter.next(), j = right_iter.next();
                keys[i] = other.keys[j];
                children[i-1] = other.children[j-1];
            }
            assert(merged_iter.end());
            children[Node::count+other.count-1] = other.children[other.count-1];
            // update counts
            Node::count += other.count;
            other.count = 0; */
        }

        /// Returns the keys.
        std::vector<KeyT> get_key_vector() {
            return std::vector<KeyT>(&keys[1], &keys[Node::count]);
        }

        std::vector<uint64_t> get_sorted_pids() {
            std::vector<uint64_t> sorted_pids;
            sorted_pids.reserve(Node::count);
            for (auto it = Iterator::begin(Node::count), end = Iterator::end(Node::count); it != end; ++it) {
                sorted_pids.push_back(children[*it].asPageID());
            }
            sorted_pids.push_back(children[0].asPageID());
            return sorted_pids;
        }

        /// Returns the keys in sorted order
        std::vector<KeyT> get_sorted_keys() {
            std::vector<KeyT> sorted_keys;
            sorted_keys.reserve(Node::count-1);
            for (auto it = Iterator::begin(Node::count), end = Iterator::end(Node::count); it != end; ++it) {
                sorted_keys.push_back(keys[*it]);
            }
            return sorted_keys;
        }

        static void rebalance(InnerNode &left, InnerNode &right, KeyT &separator) {
            if (left.count > right.count) {
                // shift right
                uint16_t to_shift = (left.count-right.count)/2; // shift half the difference
                Iterator old_right_it = Iterator::rbegin(right.count), new_right_it = Iterator::rbegin(right.count+to_shift),
                         old_left_it = Iterator::rbegin(left.count), new_left_it = Iterator::rbegin(left.count-to_shift);
                // make space for keys and children from left node
                for (int i=0; i<right.count-1; ++i, --old_right_it, --new_right_it) {
                    right.keys[*new_right_it] = right.keys[*old_right_it];
                    right.children[*new_right_it] = right.children[*old_right_it];
                }
                assert(old_right_it == Iterator::rend(right.count));
                // insert separator and last child of the left node
                right.keys[*new_right_it] = separator;
                right.children[*new_right_it] = left.children[0];
                --new_right_it;
                // insert keys and children from left node
                for (int i=0; i<to_shift-1; ++i, --old_left_it, --new_right_it) {
                    right.keys[*new_right_it] = left.keys[*old_left_it];
                    right.children[*new_right_it] = left.children[*old_left_it];
                }
                assert(new_right_it == Iterator::rend(right.count+to_shift));
                // update separator
                separator = left.keys[*old_left_it];
                left.children[0] = left.children[*old_left_it];
                --old_left_it;
                // rearrange keys and children in left node
                for (int i=0; i<left.count-to_shift-1; ++i, --old_left_it, --new_left_it) {
                    left.keys[*new_left_it] = left.keys[*old_left_it];
                    left.children[*new_left_it] = left.children[*old_left_it];
                }
                assert(old_left_it == Iterator::rend(left.count));
                assert(new_left_it == Iterator::rend(left.count-to_shift));
                // update counts
                left.count -= to_shift;
                right.count += to_shift;
            }
            else {
                // shift left
                uint16_t to_shift = (right.count-left.count)/2; // shift half the difference
                Iterator old_right_it = Iterator::begin(right.count), new_right_it = Iterator::begin(right.count-to_shift),
                         old_left_it = Iterator::begin(left.count), new_left_it = Iterator::begin(left.count+to_shift);
                // make space for keys and children from left node
                for (int i=0; i<left.count-1; ++i, ++old_left_it, ++new_left_it) {
                    left.keys[*new_left_it] = left.keys[*old_left_it];
                    left.children[*new_left_it] = left.children[*old_left_it];
                }
                assert(old_left_it == Iterator::end(left.count));
                // insert separator and last child of the left node
                left.keys[*new_left_it] = separator;
                left.children[*new_left_it] = left.children[0];
                ++new_left_it;
                // insert keys and children from right node
                for (int i=0; i<to_shift-1; ++i, ++old_right_it, ++new_left_it) {
                    left.keys[*new_left_it] = right.keys[*old_right_it];
                    left.children[*new_left_it] = right.children[*old_right_it];
                }
                assert(new_left_it == Iterator::end(left.count+to_shift));
                // update separator
                separator = right.keys[*old_right_it];
                left.children[0] = right.children[*old_right_it];
                ++old_right_it;
                // rearrange keys and children in right node
                for (int i=0; i<right.count-to_shift-1; ++i, ++old_right_it, ++new_right_it) {
                    right.keys[*new_right_it] = right.keys[*old_right_it];
                    right.children[*new_right_it] = right.children[*old_right_it];
                }
                assert(old_right_it == Iterator::end(right.count));
                assert(new_right_it == Iterator::end(right.count-to_shift));
                // update counts
                left.count += to_shift;
                right.count -= to_shift;
            }
        }
    };

};
}  // namespace guidedresearch

#endif