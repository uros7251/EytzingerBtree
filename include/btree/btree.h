#ifndef INCLUDE_BTREE_H
#define INCLUDE_BTREE_H

#include <type_traits>
#include <forward_list>
#include "segment.h"
#include "buffer_manager/page_guard.h"
#include "btree/inner_node.h"
#include "btree/leaf_node.h"

namespace guidedresearch {

//------------------------------------------------------------------------------------------------
// The algorithmic part of this BTree implementation is heavily inspired
// by prof. Neumann's article https://databasearchitects.blogspot.com/2022/06/btreeoperations.html
//------------------------------------------------------------------------------------------------

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize, NodeLayout inner_layout = NodeLayout::SORTED, NodeLayout leaf_layout = NodeLayout::SORTED>
struct BTree : public Segment {
    using SharedPage = PageGuard<false>;
    using UniquePage = PageGuard<true>;
    using Node = guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize>;
    using InnerNode = guidedresearch::InnerNode<KeyT, ValueT, ComparatorT, PageSize, inner_layout>;
    // using leaf node designed for fast insertions
    using LeafNode = guidedresearch::LeafNode<KeyT, ValueT, ComparatorT, PageSize, leaf_layout, false>;

    constexpr static size_t kPageSize = PageSize; // expose page size to the outside

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

    /// @brief Return depth of the tree
    /// @return depth of the tree
    uint16_t depth() {
        SharedPage page(buffer_manager, root);
        auto *node = reinterpret_cast<Node*>(page->get_data());
        return node->level;
    }

    /// Lookup an entry in the tree.
    /// @param[in] key      The key that should be searched.
    /// @return             Whether the key was in the tree.
    __attribute__((noinline))
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
            parent_page = std::move(page);
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
            child_page = UniquePage(buffer_manager, swip);
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

    void merge(UniquePage& parent_page, UniquePage& child_page, uint16_t child_slot, KeyT target) {
        // hold tight, this is ~200 lines function :)

        // assumes child and parent pages are already fixed
        UniquePage neighboor_page(buffer_manager);
        // BufferFrame *neighboor_page;
        Node *left_node, *right_node;
        InnerNode &parent = *reinterpret_cast<InnerNode*>(parent_page->get_data());
        uint16_t left_slot;
        // determine which sibling to merge with
        if (child_slot == parent.first_sorted()) { // NEED FIX FOR EYTZINGER, we need leftmost key in sorted order
            // only right sibling exists
            Swip &swip = parent.children[parent.next_sorted(child_slot)];
            neighboor_page.fix(swip);
            left_slot = child_slot;
            left_node = reinterpret_cast<Node*>(child_page->get_data());
            right_node = reinterpret_cast<Node*>(neighboor_page->get_data());
        }
        else if (child_slot == parent.last_sorted()) { // NEED FIX FOR EYTZINGER, we need rightmost key in sorted order
            // only left sibling exists
            left_slot = parent.prev_sorted(child_slot); // NEED FIX FOR EYTZINGER, we need previous key in sorted order
            Swip &swip = parent.children[left_slot];
            neighboor_page.fix(swip);
            left_node = reinterpret_cast<Node*>(neighboor_page->get_data());
            right_node = reinterpret_cast<Node*>(child_page->get_data());
        }
        else {
            Swip &left_swip = parent.children[parent.prev_sorted(child_slot)], // NEED FIX FOR EYTZINGER
                &right_swip = parent.children[parent.next_sorted(child_slot)]; // NEED FIX FOR EYTZINGER
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
                left_slot = parent.prev_sorted(child_slot); // NEED FIX FOR EYTZINGER, we need previous key in sorted order
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
                    KeyT old_separator = parent.keys[parent.first_sorted()]; // NEED FIX FOR EYTZINGER
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
}  // namespace guidedresearch

#endif