#ifndef INCLUDE_INNER_NODE_H
#define INCLUDE_INNER_NODE_H

#include <bit>
#include <buffer_manager/swip.h>
#include <btree/node.h>

namespace guidedresearch {

enum class NodeLayout {
    SORTED,
    EYTZINGER
};

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize, NodeLayout layout>
struct InnerNode {
    KeyT *keys;
    Swip *children;
    InnerNode();
    InnerNode(uint16_t level);
    InnerNode(const InnerNode &other);
    
    /// Assignment operator
    InnerNode& operator=(const InnerNode &other);
    bool merge_needed();
    void insert_split(const KeyT &key, Swip child);
    void erase(uint16_t index);
    void merge(InnerNode &other, const KeyT &old_separator);
    std::vector<KeyT> get_key_vector();

    uint32_t first_sorted();
    uint32_t last_sorted();
    uint32_t next_sorted(uint32_t);
    uint32_t prev_sorted(uint32_t);

    static void rebalance(InnerNode &left, InnerNode &right, KeyT &separator);
};

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct InnerNode<KeyT, ValueT, ComparatorT, PageSize, NodeLayout::SORTED> : public Node<KeyT, ValueT, ComparatorT, PageSize> {
    using Node = guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize>;
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

    uint32_t first_sorted() { return 0; }
    uint32_t last_sorted() { return Node::count-1; }
    uint32_t next_sorted(uint32_t i) { return i+1; }
    uint32_t prev_sorted(uint32_t i) { return i-1; }

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

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct InnerNode<KeyT, ValueT, ComparatorT, PageSize, NodeLayout::EYTZINGER> : public Node<KeyT, ValueT, ComparatorT, PageSize> {
    using Node = guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize>;

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
    InnerNode(const InnerNode &other) = default;

    InnerNode& operator=(const InnerNode& other) = default;

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

        // idea is to logically insert a (key, child) pair into the node and then move it to the right place.
        // we start at the index Node::count (which was previously empty) and first determine if to-be-inserted
        // key should be placed to the left or to the right of it in ascending order (note: if i<j, it might be that 
        // in_order(i) > in_order(j))

        Iterator lead(Node::count, Node::count+1);
        Iterator lag = lead;
        ComparatorT less;
        bool forward;
        // determine if we should go forward or backward in in-order traversal
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
    }

    /// @brief Erases a key and a child associated with the next key
    /// @param index Index of the key to be erased
    void erase(uint16_t index) {
        // for example, erase(1)
        // BEFORE: keys = [1, 3, 4, 9], children = [8, 3, 40, 13, 7], count = 5
        // AFTER: keys = [1, 4, 9, 9], children = [8, 3, 13, 7, 7], count = 4

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
    }

    void merge(InnerNode &other, const KeyT &old_separator) {
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
    }

    uint32_t first_sorted() { return *Iterator::begin(Node::count); }
    uint32_t last_sorted() { return 0u; }
    uint32_t next_sorted(uint32_t i) {
        Iterator it(i, Node::count);
        ++it;
        return *it;
    }
    uint32_t prev_sorted(uint32_t i) {
        Iterator it(i, Node::count);
        --it;
        return *it;
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
}
#endif