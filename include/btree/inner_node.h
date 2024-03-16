#ifndef INCLUDE_INNER_NODE_H
#define INCLUDE_INNER_NODE_H

#include <bit>
#include <buffer_manager/swip.h>
#include <btree/node.h>
#include <immintrin.h>
#include <avx2intrin.h>

namespace guidedresearch {

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize, NodeLayout layout>
struct InnerNode: public Node<KeyT, ValueT, ComparatorT, PageSize> {
    
    static_assert(PageSize % CACHELINE == 0); // PageSize should be multiple of cacheline size
    
    using Node = guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize>;
    using Iterator = typename decltype(IteratorPicker<layout, KeyT>())::type;

    // NOTE: two functions below assume cacheline is 64B!
    static consteval uint32_t GetKeysOffset() {
        // keys array has to be aligned on cacheline ie. 64B
        auto block_offset = (PageSize+sizeof(KeyT))/(64+8*sizeof(KeyT));
        if (8*block_offset-1 < (PageSize-64*(block_offset+1))/sizeof(KeyT)) ++block_offset;
        return block_offset*64;
    }
    static consteval uint32_t GetKCapacity() {
        auto block_offset = (PageSize+sizeof(KeyT))/(64+8*sizeof(KeyT));
        return std::max(8*block_offset-1, (PageSize-64*(block_offset+1))/sizeof(KeyT))-1;
    }

    /// The capacity of a node.
    static constexpr uint32_t kCapacity = GetKCapacity();

    static consteval uint32_t max_keys() { return kCapacity; }
    static consteval uint32_t max_children() { return kCapacity+1; }

    /// The children. Child associated with the key at index i is stored at position i. The remaining child is stored at the 0th position.
    Swip children[InnerNode::max_children()]; 
    /// The keys. Stored starting at index 1. Aligned to cache line size (64 bytes)
    alignas(64) KeyT keys[(PageSize-GetKeysOffset())/sizeof(KeyT)]; // we leave first empty 

    /// Constructor.
    InnerNode() : Node(0, 0) { init(); }
    InnerNode(uint16_t level) : Node(level, 0) { init(); }
    InnerNode(const InnerNode &other) = default;

    InnerNode& operator=(const InnerNode& other) = default;

    /// this method can be virtual, but we avoid that on purpose
    /// @brief check if node is less than half full
    /// @return 
    bool merge_needed() const { return Node::count < max_children()/2; } 

    std::pair<uint32_t, bool> lower_bound(const KeyT &target) {
        if constexpr (layout == NodeLayout::EYTZINGER) {
            // explained at https://en.algorithmica.org/hpc/data-structures/binary-search/
            if (Node::count <= 1) return {0u, false}; // no keys
            
            ComparatorT less;
            uint32_t i=1;
            while (i < Node::count) {
                //__builtin_prefetch(&keys[(CACHELINE/sizeof(KeyT))*i]);
                i = 2*i + less(keys[i], target);
            }

            // recover index
            i >>= std::countr_one(i)+1;

            return i ? // check if last key is less than key
                std::make_pair(i, keys[i] == target) :
                std::make_pair(0u, false); 
        }
        else if constexpr (layout == NodeLayout::EYTZINGER_SIMD) {
            // explained at https://en.algorithmica.org/hpc/data-structures/s-tree/
            static_assert(std::is_same_v<KeyT, int32_t>); // SIMD support provided only for int32_t

            if (Node::count <= 1) return {0u, false}; // no keys

            constexpr auto B = CACHELINE/sizeof(KeyT); // block size
            
            #if defined(__AVX512F__) && defined(__AVX512VL__)
            __m512i key_vec = _mm512_set1_epi32(target);
            #elif defined(__AVX__) && defined(__AVX2__)
            __m256i key_vec = _mm256_set1_epi32(target);
            #else
            ComparatorT less;
            #endif

            uint32_t N = (Node::count-1u)/B+1; // block_count = ceil(count/B)
            // save last not-less-than target
            // ComparatorT less;
            auto k=0u;
            for (uint32_t i=0, j=1, mask=1; i < N; i+=i*B+j, j=0, mask=0) {
                // comparison
                #if defined(__AVX512F__) && defined(__AVX512VL__)
                __m512i y_vec = _mm512_load_si512(&keys[i*B]);
                mask |= _mm512_cmplt_epi32_mask(y_vec, key_vec);
                #elif defined(__AVX__) && defined(__AVX2__)
                mask |= less256(key_vec, &keys[i*B]) + (less256(key_vec, &keys[i*B+8]) << 8);
                #else
                for (; j<B; ++j) {
                    mask |= (less(keys[i*B+j], target) << j);
                }
                #endif
                j = std::countr_one(mask);
                if (j < B) {
                    k = i*B+j;
                }
            }

            if (k >= Node::count) k/=B;

            return k ? // check if last key is less than key
                std::make_pair(k, keys[k] == target) :
                std::make_pair(0u, false); 
        }
        else if constexpr (layout == NodeLayout::SORTED) {
            // explained at https://en.algorithmica.org/hpc/data-structures/binary-search/
            if (Node::count <= 1) return {0u, false}; // no keys
            
            ComparatorT less;
            uint32_t i=1, n = Node::count-1u;
            // branchless binary search
            while (n>1) {
                auto half = n/2;
                n -= half; // ceil(n/2)
                i += less(keys[i+half-1], target) * half; // hopefully compiler translates this to cmov
            }

            return (i==Node::count-1u) && less(keys[i], target) ? // check if last key is less than key
                std::make_pair(0u, false) : // return index one after last key
                std::make_pair(i, keys[i] == target);   
        }
        else return std::make_pair(0u, false);
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
        // for example, erase(3)
        // BEFORE: keys = [1, 3, 4, 9], children = [8, 3, 40, 13, 7], count = 5
        // AFTER: keys = [1, 4, 9, 9], children = [8, 3, 13, 7, 7], count = 4

        Iterator lag(index, Node::count);
        Iterator lead = lag;
        ComparatorT less;
        
        // prepare children array for later shifting
        if (lead == Iterator::rbegin(Node::count)) {
            // if index comes last during in-order traversal
            children[0] = children[index];
            // it's guaranteed that we go left after (take else branch)
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
        if (Node::count) keys[Node::count] = std::numeric_limits<KeyT>::max();
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
            keys[*old_it] = std::numeric_limits<KeyT>::max();
        }
        assert(right_it == Iterator::rend(right_node_count));
        
        auto separator = keys[*old_it];
        keys[*old_it] = std::numeric_limits<KeyT>::max();
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

    std::vector<uint64_t> get_pid_vector() {
        std::vector<uint64_t> sorted_pids;
        sorted_pids.reserve(Node::count);
        for (auto it = Iterator::begin(Node::count), end = Iterator::end(Node::count); it != end; ++it) {
            sorted_pids.push_back(children[*it].asPageID());
        }
        sorted_pids.push_back(children[0].asPageID());
        return sorted_pids;
    }

    /// Returns the keys in sorted order
    std::vector<KeyT> get_key_vector() {
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
            for (unsigned i=left.count-to_shift; i<sizeof(left.keys)/sizeof(KeyT); ++i) {
                left.keys[i] = std::numeric_limits<KeyT>::max();
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
            for (unsigned i=right.count-to_shift; i<sizeof(right.keys)/sizeof(KeyT); ++i) {
                right.keys[i] = std::numeric_limits<KeyT>::max();
            }
            assert(old_right_it == Iterator::end(right.count));
            assert(new_right_it == Iterator::end(right.count-to_shift));
            // update counts
            left.count += to_shift;
            right.count -= to_shift;
        }
    }

    private:
    void init() {
        assert((sizeof(keys)/sizeof(KeyT))%8==0);
        // keys[0] = std::numeric_limits<KeyT>::min();
        for (auto i=1u; i<sizeof(keys)/sizeof(KeyT); ++i) {
            keys[i] = std::numeric_limits<KeyT>::max();
        }
    }

    #if !defined(__AVX512F__) && defined(__AVX__) && defined(__AVX2__)
    static int less256(__m256i x_vec, int* y_ptr) {
        __m256i y_vec = _mm256_load_si256((__m256i*) y_ptr); // load 8 sorted elements
        __m256i mask = _mm256_cmpgt_epi32(x_vec, y_vec); // compare against the key
        return _mm256_movemask_ps((__m256) mask);    // extract the 8-bit mask
    }
    #endif
};


template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct InnerNode<KeyT, ValueT, ComparatorT, PageSize, NodeLayout::SORTED> : public Node<KeyT, ValueT, ComparatorT, PageSize> {
    using Node = guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize>;
    /// Remark: it might be better to store the rightmost key to make nodes more uniform and algorithms more concise

    /// The capacity of a node.
    static constexpr uint32_t kCapacity = (PageSize-sizeof(Node)-sizeof(Swip))/(sizeof(KeyT)+sizeof(Swip));

    inline static consteval uint32_t max_keys() { return kCapacity; }
    inline static consteval uint32_t max_children() { return kCapacity+1; }

    /// The keys.
    KeyT keys[InnerNode::max_keys()]; 
    /// The values.
    Swip children[InnerNode::max_children()]; 

    /// Constructor.
    InnerNode() : Node(0, 0) { }
    InnerNode(uint16_t level) : Node(level, 0) {}
    InnerNode(const InnerNode &other) : Node(other) {
        copy();
    }
    
    /// Assignment operator
    InnerNode& operator=(const InnerNode &other) {
        copy(other);
        return *this;
    }

    /// @brief check if node is less than half full (this method can be virtual, but we avoid that on purpose)
    /// @return 
    bool merge_needed() const { return Node::count < max_children()/2; } 

    /// Get the index of the first key that is not less than than a provided key.
    /// @param[in] key          The key that should be inserted.
    /// @return                 The index of the first key that is not less than the provided key and boolean indicating if the key is equal.
    std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
        if (Node::count <= 1) return {0u, false}; // no keys
        
        ComparatorT less;
        uint32_t i=0, n = Node::count-1u;
        // branchless binary search
        while (n>1) {
            auto half = n/2;
            n -= half; // ceil(n/2)
            // __builtin_prefetch(&keys[i+n/2-1]); // prefetch left
            // __builtin_prefetch(&keys[i+half+n/2-1]); // prefetch right
            i += less(keys[i+half-1], key) * half; // hopefully compiler translates this to cmov
        }

        return (i==Node::count-2u && less(keys[i], key)) ? // check if last key is less than key
            std::make_pair(i+1, false) : // return index one after last key
            std::make_pair(i, keys[i] == key); 
    }

    /// Insert a key.
    /// @param[in] key          The key that should be inserted.
    /// @param[in] child        The child that should be inserted.
    void insert_split(const KeyT &key, Swip child) {
        // for example, insert_split(3, 40)
        // BEFORE: keys = [1, 4, 9, ?], children = [8, 3, 13, 7, ?], count = 4
        // AFTER: keys = [1, 3, 4, 9], children = [8, 3, 40, 13, 7], count = 5
        
        assert(Node::count > 0);
        assert(Node::count-1 != InnerNode::kCapacity);
        uint32_t i; // current contested position (keys[i-1] and key are competing for ith position)
        for (i=Node::count-1; i > 0 && keys[i-1]>key; --i) {
            keys[i] = keys[i-1];
            children[i+1] = children[i];
        }
        keys[i] = key; // insert new separator
        children[i+1] = child; // insert child
        ++Node::count;
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
        for (unsigned i=0u, pos=(Node::count+1)/2; i<Node::count/2-1u; ++i) {
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

    std::vector<uint64_t> get_pid_vector() {
        std::vector<uint64_t> pids;
        pids.reserve(Node::count);
        for (uint32_t i=0; i<Node::count; ++i) {
            pids.push_back(children[i].asPageID());
        }
        return pids;
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

}
#endif