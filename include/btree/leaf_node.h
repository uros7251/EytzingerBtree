#ifndef INCLUDE_LEAF_NODE_H_
#define INCLUDE_LEAF_NODE_H_

#include "btree/node.h"

namespace guidedresearch {

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize, NodeLayout layout, bool FastInsertion>
struct LeafNode;

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize, NodeLayout layout>
struct LeafNode<KeyT, ValueT, ComparatorT, PageSize, layout, false> : public guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize> {
    
    static_assert(PageSize % CACHELINE == 0); // PageSize should be multiple of cacheline size
    
    using Node = guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize>;
    
    // NOTE: two functions below assume cacheline is 64B!
    static consteval uint32_t GetKeysOffset() {
        // keys array has to be aligned on cacheline ie. 64B
        auto A = std::max(sizeof(Node), alignof(ValueT)); 
        auto block_offset = (PageSize*sizeof(ValueT)+A*sizeof(KeyT))/(64*(sizeof(ValueT)+sizeof(KeyT)));
        if ((64*block_offset-A)/sizeof(ValueT) < (PageSize-64*(block_offset+1))/sizeof(KeyT)) ++block_offset;
        return block_offset*64;
    }
    static consteval uint32_t GetKCapacity() {
        auto A = std::max(sizeof(Node), alignof(ValueT)); 
        auto block_offset = (PageSize*sizeof(ValueT)+A*sizeof(KeyT))/(64*(sizeof(ValueT)+sizeof(KeyT)));
        return std::max((64*block_offset-A)/sizeof(ValueT), (PageSize-64*(block_offset+1))/sizeof(KeyT))-1;
    }

    /// The capacity of a node.
    static constexpr uint32_t kCapacity = GetKCapacity();

    static consteval uint32_t max_keys() { return kCapacity; }
    static consteval uint32_t max_values() { return kCapacity; }

    /// The values. Child associated with the key at index i is stored at position i.
    ValueT values[LeafNode::max_values()+1]; 
    /// The keys. Stored in index range [1, Node::count]. Aligned to cache line size (64 bytes)
    alignas(64) KeyT keys[(PageSize-GetKeysOffset())/sizeof(KeyT)]; // we leave first empty 

    
    using Iterator = typename decltype(IteratorPicker<layout, KeyT, kCapacity+1>())::type;

    /// Constructor.
    LeafNode() : Node(0, 0) { init(); }
    LeafNode(const LeafNode &other) = default;
    
    /// Assignment operator
    LeafNode& operator=(const LeafNode &other) = default;

    /// @brief check if node is less than half full
    /// @return 
    bool merge_needed() const { return Node::count < max_values()/2; } 

    /// Get the index of the first key that is not less than than a provided key.
    /// @param[in] key          The key that should be inserted.
    /// @return                 The index of the first key that is not less than the provided key and boolean indicating if the key is equal.
    std::pair<uint32_t, bool> lower_bound(const KeyT &target) {
        if constexpr (layout == NodeLayout::EYTZINGER) {
            // explained at https://en.algorithmica.org/hpc/data-structures/binary-search/
            if (Node::count == 0) return {0u, false}; // no keys
            
            ComparatorT less;
            uint32_t i=1;
            while (i <= Node::count) {
                __builtin_prefetch(&keys[(CACHELINE/sizeof(KeyT))*i]);
                i = 2*i + less(keys[i], target);
            }

            // recover index
            i >>= std::countr_one(i)+1;

            return std::make_pair(i, keys[i] == target); 
        }
        else if constexpr (layout == NodeLayout::EYTZINGER_SIMD) {
            // explained at https://en.algorithmica.org/hpc/data-structures/s-tree/
            static_assert(std::is_same_v<KeyT, int64_t>); // SIMD support provided only for int64_t (can be extended to other integer types)

            if (Node::count == 0u) return {0u, false}; // no keys

            constexpr auto B = CACHELINE/sizeof(KeyT); // block size
            
            #if defined(__AVX512F__) && defined(__AVX512VL__)
            __m512i key_vec = _mm512_set1_epi64(target);
            #elif defined(__AVX__) && defined(__AVX2__)
            __m256i key_vec = _mm256_set1_epi64x(target);
            #else
            ComparatorT less;
            #endif

            uint32_t N = Node::count/B+1u; // block_count = block index of the last element + 1
            // save last not-less-than target
            // ComparatorT less;
            auto k=0u;
            for (uint32_t i=0, j=1, mask=1; i < N; i+=i*B+j, j=0, mask=0) {
                assert(reinterpret_cast<uintptr_t>(&keys[i*B]) % CACHELINE == 0);
                // comparison
                #if defined(__AVX512F__) && defined(__AVX512VL__)
                __m512i y_vec = _mm512_load_si512(&keys[i*B]);
                mask = _mm512_cmplt_epi64_mask(y_vec, key_vec); // we assume target > kNegInf, replace '=' with '|=' for any target
                j = std::countr_one(mask);
                #elif defined(__AVX__) && defined(__AVX2__)
                mask = less256(key_vec, &keys[i*B]) + (less256(key_vec, &keys[i*B+4]) << 4);
                j = std::countr_one(mask);
                #else
                for (; j<B && less(keys[i*B+j], target); ++j);
                #endif
                k = j < B ? i*B+j : k; // save descent to the left
            }

            return std::make_pair(k, keys[k] == target);
        }
        else if constexpr (layout == NodeLayout::ORDERED) {
            // explained at https://en.algorithmica.org/hpc/data-structures/binary-search/
            if (Node::count==0) return {0u, false}; // no keys
        
            ComparatorT less;
            uint32_t i=1, n = Node::count;
            // branchless binary search
            while (n>1) {
                auto half = n/2;
                n -= half; // ceil(n/2)
                __builtin_prefetch(&keys[i+n/2-1]); // prefetch left
                __builtin_prefetch(&keys[i+half+n/2-1]); // prefetch right
                i += less(keys[i+half-1], target) * half; // hopefully compiler translates this to cmov
            }

            return (i==Node::count && less(keys[i], target)) ? // check if last key is less than key
                std::make_pair(0u, false) : // return index one after last key
                std::make_pair(i, keys[i] == target);   
        }
        else __builtin_unreachable();
    }

    /// Insert a key.
    /// @param[in] key          The key that should be inserted.
    /// @param[in] value        The value that should be inserted.
    void insert(const KeyT &key, const ValueT &value) {
        // for example, insert(3, 40)
        // BEFORE: in_order(keys) = [1, 4, 9], values = [8, 3, 13, 7], count = 3
        // AFTER: in_order(keys) = [1, 3, 4, 9], values = [8, 4o, 13, 7], count = 4

        assert(Node::count < LeafNode::kCapacity); // more place available

        // trivial case with no keys
        if (Node::count == 0u) {
            // no keys
            keys[1] = key;
            values[1] = value;
            ++Node::count;
            return;
        }

        auto [index,found] = lower_bound(key);
        if (found) {
            values[index] = value;
            return;
        }

        // the idea is to logically insert a (key, child) pair into the node and then move it to the right place.
        // we start at the index Node::count+1 (which is the next empty slot) and first determine if to-be-inserted
        // key should be placed to the left or to the right of it in ascending order (note: if i<j, it might be that 
        // in_order(i) > in_order(j))
        
        ++Node::count;
        Iterator lead(Node::count, Node::count+1);
        Iterator lag = lead;
        ComparatorT less;
        bool forward;
        // determine if we should go forward or backward in in-order traversal
        if (!index || lead == Iterator::begin(Node::count+1)) forward = true;
        else if (lead == Iterator::rbegin(Node::count+1)) forward = false;
        else {
            ++lead;
            if (*lead == index) goto end; // we insert key into new slot at index Node::count
            if (less(keys[*lead], key)) {
                forward = true;
                keys[*lag] = keys[*lead];
                values[*lag] = values[*lead];
                lag = lead;
            }
            else {
                forward = false;
                lead = lag;
            }
        }
        assert(lag == lead);
        if (forward) {
            ++lead;
            for (; *lead != index; lag=lead, ++lead) {
                keys[*lag] = keys[*lead];
                values[*lag] = values[*lead];
            }
        }
        else {
            --lead;
            for (; *lag != index; lag=lead, --lead) {
                keys[*lag] = keys[*lead];
                values[*lag] = values[*lead];
            }
        }
        end:
        keys[*lag] = key;
        values[*lag] = value;
    }

    /// Erase a key.
    bool erase(const KeyT &key) {
        // for example, erase(3)
        // BEFORE: keys = [1, 3, 4, 9], values = [8, 3, 40, 13], count = 4
        // AFTER: keys = [1, 4, 9], values = [8, 40, 13], count = 3

        auto [index, found] = lower_bound(key);
        if (!found) return false; // key not found

        Iterator lag(index, Node::count+1);
        Iterator lead = lag;
        ComparatorT less;
        
        if (less(key, keys[Node::count])) {
            // go right
            ++lead;
            for (; *lag != Node::count; lag=lead, ++lead) {
                keys[*lag] = keys[*lead];
                values[*lag] = values[*lead];
            }
        }
        else {
            // go left
            --lead;
            for (; *lag != Node::count; lag=lead, --lead) {
                keys[*lag] = keys[*lead];
                values[*lag] = values[*lead];
            }
        }
        keys[Node::count] = kNegInf;
        --Node::count;
        return true;
    }

    /// Split the node.
    /// @param[in] buffer       The buffer for the new page.
    /// @return                 The separator key.
    KeyT split(char* buffer) {
        LeafNode* new_node = new (buffer) LeafNode();
        // new_node->count = Node::count;
        // old_node retains ceil(count/2) keys
        // new_node gets floor(count/2) keys
        uint16_t left_node_count = (Node::count+1u)/2, right_node_count = Node::count/2;
        assert(left_node_count+right_node_count==Node::count);
        
        Iterator old_it = Iterator::rbegin(Node::count+1u), left_it = Iterator::rbegin(left_node_count+1u), right_it = Iterator::rbegin(right_node_count+1u);
        for (auto i=0; i<right_node_count; ++i, --old_it, --right_it) {
            new_node->keys[*right_it] = keys[*old_it];
            new_node->values[*right_it] = values[*old_it];
        }
        assert(right_it == Iterator::rend(right_node_count+1));
        
        auto separator = keys[*old_it];
        for (auto i=0; i<left_node_count; ++i, --old_it, --left_it) {
            keys[*left_it] = keys[*old_it];
            values[*left_it] = values[*old_it];
        }
        assert(left_it == Iterator::rend(left_node_count+1u));
        assert(old_it == Iterator::rend(Node::count+1u));
        for (auto i=left_node_count+1; i<=Node::count; ++i) {
            keys[i] = kNegInf;
        }
        Node::count = left_node_count;
        new_node->count = right_node_count;
        return separator;
    }

    void merge(LeafNode &other) {
        Iterator left_it = Iterator::begin(Node::count+1u), right_it = Iterator::begin(other.count+1u), merge_it = Iterator::begin(Node::count+other.count+1u);
        // rearrange values in left node
        for (auto i=0; i<Node::count; ++i, ++left_it, ++merge_it) {
            keys[*merge_it] = keys[*left_it];
            values[*merge_it] = values[*left_it];
        }
        assert(left_it == Iterator::end(Node::count+1u));
        // insert values from right node
        for (auto i=0; i<other.count; ++i, ++right_it, ++merge_it) {
            keys[*merge_it] = other.keys[*right_it];
            values[*merge_it] = other.values[*right_it];
        }
        assert(right_it == Iterator::end(other.count+1u));
        // update counts
        Node::count += other.count;
        other.count = 0;
    }

    template<typename Function>
    void for_each(Function &f, KeyT lower_bound = kNegInf, KeyT upper_bound = kInf) {
        // if possible, circumvent the iterator and directly access the keys and values
        if (lower_bound == kNegInf && upper_bound == kInf) {
            for (auto i=1u; i<=Node::count; ++i) {
                f(keys[i], values[i]);
            }
            return;
        }

        Iterator it = 
            lower_bound == kNegInf ? Iterator::begin(Node::count+1) : // needed because lower_bound implementation assumes all keys are > kNegInf
            Iterator(this->lower_bound(lower_bound).first, Node::count+1);
        auto [end, include_end] = upper_bound == kInf ? std::make_pair(0u, false) : this->lower_bound(upper_bound);
        for (; *it != end; ++it) {
            f(keys[*it], values[*it]);
        }
        if (include_end) {
            f(keys[end], values[end]);
        }
    }

    /// Returns the keys.
    std::vector<KeyT> get_key_vector() {
        std::vector<KeyT> sorted_keys;
        sorted_keys.reserve(Node::count);
        for (auto it = Iterator::begin(Node::count+1u), end = Iterator::end(Node::count+1u); it != end; ++it) {
            sorted_keys.emplace_back(keys[*it]);
        }
        return sorted_keys;
    }

    /// Returns the values.
    std::vector<ValueT> get_value_vector() {
        std::vector<ValueT> sorted_values;
        sorted_values.reserve(Node::count);
        for (auto it = Iterator::begin(Node::count+1u), end = Iterator::end(Node::count+1u); it != end; ++it) {
            sorted_values.emplace_back(values[*it]);
        }
        return sorted_values;
    }

    static void rebalance(LeafNode &left, LeafNode &right, KeyT &separator) {
        if (left.count > right.count) {
            // shift right
            uint16_t to_shift = (left.count-right.count)/2; // shift half the difference
            Iterator old_right_it = Iterator::rbegin(right.count+1), new_right_it = Iterator::rbegin(right.count+to_shift+1),
                        old_left_it = Iterator::rbegin(left.count+1), new_left_it = Iterator::rbegin(left.count-to_shift+1);
            // make space for keys and values from left node
            for (int i=0; i<right.count; ++i, --old_right_it, --new_right_it) {
                right.keys[*new_right_it] = right.keys[*old_right_it];
                right.values[*new_right_it] = right.values[*old_right_it];
            }
            assert(old_right_it == Iterator::rend(right.count+1));
            // insert keys and values from left node
            for (int i=0; i<to_shift; ++i, --old_left_it, --new_right_it) {
                right.keys[*new_right_it] = left.keys[*old_left_it];
                right.values[*new_right_it] = left.values[*old_left_it];
            }
            assert(new_right_it == Iterator::rend(right.count+to_shift+1));
            // update separator
            separator = left.keys[*old_left_it];
            // rearrange keys and values in left node
            for (int i=0; i<left.count-to_shift; ++i, --old_left_it, --new_left_it) {
                left.keys[*new_left_it] = left.keys[*old_left_it];
                left.values[*new_left_it] = left.values[*old_left_it];
            }
            for (unsigned i=left.count-to_shift+1; i<=left.count; ++i) {
                left.keys[i] = kNegInf;
            }
            assert(old_left_it == Iterator::rend(left.count+1));
            assert(new_left_it == Iterator::rend(left.count-to_shift+1));
            // update counts
            left.count -= to_shift;
            right.count += to_shift;
        }
        else {
            // shift left
            uint16_t to_shift = (right.count-left.count)/2; // shift half the difference
            Iterator old_right_it = Iterator::begin(right.count+1), new_right_it = Iterator::begin(right.count-to_shift+1),
                        old_left_it = Iterator::begin(left.count+1), new_left_it = Iterator::begin(left.count+to_shift+1);
            // make space for keys and values from left node
            for (int i=0; i<left.count; ++i, ++old_left_it, ++new_left_it) {
                left.keys[*new_left_it] = left.keys[*old_left_it];
                left.values[*new_left_it] = left.values[*old_left_it];
            }
            assert(old_left_it == Iterator::end(left.count+1));
            // insert keys and values from right node
            for (int i=0; i<to_shift-1; ++i, ++old_right_it, ++new_left_it) {
                left.keys[*new_left_it] = right.keys[*old_right_it];
                left.values[*new_left_it] = right.values[*old_right_it];
            }
            assert(new_left_it == Iterator::rbegin(left.count+to_shift+1));
            // we take out the last iteration out of loop to update separator as well
            left.keys[*new_left_it] = right.keys[*old_right_it];
            left.values[*new_left_it] = right.values[*old_right_it];
            separator = right.keys[*old_right_it];
            ++old_right_it;
            // rearrange keys and values in right node
            for (int i=0; i<right.count-to_shift; ++i, ++old_right_it, ++new_right_it) {
                right.keys[*new_right_it] = right.keys[*old_right_it];
                right.values[*new_right_it] = right.values[*old_right_it];
            }
            for (unsigned i=right.count-to_shift+1; i<=right.count; ++i) {
                right.keys[i] = kNegInf;
            }
            assert(old_right_it == Iterator::end(right.count+1));
            assert(new_right_it == Iterator::end(right.count-to_shift+1));
            // update counts
            left.count += to_shift;
            right.count -= to_shift;
        }
    }
    
private:
    private:
    void init() {
        assert((sizeof(keys)/sizeof(KeyT))%8==0);
        for (auto i=0u; i<sizeof(keys)/sizeof(KeyT); ++i) {
            keys[i] = kNegInf;
        }
    }

    #if !defined(__AVX512F__) && defined(__AVX__) && defined(__AVX2__)
    static int less256(__m256i x_vec, KeyT* y_ptr) {
        __m256i y_vec = _mm256_load_si256((__m256i*) y_ptr); // load 4 sorted elements
        __m256i mask = _mm256_cmpgt_epi64(x_vec, y_vec); // compare against the key
        return _mm256_movemask_pd((__m256d) mask);    // extract the 4-bit mask
    }
    #endif
};

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
struct LeafNode<KeyT, ValueT, ComparatorT, PageSize, NodeLayout::SORTED, false> : public guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize> {

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
    /// @param[in] target       The key to look for.
    /// @return                 The index of the first key that is not less than the provided key and boolean indicating if the key is equal.
    std::pair<uint32_t, bool> lower_bound(const KeyT &target) {
        if (Node::count==0) return {0u, false}; // no keys
        
        ComparatorT less;

        uint32_t i=0, n=Node::count;
        // branchless binary search
        while (n>1) {
            auto half = n/2;
            n -= half; // ceil(n/2)
            // __builtin_prefetch(&keys[i+n/2-1]); // prefetch left
            // __builtin_prefetch(&keys[i+half+n/2-1]); // prefetch right
            i += less(keys[i+half-1], target) * half; // hopefully compiler translates this to cmov
        }

        return (i==Node::count-1u && less(keys[i], target)) ? // check if last key is less than key
            std::make_pair(static_cast<uint32_t>(Node::count), false) : // return index one after last key
            std::make_pair(i, keys[i] == target);
    }

    /// Insert a key.
    /// @param[in] key          The key that should be inserted.
    /// @param[in] value        The value that should be inserted.
    void insert(const KeyT &key, const ValueT &value) {
        assert(Node::count < kCapacity);
        auto [index, found] = lower_bound(key);
        if (found) {
            values[index] = value;
            return;
        }

        for (auto i=Node::count; i>index; --i) {
            keys[i] = keys[i-1];
            values[i] = values[i-1];
        }
        keys[index] = key;
        values[index] = value;
        ++Node::count;
    }

    /// Erase a key.
    bool erase(const KeyT &key) {
        // try to find key
        auto [index, found] = lower_bound(key);
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

    template<typename Function>
    void for_each(Function &f, KeyT lower_bound = kNegInf, KeyT upper_bound = kInf) {
        // circumvent calls to lower_bound if possible
        auto start = lower_bound == kNegInf ? 0u : this->lower_bound(lower_bound).first;
        auto [end, include_end] = upper_bound == kInf ? std::make_pair(static_cast<uint32_t>(Node::count), false) : this->lower_bound(upper_bound);
        end = include_end ? end+1u : end;
        for (uint16_t i=start; i<end; ++i) {
            f(keys[i], values[i]);
        }
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
struct LeafNode<KeyT, ValueT, ComparatorT, PageSize, NodeLayout::SORTED, true> : public guidedresearch::Node<KeyT, ValueT, ComparatorT, PageSize> {

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
    /// @tparam indirect        If false (default), return the index of the key in the keys array, otherwise return the index i such that keys[indices[i]] == key.
    /// @return                 The index of the first key that is not less than the provided key and boolean indicating if the key is equal.
    template<bool indirect = false>
    std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
        if (Node::is_leaf() && Node::count==0) return {0u, false}; // no keys
        
        ComparatorT less;
        uint32_t i=0, n = Node::count;
        // branchless binary search
        while (n>1) {
            auto half = n/2;
            n -= half; // ceil(n/2)
            // __builtin_prefetch(&indices[i+n/2-1]); // prefetch left
            // __builtin_prefetch(&indices[i+half+n/2-1]); // prefetch right
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
        auto [index, found] = lower_bound<true>(key);
        if (found) {
            values[indices[index]] = value;
            return;
        }
        keys[Node::count] = key;
        values[Node::count] = value;
        for (auto i=Node::count; i>index; --i) {
            indices[i] = indices[i-1];
        }
        indices[index] = Node::count;
        ++Node::count;
    }

    /// Erase a key.
    bool erase(const KeyT &key) {
        // try to find key
        auto [index, found] = lower_bound<true>(key);
        if (!found) return false; // key not found
        // find i such that indices[i] == Node::count-1
        auto [index_last, found_last] = lower_bound<true>(keys[Node::count-1]);
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

    template<typename Function>
    void for_each(Function &f, KeyT lower_bound = kNegInf, KeyT upper_bound = kInf) {
        // circumvent calls to lower_bound if possible
        auto start = lower_bound == kNegInf ? 0u : this->lower_bound<true>(lower_bound).first;
        auto [end, include_end] = upper_bound == kInf ? std::make_pair(static_cast<uint32_t>(Node::count), false) : this->lower_bound<true>(upper_bound);
        end = include_end ? end+1u : end;
        for (uint16_t i=start; i<end; ++i) {
            f(keys[indices[i]], values[indices[i]]);
        }
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