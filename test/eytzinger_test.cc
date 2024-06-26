#include <gtest/gtest.h>
#include <random>

#include "btree/btree.h"

#define PageSize 4096

using BufferFrame = guidedresearch::BufferFrame;
using BufferManager = guidedresearch::BufferManager;
using AlignedVector = guidedresearch::AlignedVector;
using KeyT = int64_t;
using ValueT = int64_t;
using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::EYTZINGER_SIMD>;
using Swip = guidedresearch::Swip;

namespace {

TEST(EytzingerTest, StaticBlockIterator) {
    using Iterator = guidedresearch::StaticBlockIterator<64, 64/sizeof(KeyT)>;

    auto it = Iterator::begin(20);
    ASSERT_EQ(*it, 8);
    std::vector<uint32_t> expected = {8,9,10,11,12,13,14,15,1,16,17,18,19,2,3,4,5,6,7};
    for (auto i = 0u; i < expected.size(); ++i) {
        ASSERT_EQ(*it, expected[i]);
        ++it;
    }
    ASSERT_TRUE(it == Iterator::end(20));

    it = Iterator::rbegin(20);
    ASSERT_EQ(*it, 7);
    for (auto i = expected.size(); i > 0;) {
        ASSERT_EQ(*it, expected[--i]);
        --it;
    }
    ASSERT_TRUE(it == Iterator::rend(20));

    auto it2 = Iterator::begin(1);
    ASSERT_TRUE(it2 == Iterator::end(1));

    it2 = Iterator::rbegin(1);
    ASSERT_TRUE(it2 == Iterator::rend(1));
}

TEST(EytzingerTest, BlockIterator) {
    using Iterator = guidedresearch::EytzingerIterator<64/sizeof(KeyT)>;
    auto it = Iterator::begin(20);
    ASSERT_EQ(*it, 8);
    std::vector<uint32_t> expected = {8,9,10,11,12,13,14,15,1,16,17,18,19,2,3,4,5,6,7};
    for (auto i = 0u; i < expected.size(); ++i) {
        ASSERT_EQ(*it, expected[i]);
        ++it;
    }
    ASSERT_TRUE(it == Iterator::end(20));

    it = Iterator::rbegin(20);
    ASSERT_EQ(*it, 7);
    for (auto i = expected.size(); i > 0;) {
        ASSERT_EQ(*it, expected[--i]);
        --it;
    }
    ASSERT_TRUE(it == Iterator::rend(20));

    auto it2 = Iterator::begin(1);
    ASSERT_TRUE(it2 == Iterator::end(1));

    it2 = Iterator::rbegin(1);
    ASSERT_TRUE(it2 == Iterator::rend(1));
}

TEST(EytzingerTest, InOrderIterator) {
    using Iterator = guidedresearch::EytzingerIterator<1u>;
    auto it = Iterator::begin(12);
    ASSERT_EQ(*it, 8);
    std::vector<uint32_t> expected = {8, 4, 9, 2, 10, 5, 11, 1, 6, 3, 7};
    for (auto i = 0u; i < expected.size(); ++i) {
        ASSERT_EQ(*it, expected[i]);
        ++it;
    }
    ASSERT_TRUE(it == Iterator::end(12));

    it = Iterator::begin(16);
    ASSERT_EQ(*it, 8);

    it = Iterator::rbegin(12);
    ASSERT_EQ(*it, 7);
    for (auto i = expected.size(); i > 0;) {
        ASSERT_EQ(*it, expected[--i]);
        --it;
    }
    ASSERT_TRUE(it == Iterator::rend(12));

    auto it2 = Iterator::begin(1);
    ASSERT_TRUE(it2 == Iterator::end(1));

    it2 = Iterator::rbegin(1);
    ASSERT_TRUE(it2 == Iterator::rend(1));
}

TEST(EytzingerTest, InnerNodeInsert) {
    AlignedVector buffer;
    buffer.resize(PageSize, 1024);
    // init inner node
    auto node = new (buffer.data()) BTree::InnerNode();
    node->children[0] = Swip::fromPID(0);
    node->count = 1;

    auto n = BTree::InnerNode::kCapacity;

    // Insert children into the leaf nodes
    for (uint32_t i = 1, j = 2; i <= n; ++i, j = i * 2) {
        Swip s = Swip::fromPID(j);
        node->insert_split(i, s);

        ASSERT_EQ(node->count, i + 1)
            << "LeafNode::insert did not increase the the child count";
    }

    // Check the number of keys & children
    auto keys = node->get_key_vector();
    auto pids = node->get_pid_vector();
    for (uint32_t i = 0, len = keys.size(); i < len; ++i) {
        ASSERT_EQ(keys[i], i+1);
        ASSERT_EQ(2*keys[i], pids[i+1]);
    }

    ASSERT_EQ(keys.size(), n)
        << "leaf node must contain " << n << " keys for " << n+1 << " values";
}

TEST(EytzingerTest, InnerNodeInsertDecreasing) {
    AlignedVector buffer;
    buffer.resize(PageSize, 1024);
    // init inner node
    auto node = new (buffer.data()) BTree::InnerNode();
    node->children[0] = Swip::fromPID(0);
    node->count = 1;

    auto n = BTree::InnerNode::kCapacity;

    // Insert children into the leaf nodes
    for (uint32_t i = n, j = n*2; i > 0; --i, j = i * 2) {
        Swip s = Swip::fromPID(j);
        node->insert_split(i, s);

        ASSERT_EQ(node->count, n - i + 2)
            << "LeafNode::insert did not increase the the child count";
    }

    // Check the number of keys & children
    auto keys = node->get_key_vector();
    auto pids = node->get_pid_vector();
    for (uint32_t i = 0, len = keys.size(); i < len; ++i) {
        ASSERT_EQ(keys[i], i+1) << "i=" << i;
        ASSERT_EQ(2*keys[i], pids[i+1]);
    }

    ASSERT_EQ(keys.size(), n)
        << "leaf node must contain " << n << " keys for " << n+1 << " values";
}

TEST(EytzingerTest, InnerNodeErase) {
    AlignedVector buffer;
    buffer.resize(PageSize, 1024);
    // init inner node
    auto node = new (buffer.data()) BTree::InnerNode();
    node->children[0] = Swip::fromPID(0);
    node->count = 1;

    auto n = BTree::InnerNode::kCapacity;

    // Insert children into the leaf nodes
    for (uint32_t i = 1, j = 2; i <= n; ++i, j = i * 2) {
        Swip s = Swip::fromPID(j);
        node->insert_split(i, s);
    }

    for (uint32_t i = 1; i <= n; ++i) {
        auto [index, found] = node->lower_bound(i);
        ASSERT_TRUE(found)
            << "k=" << i << " was not in the tree";
        node->erase(index);
        std::tie(index, found) = node->lower_bound(i);
        ASSERT_FALSE(found)
            << "k=" << i << " was not removed from the tree";
    }
}

TEST(EytzingerTest, InnerNodeEraseDecreasing) {
    AlignedVector buffer;
    buffer.resize(PageSize, 1024);
    // init inner node
    auto node = new (buffer.data()) BTree::InnerNode();
    node->children[0] = Swip::fromPID(0);
    node->count = 1;

    auto n = BTree::InnerNode::kCapacity;

    // Insert children into the leaf nodes
    for (uint32_t i = 1, j = 2; i <= n; ++i, j = i * 2) {
        Swip s = Swip::fromPID(j);
        node->insert_split(i, s);
    }

    for (uint32_t i = n; i > 0; --i) {
        auto [index, found] = node->lower_bound(i);
        ASSERT_TRUE(found)
            << "k=" << i << " was not in the tree";
        node->erase(index);
        std::tie(index, found) = node->lower_bound(i);
        ASSERT_FALSE(found)
            << "k=" << i << " was not removed from the tree";
    }
}

TEST(EytzingerTest, InnerNodeEraseRandomPermutation) {
    AlignedVector buffer;
    buffer.resize(PageSize, 1024);
    // init inner node
    auto node = new (buffer.data()) BTree::InnerNode();
    node->children[0] = Swip::fromPID(0);
    node->count = 1;

    auto n = BTree::InnerNode::kCapacity;

    std::vector<uint32_t> keys;
    keys.reserve(n);

    // Insert children into the leaf nodes
    for (uint32_t i = 1, j = 2; i <= n; ++i, j = i * 2) {
        Swip s = Swip::fromPID(j);
        node->insert_split(i, s);
        keys.push_back(i);
    }

    // std::random_device rd;
    std::mt19937 g(7251);
    std::shuffle(keys.begin(), keys.end(), g);

    for (auto i : keys) {
        auto [index, found] = node->lower_bound(i);
        ASSERT_TRUE(found)
            << "k=" << i << " was not in the tree";
        node->erase(index);
        std::tie(index, found) = node->lower_bound(i);
        ASSERT_FALSE(found)
            << "k=" << i << " was not removed from the tree";
    }
    
}

TEST(EytzingerTest, InnerNodeSplit) {
    AlignedVector buffer_left;
    AlignedVector buffer_right;
    buffer_left.resize(PageSize, 1024);
    buffer_right.resize(PageSize, 1024);

    auto left_node = new (buffer_left.data()) BTree::InnerNode();
    left_node->children[0] = Swip::fromPID(0u);
    left_node->count = 1;
    auto right_node = reinterpret_cast<BTree::InnerNode*>(buffer_right.data());

    auto n = BTree::InnerNode::kCapacity;

    // Fill the left node
    for (uint32_t i = 1, j = 2; i <= n; ++i, j = i * 2) {
        left_node->insert_split(i, j);
    }

    // Check the number of keys & children
    auto left_keys = left_node->get_key_vector();
    auto left_values = left_node->get_pid_vector();
    ASSERT_EQ(left_keys.size(), n);
    ASSERT_EQ(left_values.size(), n+1);

    // Now split the left node
    auto separator = left_node->split(buffer_right.data());
    ASSERT_EQ(left_node->count, n / 2 + 1);
    ASSERT_EQ(right_node->count, (n + 1) / 2);
    ASSERT_EQ(separator, n / 2 + 1);

    // Check keys & children of the left node
    left_keys = left_node->get_key_vector();
    left_values = left_node->get_pid_vector();
    ASSERT_EQ(left_keys.size(), left_node->count-1u);
    ASSERT_EQ(left_values.size(), left_node->count);
    for (auto i = 0; i < left_node->count-1; ++i) {
        ASSERT_EQ(left_keys[i], i+1);
    }
    for (auto i = 0; i < left_node->count; ++i) {
        ASSERT_EQ(left_values[i], i * 2);
    }

    // Check keys & children of the right node
    auto right_keys = right_node->get_key_vector();
    auto right_values = right_node->get_pid_vector();
    ASSERT_EQ(right_keys.size(), right_node->count-1u);
    ASSERT_EQ(right_values.size(), right_node->count);
    for (auto i = 0; i < right_node->count-1; ++i) {
        ASSERT_EQ(right_keys[i], left_node->count + i + 1);
    }
    for (auto i = 0; i < right_node->count; ++i) {
        ASSERT_EQ(right_values[i], (left_node->count + i) * 2);
    }
}

TEST(EytzingerTest, InnerNodeMerge) {
    AlignedVector buffer_left;
    AlignedVector buffer_right;
    buffer_left.resize(PageSize, 1024);
    buffer_right.resize(PageSize, 1024);

    auto n = BTree::InnerNode::kCapacity;
    // init left node
    auto left_node = new (buffer_left.data()) BTree::InnerNode();
    left_node->children[0] = Swip::fromPID(0);
    left_node->count = 1;
    // init right node
    auto right_node = new (buffer_right.data()) BTree::InnerNode();
    right_node->children[0] = Swip::fromPID(n & (~0x1));
    right_node->count = 1;

    // insert into left node
    for (uint32_t i = 1, j = 2; i < n/2; ++i, j = i * 2) {
        left_node->insert_split(i, j);
    }
    ASSERT_EQ(left_node->count, n/2);
    int32_t separator = n/2;
    // insert into right node
    for (uint32_t i = 1, j = 2; i < n/2; ++i, j = i * 2) {
        right_node->insert_split(n/2 + i, (n & (~0x1)) + j);
    }
    ASSERT_EQ(right_node->count, n/2);
    // merge
    left_node->merge(*right_node, separator);
    ASSERT_EQ(left_node->count, n & (~0x1));
    auto keys = left_node->get_key_vector();
    auto values = left_node->get_pid_vector();
    for (uint32_t i = 0; i < left_node->count-1u; ++i) {
        ASSERT_EQ(keys[i], i+1); 
    }
    for (auto i = 0; i < left_node->count; ++i) {
        ASSERT_EQ(values[i], i * 2);
    }
}

TEST(EytzingerTest, InnerNodeRebalance) {
    AlignedVector buffer_left;
    AlignedVector buffer_right;
    buffer_left.resize(PageSize, 1024);
    buffer_right.resize(PageSize, 1024);

    auto n = BTree::InnerNode::kCapacity;
    // init left node
    auto left_node = new (buffer_left.data()) BTree::InnerNode();
    left_node->children[0] = Swip::fromPID(0);
    left_node->count = 1;
    // insert into left node
    for (uint32_t i = 1, j = 2; i < n; ++i, j = i * 2) {
        left_node->insert_split(i, j);
    }
    ASSERT_EQ(left_node->count, n);
    // init separator
    KeyT separator = n;
    // init right node
    auto right_node = new (buffer_right.data()) BTree::InnerNode();
    right_node->children[0] = Swip::fromPID(2*n);
    right_node->count = 1;
    // insert into right node
    for (uint32_t i = 1, j = 2; i < n/2; ++i, j = i * 2) {
        right_node->insert_split(n + i, 2*n + j);
    }
    ASSERT_EQ(right_node->count, n/2);
    // rebalance
    BTree::InnerNode::rebalance(*left_node, *right_node, separator);
    int to_shift = (n+1)/4;
    ASSERT_EQ(left_node->count, n - to_shift);
    auto left_keys = left_node->get_key_vector();
    auto left_values = left_node->get_pid_vector();
    // check left node
    for (auto i = 0; i < left_node->count-1; ++i) {
        ASSERT_EQ(left_keys[i], i+1); 
    }
    for (auto i = 0; i < left_node->count; ++i) {
        ASSERT_EQ(left_values[i], i * 2);
    }
    // check separator
    ASSERT_EQ(separator, left_node->count);
    // check right node
    ASSERT_EQ(right_node->count, n/2 + to_shift);
    auto right_keys = right_node->get_key_vector();
    auto right_values = right_node->get_pid_vector();
    for (auto i = 0; i < right_node->count-1; ++i) {
        ASSERT_EQ(right_keys[i], i + left_node->count + 1); 
    }
    for (auto i = 0; i < right_node->count; ++i) {
        ASSERT_EQ(right_values[i], (i + left_node->count) * 2);
    }
}

TEST(EytzingerTest, InnerNodeRebalanceRightToLeft) {
    AlignedVector buffer_left;
    AlignedVector buffer_right;
    buffer_left.resize(PageSize, 1024);
    buffer_right.resize(PageSize, 1024);

    auto n = BTree::InnerNode::kCapacity;
    // init left node
    auto left_node = new (buffer_left.data()) BTree::InnerNode();
    left_node->children[0] = Swip::fromPID(0);
    left_node->count = 1;
    // insert into left node
    for (uint32_t i = 1, j = 2; i < n/2; ++i, j = i * 2) {
        left_node->insert_split(i, j);
    }
    ASSERT_EQ(left_node->count, n/2);
    // init separator
    KeyT separator = n/2;
    // init right node
    auto right_node = new (buffer_right.data()) BTree::InnerNode();
    right_node->children[0] = Swip::fromPID(n & (~0x1));
    right_node->count = 1;
    // insert into right node
    for (uint32_t i = 1, j = 2; i < n; ++i, j = i * 2) {
        right_node->insert_split(n/2 + i, (n & (~0x1)) + j);
    }
    ASSERT_EQ(right_node->count, n);
    // rebalance
    BTree::InnerNode::rebalance(*left_node, *right_node, separator);
    int to_shift = (n+1)/4;
    ASSERT_EQ(left_node->count, n/2 + to_shift);
    auto left_keys = left_node->get_key_vector();
    auto left_values = left_node->get_pid_vector();
    // check left node
    for (auto i = 0; i < left_node->count-1; ++i) {
        ASSERT_EQ(left_keys[i], i+1); 
    }
    for (auto i = 0; i < left_node->count; ++i) {
        ASSERT_EQ(left_values[i], i * 2);
    }
    // check separator
    ASSERT_EQ(separator, left_node->count);
    // check right node
    ASSERT_EQ(right_node->count, n - to_shift);
    auto right_keys = right_node->get_key_vector();
    auto right_values = right_node->get_pid_vector();
    for (auto i = 0; i < right_node->count-1; ++i) {
        ASSERT_EQ(right_keys[i], i + left_node->count + 1); 
    }
    for (auto i = 0; i < right_node->count; ++i) {
        ASSERT_EQ(right_values[i], (i + left_node->count) * 2);
    }
}
}
