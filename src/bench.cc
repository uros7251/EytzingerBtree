#include <cstdint>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
#include "btree/btree.h"

using BufferManager = guidedresearch::BufferManager;
using KeyT = int32_t;
using ValueT = int64_t;

uint64_t rdtsc() {
  unsigned int lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
}

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize, guidedresearch::NodeLayout inner_layout, guidedresearch::NodeLayout leaf_layout = guidedresearch::NodeLayout::SORTED>
std::tuple<uint64_t,uint64_t, uint64_t> bench() {
    using BTree = guidedresearch::BTree<KeyT, ValueT, ComparatorT, PageSize, inner_layout, leaf_layout>;

    constexpr auto leaf_node_count = 1<<10; // ~ 1k leaf nodes
    auto n = leaf_node_count*BTree::LeafNode::kCapacity; // ~ 1M keys
    BufferManager buffer_manager(PageSize, 2*leaf_node_count);
    BTree tree(0, buffer_manager);

    // Generate random non-repeating key sequence
    std::vector<int32_t> keys(n);
    std::iota(keys.begin(), keys.end(), n);
    std::mt19937_64 engine(0);
    std::shuffle(keys.begin(), keys.end(), engine);

    auto start {rdtsc()};
    // Insert values
    for (auto i = 0ul; i < n; ++i) {
        tree.insert(keys[i], 2 * keys[i]);
    }
    auto end {rdtsc()};
    auto duration = end-start;
    auto insert_time = duration/n;

    // Lookup all values
    uint64_t checksum = 0ul;
    auto m = n*10u;
    start = rdtsc();
    for (auto i = 0ul; i < m; ++i) {
        checksum ^= *tree.lookup(keys[i%n]);
    }
    end = rdtsc();
    duration = end-start;
    auto lookup_time = duration/m;

    start = rdtsc();
    for (auto i = 0ul; i < n; ++i) {
        tree.erase(keys[i]);
    }
    end = rdtsc();
    duration = end-start;

    return {insert_time, lookup_time, duration/m};
}

constexpr size_t page_sizes[] = {1024, 4096, 16384, 65536};

template<unsigned N>
void page_size_comparison() {
    constexpr unsigned n = std::min(4u, N);

    std::cout << "Page size: " << page_sizes[4u-n] << "\n";
    
    auto [insert, lookup, erase] = bench<KeyT, ValueT, std::less<KeyT>, page_sizes[4u-n], guidedresearch::NodeLayout::SORTED>();
    std::cout << "Standard BTree:\n"
        << insert << "cyc/insert, " << lookup << "cyc/lookup, " << erase << "cyc/erase\n";

    std::tie(insert, lookup, erase) = bench<KeyT, ValueT, std::less<KeyT>, page_sizes[4u-n], guidedresearch::NodeLayout::EYTZINGER>();
    std::cout << "Eytzinger BTree:\n"
        << insert << "cyc/insert, " << lookup << "cyc/lookup, " << erase << "cyc/erase\n";

    std::tie(insert, lookup, erase) = bench<KeyT, ValueT, std::less<KeyT>, page_sizes[4u-n], guidedresearch::NodeLayout::EYTZINGER_SIMD>();
    std::cout << "Eytzinger+SIMD BTree:\n"
        << insert << "cyc/insert, " << lookup << "cyc/lookup, " << erase << "cyc/erase\n\n";
    
    page_size_comparison<n-1>();
}

template<>
void page_size_comparison<0u>() {
    // break recursion
}

void inner_node_comparison() {
    constexpr auto PageSize = 1<<14; // 16KB 

    using AlignedVector = guidedresearch::AlignedVector;
    using KeyT = int32_t;
    using ValueT = int64_t;
    using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::SORTED>;
    using Swip = guidedresearch::Swip;

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
    
    auto time=rdtsc();
    int checksum = 0;
    constexpr auto N = 1<<10;
    for (auto j=0; j<N; ++j) {
    for (auto i : keys) {
        auto [index, found] = node->lower_bound(i);
        if (found) checksum ^= index;
    }
    }
    time = rdtsc()-time;
    if (checksum) std::cerr << "checksum: " << checksum << "!\n";
    std::cout << "Time:" << time/(n*N) << " cycles\n";
}

void layout_comparison() {
    constexpr size_t page_size = 4096<<1;
    uint64_t ref_insert, ref_lookup, ref_erase,
            eytzinger_insert, eytzinger_lookup, eytzinger_erase,
            simd_insert, simd_lookup, simd_erase;
    {
        std::tie(ref_insert, ref_lookup, ref_erase) = bench<KeyT, ValueT, std::less<KeyT>, page_size, guidedresearch::NodeLayout::SORTED>();
        std::cout << "Standard BTree:\n"
            << ref_insert << "cyc/insert, " << ref_lookup << "cyc/lookup, " << ref_erase << "cyc/erase\n";
    }
    {
        std::tie(eytzinger_insert, eytzinger_lookup, eytzinger_erase) = bench<KeyT, ValueT, std::less<KeyT>, page_size, guidedresearch::NodeLayout::EYTZINGER, guidedresearch::NodeLayout::SORTED>();
        std::cout << "Eytzinger BTree:\n"
            << eytzinger_insert << "cyc/insert, " << eytzinger_lookup << "cyc/lookup, " << eytzinger_erase << "cyc/erase\n";
    }
    {
        std::tie(simd_insert, simd_lookup, simd_erase) = bench<KeyT, ValueT, std::less<KeyT>, page_size, guidedresearch::NodeLayout::EYTZINGER_SIMD, guidedresearch::NodeLayout::SORTED>();
        std::cout << "Eytzinger+SIMD BTree:\n"
            << simd_insert << "cyc/insert, " << simd_lookup << "cyc/lookup, " << simd_erase << "cyc/erase\n";
    }
    std::cout << "Eytzinger speedup:\n";
    std::cout << "insert: " << static_cast<double>(ref_insert)/eytzinger_insert
            << "x, lookup: " << static_cast<double>(ref_lookup)/eytzinger_lookup
            << "x, erase: " << static_cast<double>(ref_erase)/eytzinger_erase
            << "\n";

    std::cout << "Eytzinger+SIMD speedup:\n";
    std::cout << "insert: " << static_cast<double>(ref_insert)/simd_insert
            << "x, lookup: " << static_cast<double>(ref_lookup)/simd_lookup
            << "x, erase: " << static_cast<double>(ref_erase)/simd_erase
            << "\n";
}

int main() {
    //page_size_comparison<4>();
    layout_comparison();
    //inner_node_comparison();
    return 0;
}