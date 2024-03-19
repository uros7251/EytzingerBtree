#include <cstdint>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
#include "btree/btree.h"
#include "utils/load.h"

using BufferManager = guidedresearch::BufferManager;
using KeyT = int32_t;
using ValueT = int64_t;

uint64_t rdtsc() {
  unsigned int lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
}

template<typename BTree>
std::tuple<uint64_t,uint64_t, uint64_t> bench(guidedresearch::IntLoad<KeyT>& load) {
    BufferManager buffer_manager(BTree::kPageSize, 2*(load.size/BTree::LeafNode::kCapacity+1));
    BTree tree(0, buffer_manager);

    auto& keys = load.keys;
    auto n = load.size;
    auto m = load.size*load.lookup_factor;

    auto start {rdtsc()};
    // Insert values
    for (auto key: keys) {
        tree.insert(key, 2 * key);
    }
    std::cerr << "Tree depth: " << tree.depth() << "\n";

    auto end {rdtsc()};
    auto duration = end-start;
    auto insert_time = duration/n;

    // Lookup all values
    uint64_t checksum = 0ul;
    start = rdtsc();
    for (auto i = 0ul; i < m; ++i) {
        auto key = load.lookup(i);
        auto result = tree.lookup(key);
        if (result) checksum ^= *result;
        else {
            std::cerr << "lookup failed for key " << key << "!\n";
            return {0, 0, 0};
        }
    }
    if (checksum) std::cerr << "checksum: " << checksum << "!\n";
    end = rdtsc();
    duration = end-start;
    auto lookup_time = duration/m;

    start = rdtsc();
    for (auto i = 0ul; i < n; ++i) {
        tree.erase(keys[i]);
    }
    end = rdtsc();
    duration = end-start;

    return {insert_time, lookup_time, duration/n};
}

template<size_t PageSize = 4096>
void layout_comparison(guidedresearch::IntLoad<KeyT>& load) {
    uint64_t ref_insert, ref_lookup, ref_erase,
            eytzinger_insert, eytzinger_lookup, eytzinger_erase,
            simd_insert, simd_lookup, simd_erase;
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::SORTED>;
        std::tie(ref_insert, ref_lookup, ref_erase) = bench<BTree>(load);
        std::cout << load.name() << "|" << PageSize << "|Standard|" << ref_insert << "|" << ref_lookup << "|" << ref_erase << "|\n";
    }
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::EYTZINGER, guidedresearch::NodeLayout::EYTZINGER>;
        std::tie(eytzinger_insert, eytzinger_lookup, eytzinger_erase) = bench<BTree>(load);
        std::cout << load.name() << "|" << PageSize << "|Eytzinger|" << eytzinger_insert << "|" << eytzinger_lookup << "|" << eytzinger_erase << "|\n";
    }
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::EYTZINGER_SIMD, guidedresearch::NodeLayout::EYTZINGER_SIMD>;
        std::tie(simd_insert, simd_lookup, simd_erase) = bench<BTree>(load);
        std::cout << load.name() << "|" << PageSize << "|Eytzinger+SIMD|" << simd_insert << "|" << simd_lookup << "|" << simd_erase << "|\n";
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

constexpr size_t page_sizes[] = {1024, 4096, 16384};
constexpr unsigned ARRAY_SIZE = sizeof(page_sizes)/sizeof(size_t);

template<unsigned N = ARRAY_SIZE>
void page_size_comparison(guidedresearch::IntLoad<KeyT>& load) {
    constexpr unsigned n = std::min(ARRAY_SIZE, N);
    layout_comparison<page_sizes[ARRAY_SIZE-n]>(load);
    page_size_comparison<n-1>(load);
}

template<>
void page_size_comparison<0u>([[maybe_unused]] guidedresearch::IntLoad<KeyT>& load) {
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

void load_comparison() {
    {
        guidedresearch::ZipfianLoad<KeyT> load(1<<20);
        page_size_comparison(load);
    }
    {
        guidedresearch::UniformLoad<KeyT> load(1<<20);
        page_size_comparison(load);
    }
}
int main() {
    // load_comparison();
    guidedresearch::UniformLoad<KeyT> load(1<<20);
    layout_comparison(load);
    //inner_node_comparison(load);
    return 0;
}