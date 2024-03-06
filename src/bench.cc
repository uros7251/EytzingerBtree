#include <cstdint>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
#include "btree/btree.h"

using BufferManager = guidedresearch::BufferManager;

uint64_t rdtsc() {
  unsigned int lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
}

template<typename T>
std::pair<uint64_t,uint64_t> bench() {
    BufferManager buffer_manager(1024, 100000);
    T tree(0, buffer_manager);
    auto n = 10000 * T::LeafNode::kCapacity;

    // Generate random non-repeating key sequence
    std::vector<uint64_t> keys(n);
    std::iota(keys.begin(), keys.end(), n);
    std::mt19937_64 engine(10);
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
        checksum ^= *tree.lookup(keys[i % n]);
    }
    end = rdtsc();
    duration = end-start;
    return {insert_time, duration/m};
}

int main() {
    uint64_t ref_insert, ref_lookup, eytzinger_insert, eytzinger_lookup, simd_insert, simd_lookup;
    {
        using BTree = guidedresearch::BTree<uint64_t, uint64_t, std::less<uint64_t>, 1024, guidedresearch::NodeLayout::SORTED>;
        std::tie(ref_insert, ref_lookup) = bench<BTree>();
        std::cout << "Standard BTree:\n"
            << ref_insert << "cyc/insert, " << ref_lookup << "cyc/lookup\n";
    }
    {
        using BTree = guidedresearch::BTree<uint64_t, uint64_t, std::less<uint64_t>, 1024, guidedresearch::NodeLayout::EYTZINGER>;
        std::tie(eytzinger_insert, eytzinger_lookup) = bench<BTree>();
        std::cout << "Eytzinger BTree:\n"
            << eytzinger_insert << "cyc/insert, " << eytzinger_lookup << "cyc/lookup\n";
    }
    {
        using BTree = guidedresearch::BTree<uint64_t, uint64_t, std::less<uint64_t>, 1024, guidedresearch::NodeLayout::EYTZINGER_SIMD>;
        std::tie(simd_insert, simd_lookup) = bench<BTree>();
        std::cout << "Eytzinger+SIMD BTree:\n"
            << simd_insert << "cyc/insert, " << simd_lookup << "cyc/lookup\n";
    }
    std::cout << "Eytzinger speedup:\n";
    std::cout << "insert: " << static_cast<double>(ref_insert)/eytzinger_insert
            << "x, lookup: " << static_cast<double>(ref_lookup)/eytzinger_lookup
            << "x\n";

    std::cout << "Eytzinger+SIMD speedup:\n";
    std::cout << "insert: " << static_cast<double>(ref_insert)/simd_insert
            << "x, lookup: " << static_cast<double>(ref_lookup)/simd_lookup
            << "x\n";

}