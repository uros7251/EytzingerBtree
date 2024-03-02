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
    {
        using BTree = guidedresearch::BTree<uint64_t, uint64_t, std::less<uint64_t>, 1024, guidedresearch::NodeLayout::SORTED>;
        auto cpl = bench<BTree>();
        std::cout << "Standard BTree:\n"
            << cpl.first << "cyc/insert, " << cpl.second << "cyc/lookup\n";
    }
    {
        using BTree = guidedresearch::BTree<uint64_t, uint64_t, std::less<uint64_t>, 1024, guidedresearch::NodeLayout::EYTZINGER>;
        auto cpl = bench<BTree>();
        std::cout << "Eytzinger BTree:\n"
            << cpl.first << "cyc/insert, " << cpl.second << "cyc/lookup\n";
    }

}