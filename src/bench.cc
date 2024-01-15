#include "btree/btree.h"
#include <cstdint>
#include <iostream>

using BTree = guidedresearch::BTree<uint64_t, uint64_t, std::less<uint64_t>, 1024, guidedresearch::NodeLayout::EYTZINGER>;
using Swip = guidedresearch::Swip;

uint64_t rdtsc() {
  unsigned int lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return ((uint64_t)hi << 32) | lo;
}

int main() {
    std::vector<std::byte> buffer;
    buffer.resize(1024);

    auto n = BTree::InnerNode::kCapacity;

    auto start = rdtsc();

    int iters = 100;

    for (int i=0; i<iters; ++i) {
        // init inner node
        auto node = new (buffer.data()) BTree::InnerNode();
        node->children[0] = Swip::fromPID(0);
        node->count = 1;
        // Insert children into the leaf nodes
        for (uint32_t i = n, j = n*2; i > 0; --i, j = i * 2) {
            Swip s = Swip::fromPID(j);
            node->insert_split(i, s);
        }
    }

    auto end = rdtsc();

    std::cout << "Cycles per insertion: " << (end-start)/(iters*n) << "\n";
}