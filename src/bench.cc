#include <cstdint>
#include <iostream>
#include <random>
#include <algorithm>
#include <chrono>
#include "btree/btree.h"
#include "utils/load.h"

using BufferManager = guidedresearch::BufferManager;
using KeyT = int64_t;

struct ValueT {
    static constexpr int size=5;
    int64_t a[size];
    ValueT() = default;
    ValueT(int64_t a) :a{} { this->a[0]=a; }
};

uint64_t rdtsc() {
  unsigned int lo, hi;
  __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
  return (static_cast<uint64_t>(hi) << 32) | lo;
}

template<typename BTree>
std::tuple<uint64_t,uint64_t, uint64_t, uint64_t> bench(guidedresearch::IntLoad<KeyT>& load) {
    BufferManager buffer_manager(BTree::kPageSize, 2*(load.size/BTree::LeafNode::kCapacity+1));
    BTree tree(0, buffer_manager);

    auto& keys = load.keys;
    const auto n = load.size;
    const auto m = load.size*load.lookup_factor;

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
    for (auto i = 0u; i < m; ++i) {
        auto key = load.lookup(i);
        auto result = tree.lookup(key);
        if (result) checksum ^= (*result).a[0];
        else {
            std::cerr << "lookup failed for key " << key << "!\n";
            return {0, 0, 0, 0};
        }
    }

    if (checksum) std::cerr << "checksum: " << checksum << "!\n";
    end = rdtsc();
    duration = end-start;
    auto lookup_time = duration/m;

    int64_t total = 0;
    KeyT min = static_cast<KeyT>(n/64), length = 0;
    auto num_iters = n>>10; // n/1024
    start = rdtsc();
    for (auto i=0u; i<num_iters; ++i) {
        length += static_cast<KeyT>(n/(2*num_iters));
        tree.traverse([&total]([[maybe_unused]] const KeyT &k, const ValueT &v) {
            total += v.a[0];
        }, min, min+length);
    }
    end = rdtsc();
    auto range_scan_time = (end-start)/num_iters;

    std::cerr << "rangesum: " << total << "!\n";
    // erase all values
    start = rdtsc();
    for (auto i = 0ul; i < n; ++i) {
        tree.erase(keys[i]);
    }
    end = rdtsc();
    duration = end-start;

    return {insert_time, lookup_time, range_scan_time, duration/n};
}

template<size_t PageSize = 4096>
void layout_comparison(guidedresearch::IntLoad<KeyT>& load) {
    uint64_t ref_insert, ref_lookup, ref_range_scan, ref_erase;
    uint64_t ref_insert_fast_insert, ref_lookup_fast_insert, ref_range_scan_fast_insert, ref_erase_fast_insert;
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::SORTED, guidedresearch::NodeLayout::SORTED>;
        auto label = "Standard";
        std::tie(ref_insert, ref_lookup, ref_range_scan, ref_erase) = bench<BTree>(load);
        std::cout << load.name() << "," << PageSize << "," << label << ",0," << ref_insert << "," << ref_lookup << "," << ref_range_scan << "," << ref_erase << "\n";
    }
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::SORTED, guidedresearch::NodeLayout::SORTED, true>;
        auto label = "Standard";
        std::tie(ref_insert_fast_insert, ref_lookup_fast_insert, ref_range_scan_fast_insert, ref_erase_fast_insert) = bench<BTree>(load);
        std::cout << load.name() << "," << PageSize << "," << label << ",1," << ref_insert_fast_insert << "," << ref_lookup_fast_insert << "," << ref_range_scan_fast_insert << "," << ref_erase_fast_insert << "\n";
        std::cerr << label << " speedup:\n"
            << "insert: " << static_cast<double>(ref_insert)/ref_insert_fast_insert
            << "x, lookup: " << static_cast<double>(ref_lookup)/ref_lookup_fast_insert
            << "x, range scan: " << static_cast<double>(ref_range_scan)/ref_range_scan_fast_insert
            << "x, erase: " << static_cast<double>(ref_erase)/ref_erase_fast_insert
            << "\n";
    }
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::EYTZINGER, guidedresearch::NodeLayout::EYTZINGER, false>;
        auto label = "Eytzinger";
        auto [insert, lookup, range_scan, erase] = bench<BTree>(load);
        std::cout << load.name() << "," << PageSize << "," << label << ",0," << insert << "," << lookup << "," << range_scan << "," << erase << "\n";
        std::cerr << label << " speedup:\n"
            << "insert: " << static_cast<double>(ref_insert)/insert
            << "x, lookup: " << static_cast<double>(ref_lookup)/lookup
            << "x, range scan: " << static_cast<double>(ref_range_scan)/range_scan
            << "x, erase: " << static_cast<double>(ref_erase)/erase
            << "\n";
    }
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::EYTZINGER, guidedresearch::NodeLayout::EYTZINGER, true>;
        auto label = "Eytzinger";
        auto [insert, lookup, range_scan, erase] = bench<BTree>(load);
        std::cout << load.name() << "," << PageSize << "," << label << ",1," << insert << "," << lookup << "," << range_scan << "," << erase << "\n";
        std::cerr << label << " speedup:\n"
            << "insert: " << static_cast<double>(ref_insert_fast_insert)/insert
            << "x, lookup: " << static_cast<double>(ref_lookup_fast_insert)/lookup
            << "x, range scan: " << static_cast<double>(ref_range_scan_fast_insert)/range_scan
            << "x, erase: " << static_cast<double>(ref_erase_fast_insert)/erase
            << "\n";
    }
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::EYTZINGER_SIMD, guidedresearch::NodeLayout::EYTZINGER_SIMD, false>;
        #if defined(__AVX512F__) && defined(__AVX512VL__)
        auto label = "n-ary Eytzinger + SIMD";
        #elif defined(__AVX__) && defined(__AVX2__)
        auto label = "n-ary Eytzinger + SIMD";
        #else
        auto label = "n-ary Eytzinger";
        #endif
        
        auto [insert, lookup, range_scan, erase] = bench<BTree>(load);
        std::cout << load.name() << "," << PageSize << "," << label << ",0," << insert << "," << lookup << "," << range_scan << "," << erase << "\n";
        std::cerr << label << " speedup:\n"
            << "insert: " << static_cast<double>(ref_insert)/insert
            << "x, lookup: " << static_cast<double>(ref_lookup)/lookup
            << "x, range scan: " << static_cast<double>(ref_range_scan)/range_scan
            << "x, erase: " << static_cast<double>(ref_erase)/erase
            << "\n";
    }
    {
        using BTree = guidedresearch::BTree<KeyT, ValueT, std::less<KeyT>, PageSize, guidedresearch::NodeLayout::EYTZINGER_SIMD, guidedresearch::NodeLayout::EYTZINGER_SIMD, true>;
        #if defined(__AVX512F__) && defined(__AVX512VL__)
        auto label = "n-ary Eytzinger + SIMD";
        #elif defined(__AVX__) && defined(__AVX2__)
        auto label = "n-ary Eytzinger + SIMD";
        #else
        auto label = "n-ary Eytzinger";
        #endif
        auto [insert, lookup, range_scan, erase] = bench<BTree>(load);
        std::cout << load.name() << "," << PageSize << "," << label << ",1," << insert << "," << lookup << "," << range_scan << "," << erase << "\n";
        std::cerr << label << " speedup:\n"
            << "insert: " << static_cast<double>(ref_insert_fast_insert)/insert
            << "x, lookup: " << static_cast<double>(ref_lookup_fast_insert)/lookup
            << "x, range scan: " << static_cast<double>(ref_range_scan_fast_insert)/range_scan
            << "x, erase: " << static_cast<double>(ref_erase_fast_insert)/erase
            << "\n";
    }
}

constexpr size_t page_sizes[] = {1<<12, 1<<14, 1<<16}; // 4K, 16K, 64K
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
    using KeyT = int64_t;
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
    std::ranges::shuffle(keys, g);
    
    auto time=rdtsc();
    auto checksum = 0u;
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
        guidedresearch::ZipfianLoad<KeyT> load(1<<22);
        page_size_comparison(load);
    }
    {
        guidedresearch::UniformLoad<KeyT> load(1<<22);
        page_size_comparison(load);
    }
}

void iterator() {
    using Iterator = guidedresearch::EytzingerIterator<8>;

    constexpr auto size = 2500, num_iters = 100000;
    auto checksum = 0u;
    auto start = rdtsc();
    for (int i=0; i<num_iters; ++i) {
        auto it = Iterator::begin(size);
        for (auto end = Iterator::end(size); it != end; ++it) {
            checksum ^= *it;
        }
    }
    auto avg = (rdtsc()-start)/static_cast<double>(num_iters);
    std::cerr << "checksum: " << checksum << "\n";
    std::cout << "avg time: " << avg << "\n";
}

int main() {
    #if defined(__AVX512F__) && defined(__AVX512VL__)
    std::cerr << "AVX512\n";
    #elif defined(__AVX__) && defined(__AVX2__)
    std::cerr << "AVX/AVX2\n";
    #else
    std::cerr << "NO AVX\n";
    #endif
    std::cout << "load,page_size,layout,fast_inserts,Insert,Lookup,Range_scan,Erase\n";
    // load_comparison();
    guidedresearch::UniformLoad<KeyT> load(1<<22);
    page_size_comparison(load);
    // layout_comparison<1u<<16>(load);
    
    return 0;
}