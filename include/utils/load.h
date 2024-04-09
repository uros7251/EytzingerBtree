#ifndef INCLUDE_UTILS_LOAD_H
#define INCLUDE_UTILS_LOAD_H

#include <cstdint>
#include <random>
#include <algorithm>
#include <vector>
#include "utils/zipf.h"

namespace guidedresearch {

template <typename IntType>
struct IntLoad {
    std::vector<IntType> keys;
    std::vector<int> lookup_indices;
    size_t size, lookup_factor;

    IntLoad(size_t size, size_t lookup_factor) : keys(size), size(size), lookup_factor(lookup_factor) {
        std::iota(keys.begin(), keys.end(), 1);
        std::mt19937_64 g(0);
        std::shuffle(keys.begin(), keys.end(), g);
        lookup_indices.reserve(size*lookup_factor);
    }
    virtual ~IntLoad() = default;

    IntType lookup(int i) { return keys[lookup_indices.at(i)]; }
    virtual const std::string name() const = 0; 
};

template <typename IntType>
struct UniformLoad : public IntLoad<IntType> {
    UniformLoad(size_t size, size_t lookup_factor = 10) : IntLoad<IntType>(size, lookup_factor) {
        // std::random_device rd;
        std::mt19937 g(7251);
        std::uniform_int_distribution<IntType> dis(0, size-1);
        for (size_t i = 0; i < size*lookup_factor; ++i) {
            IntLoad<IntType>::lookup_indices.push_back(dis(g));
        }
    }

    const std::string name() const override { return "Uniform"; }
};

template <typename IntType>
struct ZipfianLoad : public IntLoad<IntType> {
    ZipfianLoad(size_t size, size_t lookup_factor = 10) : IntLoad<IntType>(size, lookup_factor) {
        ZipfDistribution zipf(size, 1); // alpha = 1
        // std::random_device rd;
        std::mt19937_64 engine(7251);
        for (size_t i = 0; i < size*lookup_factor; ++i) {
            IntLoad<IntType>::lookup_indices.push_back(zipf(engine)-1);
        }
    }

    const std::string name() const override { return "Zipfian"; }
};

} // namespace guidedresearch
#endif