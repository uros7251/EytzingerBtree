#ifndef INCLUDE_UTILS_ZIPF_H
#define INCLUDE_UTILS_ZIPF_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>

namespace guidedresearch {

template <typename IntType = int>
class ZipfDistribution {
    std::uniform_real_distribution<> dis;
    std::vector<double> cdf;
public:
    ZipfDistribution(int N, double alpha) :dis(0.0, 1.0) {
        // generate the CDF
        cdf.reserve(N);
        double cdf_value = 0.0;
        for (int i = 1; i <= N; ++i) {
            cdf_value += 1 / std::pow(i, alpha);
            cdf.push_back(cdf_value);
        }
        for (int i = 0; i < N; ++i) {
            cdf[i] /= cdf_value; // normalization
        }
    }

    template<typename Generator>
    int operator()(Generator& gen) {
        // generate a random number in [0, 1]
        double random_value = dis(gen);
        // x = CDF^(-1)(random_value)
        return std::lower_bound(cdf.begin(), cdf.end(), random_value) - cdf.begin() + 1;
    }
};

} // namespace guidedresearch
#endif