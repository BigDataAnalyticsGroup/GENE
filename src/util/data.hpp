#pragma once

#include <unordered_set>
#include <random>


/*======================================================================================================================
 * Load Datasets
 *====================================================================================================================*/

/** Loads a dataset from a binary file.
 * Author: Marcel Maltry <marcel.maltry@bigdata.uni-saarlan.de>
 * */
template<typename T>
std::vector<T> load_data(const std::string &filename)
{
    using data_type = T;

    std::ifstream istrm(filename, std::ios::binary);

    if (not istrm.is_open()) throw std::runtime_error("Could not open "+filename+".");

    uint64_t n_data;
    istrm.read(reinterpret_cast<char *>(&n_data), sizeof(n_data));

    std::vector<T> data;
    data.reserve(n_data);

    data_type d;
    for (std::size_t i = 0; istrm.good(); ++i) {
        istrm.read(reinterpret_cast<char *>(&d), sizeof(d));
        if (istrm.eof()) break;
        data.push_back(d);
    }
    istrm.close();

    assert(data.size() == n_data);
    assert(std::is_sorted(data.begin(), data.end()));

    return data;
}


/** Draw uniformly at random `k` elements from `data`.
 *  The resulting distribution should match the distribution of `data`. */
template<typename T, typename RandomGenerator>
std::vector<T> sample(const std::vector<T> &data, uint64_t k, RandomGenerator &g)
{
    using vector_type = std::vector<T>;
    vector_type res;
    res.reserve(k);
#if 0
    /* Equal width sampling. */
    uint64_t width = data.size() / k;
    uint64_t i = 0;
    while (i < data.size() && res.size() < k) {
        res.push_back(data[i]);
        i += width;
    }
#else
    /* Uniformly at random. */
    std::unordered_set<T> unique_values;
    unique_values.reserve(k);
    std::uniform_int_distribution<T> distr(0, data.size() - 1);

    while (res.size() < k) {
        T v = data[distr(g)];  // draw uniformly at random some value from data
        auto unique = unique_values.insert(v).second; // true if insertion took place, i.e. unique
        if (unique) res.push_back(v);
    }
    std::sort(res.begin(), res.end());
#endif
    return res;
}
