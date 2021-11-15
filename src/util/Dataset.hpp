#pragma once

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <type_traits>
#include <unordered_set>

#include "util/data.hpp"
#include "util/enum.hpp"
#include "util/util.hpp"
#include "util/Exceptions.hpp"

template<typename Key, typename T>
struct Dataset
{
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<const Key, T>;
    using size_type = std::size_t;
    using domain_range = util::Range<key_type>;

    private:
    using entry_type = std::pair<Key, T>;

    std::vector<entry_type> entries_;   //< dataset entries
    size_type n_;                       //< dataset size
    domain_range dom_;                   //< dataset domain
    Distribution dist_;                 //< distribution
    uint64_t seed_;                     //< random generator seed
    std::mt19937_64 g_;                 //< random generator

    Dataset(size_type n, domain_range dom, Distribution dist, uint64_t seed = 42)
        : n_(n)
        , dom_(dom)
        , dist_(dist)
        , seed_(seed)
        , g_(std::mt19937_64(seed))
    { }

    Dataset(std::vector<entry_type> entries, size_type n, domain_range dom, Distribution dist, uint64_t seed = 42)
        : entries_(entries)
        , n_(n)
        , dom_(dom)
        , dist_(dist)
        , seed_(seed)
        , g_(std::mt19937_64(seed))
    {}
    public:
    static Dataset * Create_Uniform_Dense(size_type n, Key lower_bound)
    {
        static_assert(std::is_integral<key_type>::value, "Integral type for Key required.");
        static_assert(std::is_integral<T>::value, "Integral type for T required.");

        key_type upper_bound = lower_bound + n - 1;
        Dataset * data = new Dataset(n, {lower_bound, upper_bound}, Distribution::uniform_dense);
        auto & entries = data->entries();
        entries.reserve(n);
        key_type k = lower_bound;
        for (size_type i = 0; i < n; ++i, k++) entries.push_back({k, k});
        return data;
    }

    static Dataset * Create_Uniform_Sparse(size_type n, domain_range dom, uint64_t seed = 42)
    {
        static_assert(std::is_integral<key_type>::value, "Integral type for Key required.");
        static_assert(std::is_integral<T>::value, "Integral type for T required.");

        Dataset * data = new Dataset(n, dom, Distribution::uniform_sparse, seed);
        auto & entries = data->entries();
        entries.reserve(n);

        std::uniform_int_distribution<key_type> distr(dom.min, dom.max);
        std::unordered_set<key_type> unique_keys;
        while (entries.size() != n) {
            key_type k = distr(data->random_generator());
            auto unique = unique_keys.insert(k).second; // true if insertion took place, i.e. unique
            if (unique) entries.push_back({k, k});
        }
        data->sort();
        return data;
    }

    static Dataset * Create_Books(size_type n, uint64_t seed = 42)
    {
        static_assert(std::is_integral<key_type>::value, "Integral type for Key required.");
        static_assert(std::is_integral<T>::value, "Integral type for T required.");

        /* Check if wiki dataset is locally present. */
        const std::string filename("data/books_200M_uint64");
        if (not std::filesystem::exists(filename)) {
            std::cerr << "Books dataset does not exist at '" << filename
                      << "Execute download script in './data'"
                      << std::endl;
            throw FileNotFoundException();
        }

        /* Load binary data. */
        auto data = load_data<uint64_t>(filename);
        key_type min = data.front();
        key_type max = data.back();

        /* Create Dataset object. */
        Dataset * dataset = new Dataset(n, {min, max}, Distribution::books, seed);
        auto & entries = dataset->entries();
        entries.reserve(n);

        auto s = sample(data, n, dataset->random_generator());
        for (auto i = 0; i < n; ++i) entries.push_back({s[i], i});

        return dataset;
    }

    static Dataset * Create_Osm(size_type n, uint64_t seed = 42)
    {
        static_assert(std::is_integral<key_type>::value, "Integral type for Key required.");
        static_assert(std::is_integral<T>::value, "Integral type for T required.");

        /* Check if wiki dataset is locally present. */
        const std::string filename("data/osm_cellids_200M_uint64");
        if (not std::filesystem::exists(filename)) {
            std::cerr << "osm dataset does not exist at '" << filename
                      << "Execute download script in './data'"
                      << std::endl;
            throw FileNotFoundException();
        }

        /* Load binary data. */
        auto data = load_data<uint64_t>(filename);
        key_type min = data.front();
        key_type max = data.back();

        /* Create Dataset object. */
        Dataset * dataset = new Dataset(n, {min, max}, Distribution::osm, seed);
        auto & entries = dataset->entries();
        entries.reserve(n);

        auto s = sample(data, n, dataset->random_generator());
        for (auto i = 0; i < n; ++i) entries.push_back({s[i], i});

        return dataset;
    }

    static Dataset * Create_Explicit(std::vector<entry_type> entries, size_type n, domain_range dom, Distribution dist, uint64_t seed = 42)
    {
        return new Dataset(entries, n, dom, dist, seed);
    }

    /* Iterator. */
    using iterator = typename std::vector<entry_type>::iterator;
    using const_iterator = typename std::vector<entry_type>::const_iterator;

    iterator begin() { return entries_.begin(); }
    iterator end()   { return entries_.end(); }
    const_iterator cbegin() const { return entries_.cbegin(); }
    const_iterator cend()   const { return entries_.cend(); }

    void shuffle() { std::shuffle(begin(), end(), g_); }
    void sort() { std::sort(begin(), end()); }
    bool is_sorted() { return std::is_sorted(begin(), end()); }

    std::vector<entry_type> & entries() { return entries_; }
    size_type size() const { assert(entries_.size() == n_); return n_; }
    domain_range domain() const { return dom_; }
    Distribution distribution() const { return dist_; }
    uint64_t seed() const { return seed_; }
    std::mt19937_64 & random_generator() { return g_; }


    entry_type & operator[](size_type idx) { return entries_[idx]; }
    const entry_type & operator[](size_type idx) const { return entries_[idx]; }

    void toFile(const std::string& filename, const std::string& delimiter=":") {
        std::ofstream file;
        file.open(filename);
        file << "Size" << delimiter << n_ << std::endl;
        file << "DomainMin" << delimiter << dom_.min << std::endl;
        file << "DomainMax" << delimiter << dom_.max << std::endl;
        file << "Distribution" << delimiter << toString(dist_) << std::endl;
        file << "Seed" << delimiter << seed_ << std::endl;
        file << "Entries" << delimiter << std::endl;
        for (auto e : entries_)
            file << e.first << delimiter << e.second << std::endl;
        file.close();
    }

    static Dataset<Key, T>* fromFile(const std::string& filename, const std::string& delimiter=":") {
        std::ifstream file(filename);
        if (file.fail()) throw std::invalid_argument("Dataset file does not exist.");
        bool sizeFound=false, domMinFound=false, domMaxFound=false, distFound=false, seedFound=false, entriesFound=false;
        size_type n = 0;
        key_type dom_min = key_type(), dom_max = key_type();
        Distribution dist = Distribution::uniform_dense;
        uint64_t seed = 42;
        std::vector<entry_type> entries;
        for(std::string line; file.is_open() && std::getline(file, line); ) {
            std::string firstToken = line.substr(0, line.find(delimiter));
            std::string secondToken = line.erase(0, line.find(delimiter) + delimiter.length());
            if (firstToken == "Size") {
                n = std::stoul(secondToken);
                sizeFound = true;
            } else if (firstToken == "DomainMin") {
                dom_min = key_type(std::stold(secondToken));
                domMinFound = true;
            } else if (firstToken == "DomainMax") {
                dom_max = key_type(std::stold(secondToken));
                domMaxFound = true;
            } else if (firstToken == "Distribution") {
                dist = distributionTypeFromString(secondToken);
                distFound = true;
            } else if (firstToken == "Seed") {
                seed = std::stoul(secondToken);
                seedFound = true;
            } else if (firstToken == "Entries") {
                for (uint i=0; i < n; i++) {
                    if (not getline(file, line))
                        throw std::invalid_argument("file does not contain the specified amount of entries");
                    key_type key = key_type(std::stold(line.substr(0, line.find(delimiter))));
                    mapped_type value = mapped_type(std::stold(line.erase(0, line.find(delimiter) + delimiter.length())));
                    entries.push_back(std::make_pair(key, value));
                }
                entriesFound = true;
            } else
                throw std::invalid_argument("Unknown parameter found.");
        }
        if (sizeFound && domMinFound && domMaxFound && distFound && seedFound && entriesFound)
            return new Dataset<Key, T>(entries, n, {dom_min, dom_max}, dist, seed);
        else
            throw std::invalid_argument("Not all necessary parameters present in file");
    }
};
