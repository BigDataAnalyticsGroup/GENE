#pragma once

#include <algorithm>
#include <cassert>
#include <err.h>
#include <optional>

#include "nodes/Node.hpp"
#include "util/Exceptions.hpp"
#include "util/enum.hpp"

template<typename Key, typename Value>
struct HorizontalPartition;

#define IMPORT(WHAT) using WHAT = typename Base::WHAT

/** Strategy design pattern. */

/** Provides the algorithms to search for a given key in a `Node` inside a `HorizontalPartition`. */
template<typename Key, typename Value>
struct IndexSearchMethod
{
    using key_type = Key;
    using mapped_type = Value;
    using size_type = std::size_t;
    using node_type = Node<key_type, mapped_type>;

    SearchType searchType;

    IndexSearchMethod() : searchType(SearchType::UnknownSearchMethod) {}
    IndexSearchMethod(SearchType s) : searchType(s) {}

    /** Finds the corresponding mapped_type according to an equality comparison. */
    virtual std::optional<mapped_type> lookup(node_type *node, key_type &key) = 0;

    virtual std::vector<mapped_type> find_range(node_type *node, uint64_t &results, const key_type &lower, const key_type &upper) = 0;
    /** Finds the corresponding mapped_type (partition) according to the insertion semantics. */
    virtual std::optional<mapped_type> choose_subset(node_type *node, key_type &key) = 0;

    virtual ~IndexSearchMethod() { };

    virtual IndexSearchMethod<Key, Value>* clone() = 0;
    virtual void print() = 0;

    SearchType getSearchType() { return searchType; }
};

/*======================================================================================================================
 * DefaultSearch
 *====================================================================================================================*/
/** Utilizes the default `get` method of the underlying data structure. */
template<typename Key, typename Value>
struct DefaultSearch : IndexSearchMethod<Key, Value>
{
    using Base = IndexSearchMethod<Key, Value>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(node_type);

    DefaultSearch() : Base(SearchType::DefaultSearchMethod) {}

    void print() override { std::cout << "Default" << std::endl; }

    std::optional<mapped_type> lookup(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    return soa->get(key);
                    break;
                }
            case Node_type::Tree:
                {
                    auto tree = static_cast<TreeNode<key_type, mapped_type>*>(node);
                    return tree->get(key);
                    break;
                }
            case Node_type::Hash:
                {
                    auto hash = static_cast<HashNode<key_type, mapped_type>*>(node);
                    return hash->get(key);
                    break;
                }
            default:
                break;
        }
        std::cerr << "DefaultSearch strategy not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::vector<mapped_type> find_range(node_type *node, uint64_t &results, const key_type &lower, const key_type &upper) override {
        std::vector<mapped_type> res;
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    std::cerr << "DefaultSearch::find_range() for SoA not implemented yet!" << std::endl;
                    return res;
                }
            case Node_type::Tree:
                {
                    std::cerr << "DefaultSearch::find_range() for TreeNode not implemented yet!" << std::endl;
                    return res;
                }
            case Node_type::Hash:
                {
                    auto hash = static_cast<HashNode<key_type, mapped_type>*>(node);
                    auto table = hash->entries();
                    if (table.size() <= (upper - lower)) {
                        /* Hash table size smaller than range query domain -> scan HT */
                        for (auto e : table)
                            if ((e.first >= lower) and (e.first <= upper)) ++results;
                    } else {
                        /* Iterate over range domain. */
                        for (auto i = lower; i <= upper; ++i) if (hash->get(i)) ++results;
                    }
                    return res;
                }
            default: break;
        }
        std::cerr << "DefaultSearch::find_range() not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::optional<mapped_type> choose_subset(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    std::cerr << "DefaultSearch::choose_subset() for SoA not implemented yet!" << std::endl;
                    break;
                }
            case Node_type::Tree:
                {
                    std::cerr << "DefaultSearch::choose_subset() for TreeNode not implemented yet!" << std::endl;
                    break;
                }
            case Node_type::Hash:
                {
                    std::cerr << "DefaultSearch::choose_subset() for HashNode not implemented yet!" << std::endl;
                    break;
                }
            default: break;
        }
        (void) key;
        if (auto hash = dynamic_cast<HashNode<key_type, mapped_type>*>(node)) { return std::nullopt; }
        std::cerr << "DefaultSearch choose_subset not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    IndexSearchMethod<key_type, mapped_type> * clone() override { return new DefaultSearch<key_type, mapped_type>(); }
};


/*======================================================================================================================
 * LinearSearch
 *====================================================================================================================*/
template<typename Key, typename Value>
struct LinearSearch : IndexSearchMethod<Key, Value>
{
    using Base = IndexSearchMethod<Key, Value>;
    IMPORT(size_type);
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(node_type);

    LinearSearch() : Base(SearchType::LinearSearchMethod) {}

    void print() override { std::cout << "Linear" << std::endl; }

    std::optional<mapped_type> lookup(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    for (std::size_t i = 0, end = soa->size(); i < end; ++i)
                        if (key == keys[i]) return std::optional<mapped_type>(values[i]);
                    return std::nullopt;
                    break;
                }
            case Node_type::Tree:
                {
                    auto tree = static_cast<TreeNode<key_type, mapped_type>*>(node);
                    auto &entries = tree->entries();
                    for (auto e = entries.begin(), end = entries.end(); e != end; ++e)
                        if (key == e->first) return std::optional<mapped_type>(e->second);
                    return std::nullopt;
                    break;
                }
            case Node_type::Hash:
                {
                    auto hash = static_cast<HashNode<key_type, mapped_type>*>(node);
                    auto &entries = hash->entries();
                    for (auto e = entries.begin(), end = entries.end(); e != end; ++e)
                        if (key == e->first) return std::optional<mapped_type>(e->second);
                    return std::nullopt;
                    break;
                }
            default: break;
        }
        std::cerr << "LinearSearch strategy not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::vector<mapped_type> find_range(node_type *node, uint64_t &results, const key_type &lower, const key_type &upper) override {
        std::vector<mapped_type> res;
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    if constexpr (std::is_pointer_v<mapped_type>) {
                        auto &keys = soa->keys();
                        auto &values = soa->values();
                        key_type x_1 = 0;
                        for (size_type i = 0, end = soa->size(); i < end; ++i) {
                            /* Check if ranges overlap:
                             *  Query range: [lower, upper]
                             *  Key range: [x_1, x_2]
                             *      Concretely: [0, keys[0]], [keys[0]+1, keys[1]], ...
                             * -> x_1 <= upper && lower <= x_2
                             */
                            key_type x_2 = keys[i];
                            /* If false, all consecutive checks will be false as well. */
                            if (x_1 <= upper) {
                                if (lower <= x_2) {
                                    res.push_back(values[i]);
                                    x_1 = x_2 + 1;
                                }
                            } else { /* early abort. */ return res; }
                        }
                    } else {
                        auto &keys = soa->keys();
                        size_type i = 0;
                        while (keys[i] < lower) ++i; // set `i` to first element definitely in range
                        for (auto end = soa->size(); i < end; ++i) {
                            /* We know that keys[i] >= lower due to while loop. */
                            if (keys[i] <= upper) ++results;
                            else break; // early abort
                        }
                    }
                    return res;
                }
            case Node_type::Tree:
                {
                    auto tree = static_cast<TreeNode<key_type, mapped_type>*>(node);
                    if constexpr (std::is_pointer_v<mapped_type>) {
                        auto &entries = tree->entries();
                        key_type x_1 = 0;
                        for (auto it = entries.begin(), end = entries.end(); it != end; ++it) {
                            key_type x_2 = it->first;
                            if (x_1 <= upper) {
                                if (lower <= x_2) {
                                    res.push_back(it->second);
                                    x_1 = x_2 + 1;
                                }
                            } else { return res; }
                        }
                    } else {
                        auto &entries = tree->entries();
                        for (auto it = entries.begin(), end = entries.end(); it != end; ++it) {
                            key_type key = it->first;
                            if (key >= lower and key <= upper) ++results;
                        }
                    }
                    return res;
                }
            case Node_type::Hash:
                {
                    if constexpr (std::is_pointer_v<mapped_type>) {
                        /* There are no nested horizontal partitions in hash nodes. */
                        return res;
                    } else {
                        auto hash = static_cast<HashNode<key_type, mapped_type>*>(node);
                        auto &entries = hash->entries();
                        for (auto it = entries.begin(), end = entries.end(); it != end; ++it) {
                            key_type key = it->first;
                            if (key >= lower and key <= upper) ++results;
                        }
                    }
                    return res;
                }
            default: break;
        }
        std::cerr << "LinearSearch::find_range() not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::optional<mapped_type> choose_subset(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    for (std::size_t i = 0, end = soa->size(); i < end; ++i) {
                        auto &keys = soa->keys();
                        auto &values = soa->values();
                        if (key <= keys[i]) return std::optional<mapped_type>(values[i]);
                    }
                    return std::nullopt;
                    break;
                }
            case Node_type::Tree:
                {
                    auto tree = static_cast<TreeNode<key_type, mapped_type>*>(node);
                    auto &entries = tree->entries();
                    for (auto e = entries.begin(), end = entries.end(); e != end; ++e) {
                        if (key <= e->first) return std::optional<mapped_type>(e->second);
                    }
                    return std::nullopt;
                    break;
                }
            case Node_type::Hash:
                {
                    return std::nullopt;
                    break;
                }
            default: break;
        }
        std::cerr << "LinearSearch choose_subset strategy not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    IndexSearchMethod<key_type, mapped_type>* clone() override { return new LinearSearch<key_type, mapped_type>(); }
};

/*======================================================================================================================
 * BinarySearch
 *====================================================================================================================*/
template<typename Key, typename Value>
struct BinarySearch : IndexSearchMethod<Key, Value>
{
    using Base = IndexSearchMethod<Key, Value>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(node_type);
    IMPORT(size_type);

    BinarySearch() : Base(SearchType::BinarySearchMethod) {}

    void print() override { std::cout << "Binary" << std::endl; }

    std::optional<mapped_type> lookup(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    auto entry = std::lower_bound(keys.begin(), keys.end(), key);
                    if (entry != keys.end() and key == *entry) {
                        std::size_t idx = std::distance(keys.begin(), entry);
                        return std::optional<mapped_type>(values[idx]);
                    }
                    return std::nullopt;
                    break;
                }
            case Node_type::Tree:
                {
                    auto tree = static_cast<TreeNode<key_type, mapped_type>*>(node);
                    auto &entries = tree->entries();
                    auto entry = entries.lower_bound(key);
                    if (entry != entries.end() and key == entry->first)
                        return std::optional<mapped_type>(entry->second);
                    return std::nullopt;
                    break;
                }
            case Node_type::Hash:
                {
                    std::cerr << "BinarySearch for HashNode not supported!" << std::endl;
                    break;
                }
            default: break;
        }
        std::cerr << "BinarySearch strategy not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::vector<mapped_type> find_range(node_type *node, uint64_t &results, const key_type &lower, const key_type &upper) override {
        std::vector<mapped_type> res;
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    if constexpr (std::is_pointer_v<mapped_type>) {
                        auto &keys = soa->keys();
                        auto &values = soa->values();
                        auto key_it = std::lower_bound(keys.begin(), keys.end(), lower);
                        auto value_it = values.begin();
                        std::advance(value_it, std::distance(keys.begin(), key_it));
                        size_type x_1 = 0;
                        if (key_it != keys.begin()) x_1 = *(std::prev(key_it)) + 1;
                        for (auto end = keys.end(); key_it != end; ++key_it, ++value_it) {
                            if (x_1 > upper) return res; // early abort
                            res.push_back(*value_it);
                            x_1 = *key_it + 1;
                        }
                    } else {
                        auto &keys = soa->keys();
                        auto key_it = std::lower_bound(keys.begin(), keys.end(), lower);
                        for (; key_it != keys.end(); ++key_it) {
                            /* We know that *key_it >= lower due to std::lower_bound. */
                            if (*key_it <= upper) ++results;
                            else break; // early abort
                        }
                    }
                    return res;
                }
            case Node_type::Tree:
                {
                    auto tree = static_cast<TreeNode<key_type, mapped_type>*>(node);
                    if constexpr (std::is_pointer_v<mapped_type>) {
                        auto &entries = tree->entries();
                        auto it = entries.lower_bound(lower);
                        size_type x_1 = 0;
                        if (it != entries.begin()) x_1 = std::prev(it)->first + 1;
                        for (auto end = entries.end(); it != end; ++it) {
                            if (x_1 > upper) return res; // early abort
                            res.push_back(it->second);
                            x_1 = it->first + 1;
                        }
                    } else {
                        auto &entries = tree->entries();
                        auto it = entries.lower_bound(lower);
                        for (auto end = entries.end(); it != end; ++it) {
                            if (it->first > upper) return res; // early abort
                            ++results;
                        }
                    }
                    return res;
                }
            case Node_type::Hash:
                {
                    return res;
                }
            default: return res;
        }
        std::cout << "BinarySearch::find_range() not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
        return res;
    }

    std::optional<mapped_type> choose_subset(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    auto entry = std::lower_bound(keys.begin(), keys.end(), key);
                    if (entry != keys.end()) {
                        std::size_t idx = std::distance(keys.begin(), entry);
                        return std::optional<mapped_type>(values[idx]);
                    }
                    return std::nullopt;
                    break;
                }
            case Node_type::Tree:
                {
                    auto tree = static_cast<TreeNode<key_type, mapped_type>*>(node);
                    auto entries = tree->entries();
                    auto entry = entries.lower_bound(key);
                    if (entry != entries.end())
                        return std::optional<mapped_type>(entry->second);
                    return std::nullopt;
                    break;
                }
            case Node_type::Hash:
                {
                    break;
                }
            default: break;
        }
        std::cerr << "BinarySearch strategy not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    IndexSearchMethod<key_type, mapped_type>* clone() override { return new BinarySearch<key_type, mapped_type>(); }
};

/*======================================================================================================================
 * InterpolationSearch
 *====================================================================================================================*/
template<typename Key, typename Value>
struct InterpolationSearch : IndexSearchMethod<Key,Value>
{
    using Base = IndexSearchMethod<Key, Value>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(node_type);
    IMPORT(size_type);

    InterpolationSearch() : Base(SearchType::InterpolationSearchMethod) {}

    void print() override { std::cout << "Interpolation" << std::endl; }

    std::optional<mapped_type> lookup(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    size_type left = 0;
                    size_type right = soa->size() - 1;
                    // if right < left: found empty node Update: Does not work if right is an unsigned type as 0-1 will then result in a large, positive number
                    if (soa->size() == 0)
                        return std::nullopt;
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    Key current_key = keys[right];
                    if (key > current_key) return std::nullopt;
                    if (key == current_key) return std::optional<mapped_type>(values[right]);
                    Key key_left = keys[left];
                    Key key_right = keys[right];
                    while (right > left and key >= key_left and key <= key_right) {
                        double slope = (double) (right - left) / (key_right - key_left);
                        size_type expected = left + (key - key_left) * slope;
                        assert(expected >= 0);
                        assert(expected < soa->size());
                        assert(expected >= left);
                        assert(expected <= right);
                        current_key = keys[expected];
                        if (current_key < key) {
                            assert(left != expected + 1);
                            left = expected + 1;
                            key_left = keys[left];
                        }
                        else if (current_key > key) {
                            assert(right != expected - 1);
                            right = expected - 1;
                            key_right = keys[right];
                        }
                        else return std::optional<mapped_type>(values[expected]);
                    }
                    if(right == left and key == key_left)
                        return std::optional<mapped_type>(values[left]);
                    return std::nullopt;
                    break;
                }
            case Node_type::Tree:
                {
                    std::cerr << "InterpolationSearch for TreeNode not supported!" << std::endl;
                    break;
                }
            case Node_type::Hash:
                {
                    std::cerr << "InterpolationSearch for HashNode not supported!" << std::endl;
                    break;
                }
            default: break;
        }
        std::cerr << "InterpolationSearch strategy not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::vector<mapped_type> find_range(node_type *node, uint64_t &results, const key_type &lower, const key_type &upper) override {
        std::vector<mapped_type> res;
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    if (soa->empty()) return res;
                    size_type left = 0;
                    size_type n = soa->size();
                    size_type right = n - 1;
                    auto &keys = soa->keys();
                    key_type current_key = keys[right];
                    if (lower > current_key) return res; // no results
                    key_type key_left = keys[left];
                    key_type key_right = keys[right];
                    size_type expected = 0;
                    while (right > left and lower >= key_left and lower <= key_right) {
                        double slope = (double) (right - left) / (key_right - key_left);
                        expected = left + (lower - key_left) * slope;
                        assert(expected >= 0);
                        assert(expected < soa->size());
                        assert(expected >= left);
                        assert(expected <= right);
                        current_key = keys[expected];
                        if (current_key < lower) {
                            assert(left != expected + 1);
                            left = expected + 1;
                            key_left = keys[left];
                            expected = left;
                        }
                        else if (current_key > lower) {
                            assert(right != expected - 1);
                            /* If correct pos, break */
                            if (expected == 0 or lower > keys[expected - 1]) break;
                            right = expected - 1;
                            key_right = keys[right];
                            expected = right;
                        } else
                            break;
                    }
                    /* At this point, `expected` is our starting point for linear search until the keys do not match our
                     * range anymore. */
                    if constexpr (std::is_pointer_v<mapped_type>) {
                        auto &values = soa->values();
                        key_type x_1 = expected == 0 ? 0 : keys[expected - 1] + 1;
                        for (; expected < n; ++expected) {
                            key_type x_2 = keys[expected];
                            if (x_1 <= upper) {
                                if (lower <= x_2) {
                                    res.push_back(values[expected]);
                                    x_1 = x_2 + 1;
                                }
                            } else { return res; }
                        }
                    } else {
                        for (; expected < n; ++expected) {
                            /* We know that keys[i] >= lower due to while loop. */
                            key_type k = keys[expected];
                            if (k <= upper) ++results;
                            else break; // early abort
                        }
                    }
                    return res;
                }
            case Node_type::Tree:
                {
                    std::cerr << "InterploationSearch::find_range() for TreeNode not supported!" << std::endl;
                    return res;
                }
            case Node_type::Hash:
                {
                    std::cerr << "InterploationSearch::find_range() for HashNode not supported!" << std::endl;
                    return res;
                }
            default: return res;
        }
        std::cerr << "InterpolationSearch strategy not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::optional<mapped_type> choose_subset(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    if (soa->empty()) return std::nullopt;
                    size_type left = 0;
                    size_type right = soa->size() - 1;
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    Key current_key = keys[right];
                    if (key > current_key) return std::nullopt;
                    Key key_left = keys[left];
                    Key key_right = keys[right];
                    size_type expected = 0;
                    while (right > left and key >= key_left and key <= key_right) {
                        double slope = (double) (right - left) / (key_right - key_left);
                        expected = left + (key - key_left) * slope;
                        assert(expected >= 0);
                        assert(expected < soa->size());
                        assert(expected >= left);
                        assert(expected <= right);
                        current_key = keys[expected];
                        if (current_key < key) {
                            assert(left != expected + 1);
                            left = expected + 1;
                            key_left = keys[left];
                            expected = left;
                        }
                        else if (current_key > key) {
                            assert(right != expected - 1);
                            /* If correct pos, break */
                            if (expected == 0 or key > keys[expected - 1]) break;
                            right = expected - 1;
                            key_right = keys[right];
                            expected = right;
                        } else {
                            break;
                        }
                    }
                    return std::optional<mapped_type>(values[expected]);
                    break;
                }
            case Node_type::Tree:
                {
                    std::cerr << "InterpolationSearch::choose_subset for TreeNode not supported!" << std::endl;
                    break;
                }
            case Node_type::Hash:
                {
                    std::cerr << "InterpolationSearch::choose_subset for HashNode not supported!" << std::endl;
                    break;
                }
            default: break;
        }
        std::cerr << "InterpolationSearch::find_range() not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    IndexSearchMethod<key_type, mapped_type>* clone() override { return new InterpolationSearch<key_type, mapped_type>(); }
};

/*======================================================================================================================
 * ExponentialSearch
 *====================================================================================================================*/
template<typename Key, typename Value>
struct ExponentialSearch : IndexSearchMethod<Key,Value>
{
    using Base = IndexSearchMethod<Key, Value>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(node_type);
    IMPORT(size_type);

    ExponentialSearch() : Base(SearchType::ExponentialSearchMethod) {}

    void print() override { std::cout << "Exponential" << std::endl; }

    std::optional<mapped_type> lookup(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    size_type n = keys.size();
                    if (soa->empty()) return std::nullopt;
                    if (keys[0] == key) return std::optional<mapped_type>(values[0]);

                    size_type i = 1;
                    while (i < n and keys[i] <= key) i *= 2;

                    // Binary Search
                    auto start_it = keys.begin();
                    std::advance(start_it, i/2);
                    auto end = std::min(i, n);
                    auto end_it = keys.begin();
                    std::advance(end_it, end);
                    assert(start_it <= end_it);
                    auto entry = std::lower_bound(start_it, end_it, key);
                    if (entry != end_it and key == *entry) {
                        std::size_t idx = std::distance(keys.begin(), entry);
                        return std::optional<mapped_type>(values[idx]);
                    }
                    return std::nullopt;
                    break;
                }
            case Node_type::Tree:
                {
                    std::cerr << "ExponentialSearch::lookup() for TreeNode not supported!" << std::endl;
                    break;
                }
            case Node_type::Hash:
                {
                    std::cerr << "ExponentialSearch::lookup() for HashNode not supported!" << std::endl;
                    break;
                }
            default: break;
        }
        std::cerr << "ExponentialSearch::lookup() not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::vector<mapped_type> find_range(node_type *node, uint64_t &results, const key_type &lower, const key_type &upper) override {
        std::vector<mapped_type> res;
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    if (soa->empty()) return res;
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    size_type n = keys.size();
                    size_type i = 0;
                    if (lower <= keys[i]) {
                        if constexpr (std::is_pointer_v<mapped_type>)
                            res.push_back(values[0]);
                        else
                            ++results;
                    } else {
                        if (n == 1) return res;
                        ++i;
                        while (i < n and keys[i] <= lower) i *= 2;
                    }

                    /* Linear scannig of area. */
                    if constexpr (std::is_pointer_v<mapped_type>) {
                        auto key_it = keys.begin();
                        std::advance(key_it, i/2);   // starting point for linear scanning
                        while (key_it != keys.end() && *key_it < lower) ++key_it; // set `key_it` to first element def. in range
                        if (key_it == keys.end())
                            return res;
                        size_type x_1 = 0;
                        if (key_it != keys.begin()) x_1 = *(std::prev(key_it)) + 1;
                        auto value_it = values.begin();
                        size_type j = std::distance(keys.begin(), key_it);
                        std::advance(value_it, j);
                        for (auto end = keys.end(); key_it != end; ++key_it, ++value_it) {
                            if (x_1 > upper) return res; // early abort
                            res.push_back(*value_it);
                            x_1 = *key_it + 1;
                        }
                    } else {
                        size_type j = std::max<size_type>(1, i/2);
                        while (j < n && keys[j] < lower) ++j; // set `j` to first element definitely in range
                        if (j == n)
                            return res;
                        for (auto end = soa->size(); j < end; ++j) {
                            /* We know that keys[i] >= lower due to while loop. */
                            if (keys[j] <= upper) ++results;
                            else break; // early abort
                        }
                    }
                    return res;
                }
            case Node_type::Tree:
                {
                    std::cerr << "ExponentialSearch::find_range() for TreeNode not supported!" << std::endl;
                    return res;
                }
            case Node_type::Hash:
                {
                    std::cerr << "ExponentialSearch::find_range() for HashNode not supported!" << std::endl;
                    return res;
                }
            default: return res;
        }
        std::cerr << "ExponentialSearch::find_range() not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::optional<mapped_type> choose_subset(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    size_type n = keys.size();
                    if (soa->empty()) return std::nullopt;
                    if (key <= keys[0]) return std::optional<mapped_type>(values[0]);

                    size_type i = 1;
                    while (i < n and keys[i] <= key) i *= 2;

                    // Binary Search
                    auto start_it = keys.begin();
                    std::advance(start_it, i/2);
                    auto end = std::min(i + 1, n);
                    auto end_it = keys.begin();
                    std::advance(end_it, end);
                    assert(start_it <= end_it);
                    auto entry = std::lower_bound(start_it, end_it, key);
                    if (entry != end_it) {
                        std::size_t idx = std::distance(keys.begin(), entry);
                        return std::optional<mapped_type>(values[idx]);
                    }
                    return std::nullopt;
                    break;
                }
            case Node_type::Tree:
                {
                    std::cerr << "ExponentialSearch::choose_subset() for TreeNode not supported!" << std::endl;
                    break;
                }
            case Node_type::Hash:
                {
                    std::cerr << "ExponentialSearch::choose_subset() for HashNode not supported!" << std::endl;
                    break;
                }
            default: break;
        }
        std::cerr << "Exponential strategy for choose_subset not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    IndexSearchMethod<key_type, mapped_type>* clone() override { return new ExponentialSearch<key_type, mapped_type>(); }
};

/*======================================================================================================================
 * LinearRegression
 *====================================================================================================================*/
/** Computes the position of a given key inside a sorted array using linear regression.
 *  Formula: y = m * x + b, where
 *      `x`         : is the given key
 *      `y`         : is the (predicted) position inside the array
 *      `slope`     : is the slope of the linear function
 *      `intercept` : is the y intercept
 *  Additionally, we need error bounds to correct any "wrong" predictions.
 *  Based on: https://github.com/learnedsystems/RMI/blob/a0fa51858a0016d8b01f8e767f493255c510cd90/src/models/linear.rs
 */
template<typename Key, typename Value>
struct LinearRegressionSearch : IndexSearchMethod<Key,Value>
{
    using Base = IndexSearchMethod<Key, Value>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(node_type);
    IMPORT(size_type);

    private:
    double slope_; //< slope
    double intercept_; //< y intercept
    bool is_trained_ = false;

    size_type lo_;  //< highest overestimation, max. steps to the left until correct key is found
    size_type hi_; //< highest underestimation, max. steps to the right until correct key is found

    public:
    LinearRegressionSearch() : Base(SearchType::LinearRegressionSearchMethod) {}

    void print() override { std::cout << "LinearRegression" << std::endl; }
    double slope() const { return slope_; }
    double intercept() const { return intercept_; }
    bool is_trained() const { return is_trained_; }
    void untrain() { is_trained_ = false; }

    size_type predict(const key_type &k, const size_type &min, const size_type &max) const {
        return std::clamp<size_type>(slope_ * k + intercept_, min, max);
    }

    void train(node_type *node) {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();

                    /* Compute slope and intercept. */
                    if (keys.empty())       { slope_ = 0.0; intercept_ = 0.0; return; }
                    if (keys.size() == 1)   { slope_ = 0.0; intercept_ = 0.0; return; }

                    double mean_x = 0.0;
                    double mean_y = 0.0;
                    double c = 0.0;
                    double m2 = 0.0;
                    size_type n = 0;

                    for (std::size_t i = 0; i != keys.size(); ++i) {
                        n += 1;
                        key_type x = keys[i];
                        double y = i;

                        double dx = x - mean_x;
                        mean_x += dx /  n;
                        mean_y += (y - mean_y) / n;
                        c += dx * (y - mean_y);

                        double dx2 = x - mean_x;
                        m2 += dx * dx2;
                    }

                    double cov = c / (n - 1);
                    double var = m2 / (n - 1);
                    assert(var >= 0.0);

                    if (var == 0.0) { slope_ = 0.0; intercept_ = mean_y; }

                    slope_ = cov / var;
                    intercept_ = mean_y - slope_ * mean_x;

                    /* Compute error bounds. */
                    lo_ = 0;
                    hi_ = 0;
                    for (size_type i = 0; i < n; ++i) {
                        size_type pred = predict(keys[i], 0, n - 1);
                        size_type actual = i;
                        if (actual < pred)
                            lo_ = pred - actual > lo_ ? pred - actual : lo_;
                        else
                            hi_ = actual - pred > hi_ ? actual - pred : hi_;
                    }
                    is_trained_ = true;

                    return;
                    break;
                }
            case Node_type::Tree:
                {
                    std::cerr << "LinearRegressionSearch::train() for TreeNode not supported!" << std::endl;
                    break;
                }
            case Node_type::Hash:
                {
                    std::cerr << "LinearRegressionSearch::train() for HashNode not supported!" << std::endl;
                    break;
                }
            default: break;
        }
        std::cerr << "LinearRegressionSearch::train() not availabe for given node type." << std::endl;
    }

    std::optional<mapped_type> lookup(node_type *node, key_type &key) override {
        if (node->empty()) return std::nullopt;
        if (not is_trained()) {
            SearchMethodNotTrained ex;
            throw ex;
        }
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    size_type max_idx = keys.size() - 1;
                    size_type pred = predict(key, 0, max_idx);

                    /* Compute bounded area. */
                    size_type start = lo_ > pred ? 0 : pred - lo_;
                    size_type end = pred + hi_ > max_idx ? max_idx : pred + hi_;

                    /* Linear Search inside the error bounded area. */
                    for (; start <= end; ++start)
                        if (keys[start] == key) return std::optional<mapped_type>(values[start]);

                    return std::nullopt;
                    break;
                }
            case Node_type::Tree:
                {
                    std::cerr << "LinearRegressionSearch::lookup() for TreeNode not supported!" << std::endl;
                    break;
                }
            case Node_type::Hash:
                {
                    std::cerr << "LinearRegressionSearch::lookup() for HashNode not supported!" << std::endl;
                    break;
                }
            default: break;
        }
        std::cerr << "LinearRegressionSearch::lookup() not available for given node type." << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    std::optional<mapped_type> choose_subset(node_type *node, key_type &key) override {
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    size_type max_idx = keys.size() - 1;
                    size_type pred = predict(key, 0, max_idx);

                    /* Compute bounded area. */
                    size_type start = lo_ > pred ? 0 : pred - lo_;
                    size_type end = pred + hi_ > max_idx ? max_idx : pred + hi_;

                    /* Linear Search inside the error bounded area. */
                    for (; start <= end; ++start)
                        if (key <= keys[start]) return std::optional<mapped_type>(values[start]);

                    return std::nullopt;
                    break;
                }
            default: break;
        }
        std::cerr << "LinearRegression strategy for choose_subset not implemented!" << std::endl;
        return std::nullopt;
    }

    std::vector<mapped_type> find_range(node_type *node, uint64_t &results, const key_type &lower, const key_type &upper) override {
        std::vector<mapped_type> res;
        if (node->empty()) return res;
        if (not is_trained()) {
            SearchMethodNotTrained ex;
            throw ex;
        }
        switch (node->node_type()) {
            case Node_type::SoA:
                {
                    auto soa = static_cast<SoA<key_type, mapped_type>*>(node);
                    auto &keys = soa->keys();
                    auto &values = soa->values();
                    size_type max_idx = keys.size() - 1;
                    size_type pred = predict(lower, 0, max_idx);
                    if (lower > keys[max_idx]) return res; // range is out of scope

                    /* Compute bounded area containing `lower`. */
                    size_type start = lo_ > pred ? 0 : pred - lo_;
                    size_type end = pred + hi_ > max_idx ? max_idx : pred + hi_;

                    /* Set `start` to the very first element satisfying the range predicate. */
                    while (lower > keys[start] and start < end) ++start;
                    /* Linear Search. */
                    if constexpr (std::is_pointer_v<mapped_type>) {
                        auto key_it = keys.begin();
                        std::advance(key_it, start);   // starting point for linear scanning
                        size_type x_1 = 0;
                        if (key_it != keys.begin()) x_1 = *(std::prev(key_it)) + 1;
                        auto value_it = values.begin();
                        size_type j = std::distance(keys.begin(), key_it);
                        std::advance(value_it, j);
                        for (auto end = keys.end(); key_it != end; ++key_it, ++value_it) {
                            if (x_1 > upper) return res; // early abort
                            res.push_back(*value_it);
                            x_1 = *key_it + 1;
                        }
                    } else {
                        for (auto n = keys.size(); start < n; ++start) {
                            /* We know that keys[start] >= lower due to while loop. */
                            if (keys[start] <= upper) ++results;
                            else break; // early abort
                        }
                    }

                    return res;
                }
            default: return res;
        }
        std::cerr << "LinearRegression strategy for find_range not implemented!" << std::endl;
        SearchMethodNotAvailable ex;
        throw ex;
    }

    IndexSearchMethod<key_type, mapped_type>* clone() override { return new LinearRegressionSearch<key_type, mapped_type>(); }
};
