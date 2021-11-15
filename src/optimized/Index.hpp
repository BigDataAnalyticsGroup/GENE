#pragma once

#include "optimized/Partition.hpp"
#include "util/hash.hpp"
#include <cassert>
#include <iostream>

namespace opt {

template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct Index final
{
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<const key_type, T>;
    using key_compare = Compare;
    using size_type = std::size_t;
    using hp_type = HorizontalPartition<key_type, mapped_type, key_compare>;

    private:
    hp_type * root_;
    std::vector<value_type> hold_out_entries_;

    public:
    ~Index() { delete root_; }

    /* Factory methods. */
    template<typename It>
    static Index Bulkload(It begin, It end, size_type inner_slot_size, size_type leaf_slot_size);

    /* Data access. */
    /** Returns the value corresponding to the given `key`, empty optional otherwise. */
    std::optional<mapped_type> find(const key_type &key) const;
    /** Returns the value corresponding to the lower bound for the given `key`, empty optional otherwise. */
    std::optional<mapped_type> lower_bound(const key_type &key) const;
    /** Returns the number of qualifying keys in the index, i.e. keys in the range [lower, upper]. */
    size_type find(const key_type &lower, const key_type &upper) const;

    /* Visualization. */
    void dot(std::ostream &out, bool structure_only = false) const;
};

/*======================================================================================================================
 * Implementation
 *====================================================================================================================*/
template<typename Key, typename T, typename Compare>
std::optional<T> Index<Key, T, Compare>::find(const key_type &key) const
{
    hp_type *node = root_;
    while (not node->is_leaf()) {
        switch (node->partition_type) {
            case PARTITION_type::InnerSortedBinary:
            {
                const InnerSortedBinary<key_type, mapped_type> *inode = static_cast<const InnerSortedBinary<key_type, mapped_type>*>(node);
                auto child_node = inode->find(key);
                if (child_node) node = *child_node;
                else return std::nullopt;
                break;
            }
            case PARTITION_type::InnerSortedLinear:
            {
                const InnerSortedLinear<key_type, mapped_type> *inode = static_cast<const InnerSortedLinear<key_type, mapped_type>*>(node);
                auto child_node = inode->find(key);
                if (child_node) node = *child_node;
                else return std::nullopt;
                break;
            }
            default:
            {
                std::cerr << "No matching innner partition." << std::endl;
                return std::nullopt;
            }
        }
    }
    switch (node->partition_type) {
        case PARTITION_type::LeafSortedBinary:
        {
            const LeafSortedBinary<key_type, mapped_type> *leaf = static_cast<const LeafSortedBinary<key_type, mapped_type>*>(node);
            return leaf->find(key);
        }
        case PARTITION_type::LeafSortedInterpolation:
        {
            const LeafSortedInterpolation<key_type, mapped_type> *leaf = static_cast<const LeafSortedInterpolation<key_type, mapped_type>*>(node);
            return leaf->find(key);
        }
        case PARTITION_type::LeafHashUnorderedMap:
        {
            const LeafHashUnorderedMap<key_type, mapped_type> *leaf = static_cast<const LeafHashUnorderedMap<key_type, mapped_type>*>(node);
            return leaf->find(key);
        }
        case PARTITION_type::LeafRobinHoodHashTable:
        {
            const LeafRobinHoodHashTable<key_type, mapped_type> *leaf = static_cast<const LeafRobinHoodHashTable<key_type, mapped_type>*>(node);
            return leaf->find(key);
        }
        default:
        {
            std::cerr << "No matching leaf partition." << std::endl;
            return std::nullopt;
        }
    }
}
template<typename Key, typename T, typename Compare>
std::optional<T> Index<Key, T, Compare>::lower_bound(const key_type &key) const
{
    hp_type *node = root_;
    while (not node->is_leaf()) {
        switch (node->partition_type) {
            case PARTITION_type::InnerSortedBinary:
            {
                const InnerSortedBinary<key_type, mapped_type> *inode = static_cast<const InnerSortedBinary<key_type, mapped_type>*>(node);
                auto child_node = inode->find(key);
                if (child_node) node = *child_node;
                else return std::nullopt;
                break;
            }
            case PARTITION_type::InnerSortedLinear:
            {
                const InnerSortedLinear<key_type, mapped_type> *inode = static_cast<const InnerSortedLinear<key_type, mapped_type>*>(node);
                auto child_node = inode->find(key);
                if (child_node) node = *child_node;
                else return std::nullopt;
                break;
            }
            default:
            {
                std::cerr << "No matching innner partition." << std::endl;
                return std::nullopt;
            }
        }
    }
    switch (node->partition_type) {
        case PARTITION_type::LeafSortedBinary:
        {
            const LeafSortedBinary<key_type, mapped_type> *leaf = static_cast<const LeafSortedBinary<key_type, mapped_type>*>(node);
            return leaf->lower_bound(key);
        }
        case PARTITION_type::LeafSortedInterpolation:
        {
            const LeafSortedInterpolation<key_type, mapped_type> *leaf = static_cast<const LeafSortedInterpolation<key_type, mapped_type>*>(node);
            return leaf->lower_bound(key);
        }
        default:
        {
            std::cerr << "lower_bound: No matching leaf partition." << std::endl;
            return std::nullopt;
        }
    }
}

template<typename Key, typename T, typename Compare>
void Index<Key, T, Compare>::dot(std::ostream &out, bool structure_only) const
{
    out << "digraph index {\n"
        << "node [shape=plain, height=.1];\n";
    root_->dot(out, "0", structure_only);
    out << "}\n";
}

/*======================================================================================================================
 * Bulkload B+Tree index
 *====================================================================================================================*/
template<typename Key, typename T, typename Compare>
template<typename It>
Index<Key, T, Compare> Index<Key, T, Compare>::Bulkload(It begin, It end,
                                                        size_type inner_slot_size, size_type leaf_slot_size)
{
    Index idx;

    size_type num_elements = std::distance(begin, end);
    if (num_elements < 1) throw std::invalid_argument("Invalid iterators passed to bulkloading");

    /* Manually drawn cutoffs. */
    std::vector<double> cutoffs = {0.1, 0.85};

    /* Compute capacities for each leaf. */
    std::vector<size_type> capacities;
    for (auto i = 0; i < cutoffs.size(); ++i) {
        if (i == 0) { capacities.push_back(cutoffs[i] * num_elements); continue; }
        capacities.push_back(cutoffs[i] * num_elements - (cutoffs[i - 1] * num_elements));
    }
    capacities.push_back(num_elements - (cutoffs.back() * num_elements));


    auto it = begin; // data iterator
    /* Create root node. */
    InnerSortedBinary<key_type, mapped_type> *inode = new InnerSortedBinary<key_type, mapped_type>(3);

    /* Create left leaf. */
    LeafPartition<key_type, mapped_type> *left_leaf = new LeafRobinHoodHashTable<key_type, mapped_type>(capacities[0]);
    for (auto i = 0; i < left_leaf->capacity(); ++i, ++it) left_leaf->insert(it->first, it->second);
    inode->insert(std::prev(it)->first, left_leaf);

    /* Create btree 'leaf' for middle partition. */
    using inner_type = InnerSortedBinary<key_type, mapped_type>;
    using leaf_type = LeafSortedBinary<key_type, mapped_type>;

    std::vector<hp_type*> current_level;
    LeafPartition<key_type, mapped_type> *leaf = new leaf_type(leaf_slot_size);
    for (auto i = 0; i < capacities[1]; ++i, ++it) {
        if (not leaf->insert(it->first, it->second)) {
            /* Leaf is full, add to current_level and create new one. */
            current_level.push_back(leaf);
            auto new_leaf = new leaf_type(leaf_slot_size);
            leaf->set_next(new_leaf);
            leaf = new_leaf;
            leaf->insert(it->first, it->second);
        }
    }
    current_level.push_back(leaf); // Add rightmost leaf to current_level
    size_type num_leaves = std::distance(current_level.begin(), current_level.end());

    if (num_leaves > 1) {
        long chunk_size = inner_slot_size;
        while (current_level.size() > 1) {
            std::vector<hp_type*> next_level;
            /* Process current_level in chunks of size `inner_slot_size`. For each chunk create new inode, add the
             * keys/children and insert the new node into the next level. */
            auto chunk_begin = current_level.begin();
            auto chunk_end = current_level.begin();
            do {
                if (std::distance(chunk_end, current_level.end()) < chunk_size) {  // we are at the end
                    chunk_end = current_level.end();
                } else {
                    std::advance(chunk_end, inner_slot_size);
                }
                /* Process chunk. */
                InnerSortedBinary<key_type, mapped_type> *inode = new inner_type(inner_slot_size);
                for (; chunk_begin != chunk_end; ++chunk_begin) {
                    hp_type *current_node = *chunk_begin;
                    key_type key = current_node->max_key();
                    inode->insert(key, current_node);
                }
                next_level.push_back(inode);
                chunk_begin = chunk_end;
            } while (std::distance(chunk_begin, current_level.end()) > 0);
            current_level.swap(next_level);
        }
    }
    assert(current_level.size() == 1);
    inode->insert(std::prev(it)->first, current_level[0]);

    /* Create right leaf. */
    LeafPartition<key_type, mapped_type> *right_leaf = new LeafRobinHoodHashTable<key_type, mapped_type>(capacities[2]);
    for (auto i = 0; i < right_leaf->capacity(); ++i, ++it) right_leaf->insert(it->first, it->second);
    inode->insert(std::prev(it)->first, right_leaf);

    assert(it == end); // all elements inserted
    idx.root_ = inode;
    return idx;
}

}
