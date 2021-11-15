#pragma once

#include "hp/Partition.hpp"
#include <math.h>
#include <random>

using size_type = std::size_t;

/*======================================================================================================================
 * Bulkload B+Tree index
 *====================================================================================================================*/
template<typename Key, typename Value, typename Iterator>
std::unique_ptr<HorizontalPartition<Key, Value>> bulkload_btree(Iterator begin, Iterator end, size_type inner_slot_size, size_type leaf_slot_size)
{
    using key_type = Key;
    using mapped_type = Value;
    using hp_type = HorizontalPartition<key_type, mapped_type>;

    size_type num_elements = std::distance(begin, end);
    if (num_elements < 1) throw std::invalid_argument("Invalid iterators passed to bulkloading");

    /* Create leaves. */
    std::vector<hp_type*> current_level;
    hp_type* leaf = new SortedArray<key_type, mapped_type>("SoA", 0, leaf_slot_size);
    for (auto it = begin; it != end; ++it) {
        if(not leaf->insert(it->first, it->second, false).first) 
        {
            /* Leaf is full, add to current_level and create new one. */
            current_level.push_back(leaf);
            leaf = new SortedArray<key_type, mapped_type>("SoA", 0, leaf_slot_size);
            leaf->insert(it->first, it->second);
        }
    }
    current_level.push_back(leaf); // Add rightmost leaf to current_level
    size_type num_leaves = std::distance(current_level.begin(), current_level.end());

    if (num_leaves == 1) return std::unique_ptr<hp_type>(current_level[0]); // we only have one leaf, terminate

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
            hp_type *inode = new SortedArray<key_type, mapped_type>("SoA", inner_slot_size, 0);
            for (; chunk_begin != chunk_end; ++chunk_begin) {
                hp_type *current_node = *chunk_begin;
                key_type key = current_node->max_key();
                inode->insert(key, current_node, false);
            }
            next_level.push_back(inode);
            chunk_begin = chunk_end;
        } while (std::distance(chunk_begin, current_level.end()) > 0);
        current_level.swap(next_level);
    }
    assert(current_level.size() == 1);

    return std::unique_ptr<hp_type>(current_level[0]);
}


/*======================================================================================================================
 * Bulkload B+Tree like index with random node types
 *====================================================================================================================*/
template<typename Key, typename Value>
HorizontalPartition<Key, Value>* generate_random_partition(std::discrete_distribution<int>* dist, std::mt19937* gen, size_type partition_slot_size, size_type entry_slot_size) {
    HorizontalPartition<Key, Value>* hp;
    int choice = (*dist)(*gen);
    switch(choice) {
    case 0:
        hp = new SortedArray<Key, Value>("SoA", partition_slot_size, entry_slot_size);
        break;
    case 1:
        hp = new Tree<Key, Value>(partition_slot_size, entry_slot_size);
        break;
    case 2:
        hp = new Hash<Key, Value>(partition_slot_size, entry_slot_size);
        break;
    default:
        throw std::invalid_argument("Distribution should only return values between 0 and 2");
    }
    auto methods = hp->getPartitionSearchMethods();
    auto it = std::find(methods->begin(), methods->end(), "LinearRegressionSearch");
    if (it != methods->end() && methods->size() > 1) {
    *it = methods->back();
    methods->pop_back();
    }
    it = std::find(methods->begin(), methods->end(), "LinearSearch");
    if (it != methods->end() && methods->size() > 1) {
    *it = methods->back();
    methods->pop_back();
    }
    std::uniform_int_distribution<unsigned int> dPartition(0, methods->size()-1);
    hp->set_partition_search_method((*methods)[dPartition(*gen)]);
    delete(methods);
    methods = hp->getEntrySearchMethods();
    it = std::find(methods->begin(), methods->end(), "LinearRegressionSearch");
    if (it != methods->end() && methods->size() > 1) {
    *it = methods->back();
    methods->pop_back();
    }
    it = std::find(methods->begin(), methods->end(), "LinearSearch");
    if (it != methods->end() && methods->size() > 1) {
    *it = methods->back();
    methods->pop_back();
    }
    std::uniform_int_distribution<unsigned int> dEntry(0, methods->size()-1);
    hp->set_entry_search_method((*methods)[dEntry(*gen)]);
    delete(methods);
    return hp;
}

template<typename Key, typename Value, typename Iterator>
HorizontalPartition<Key, Value>* bulkload_random_tree(Iterator begin, Iterator end, size_type partition_slot_size_leaves, size_type entry_slot_size_leaves, size_type partition_slot_size_inner, size_type entry_slot_size_inner, 
        double load_factor_entries_leaves=1.0, double load_factor_partitions_inner=1.0)
{
    using key_type = Key;
    using mapped_type = Value;
    using hp_type = HorizontalPartition<key_type, mapped_type>;

    size_type num_elements = std::distance(begin, end);
    if (num_elements < 1) throw std::invalid_argument("Invalid iterators passed to bulkloading");

    /* Create leaves. */
    std::mt19937* gen = new std::mt19937(std::random_device{}());
    std::discrete_distribution<int>* dist = new std::discrete_distribution<int>{1.0, 1.0, 1.0};
    std::vector<hp_type*> current_level;
    hp_type* leaf = generate_random_partition<key_type, mapped_type>(dist, gen, partition_slot_size_leaves, entry_slot_size_leaves);
    auto max_entries = ceil(load_factor_entries_leaves * leaf->entry_capacity());
    for (auto it = begin; it != end; ++it) {
        if (leaf->get_entries()->size() >= max_entries) {
            current_level.push_back(leaf);
            leaf = generate_random_partition<key_type, mapped_type>(dist, gen, partition_slot_size_leaves, entry_slot_size_leaves);
        }
        leaf->insert(it->first, it->second);
    }
    current_level.push_back(leaf);
    size_type num_leaves = std::distance(current_level.begin(), current_level.end());
    delete(dist);

    if (num_leaves == 1) {
        delete(gen);
        return current_level[0];
    }

    dist = new std::discrete_distribution<int>{1.0, 1.0, 0.0};
    long chunk_size = load_factor_partitions_inner * partition_slot_size_inner;
    while (current_level.size() > 1) {
        std::vector<hp_type*> next_level;
        auto chunk_begin = current_level.begin();
        auto chunk_end = current_level.begin();
        do {
            if (std::distance(chunk_end, current_level.end()) < chunk_size) {
                chunk_end = current_level.end();
            } else {
                std::advance(chunk_end, chunk_size);
            }
            hp_type* inode = generate_random_partition<key_type, mapped_type>(dist, gen, partition_slot_size_inner, entry_slot_size_inner);
            for (; chunk_begin != chunk_end; ++chunk_begin) {
                hp_type* current_node = *chunk_begin;
                key_type key = current_node->max_key();
                inode->insert(key, current_node, false);
            }
            next_level.push_back(inode);
            chunk_begin = chunk_end;
        } while (std::distance(chunk_begin, current_level.end()) > 0);
        current_level.swap(next_level);
    }
    assert(current_level.size() == 1);
    delete(gen);
    delete(dist);
    return current_level[0];
}

template<typename Key, typename Value, typename Iterator>
HorizontalPartition<Key, Value>* bulkload_soa_tree(Iterator begin, Iterator end, size_type partition_slot_size_leaves, size_type entry_slot_size_leaves, size_type partition_slot_size_inner, size_type entry_slot_size_inner, 
        double load_factor_entries_leaves=1.0, double load_factor_partitions_inner=1.0)
{
    using key_type = Key;
    using mapped_type = Value;
    using hp_type = HorizontalPartition<key_type, mapped_type>;

    size_type num_elements = std::distance(begin, end);
    if (num_elements < 1) throw std::invalid_argument("Invalid iterators passed to bulkloading");

    /* Create leaves. */
    std::vector<hp_type*> current_level;
    hp_type* leaf = new SortedArray<key_type, mapped_type>("SoA", partition_slot_size_leaves, entry_slot_size_leaves);
    auto max_entries = ceil(load_factor_entries_leaves * leaf->entry_capacity());
    for (auto it = begin; it != end; ++it) {
        if (leaf->get_entries()->size() >= max_entries) {
            current_level.push_back(leaf);
            leaf = new SortedArray<key_type, mapped_type>("SoA", partition_slot_size_leaves, entry_slot_size_leaves);
        }
        leaf->insert(it->first, it->second);
    }
    current_level.push_back(leaf);
    size_type num_leaves = std::distance(current_level.begin(), current_level.end());

    if (num_leaves == 1) {
        return current_level[0];
    }

    long chunk_size = load_factor_partitions_inner * partition_slot_size_inner;
    while (current_level.size() > 1) {
        std::vector<hp_type*> next_level;
        auto chunk_begin = current_level.begin();
        auto chunk_end = current_level.begin();
        do {
            if (std::distance(chunk_end, current_level.end()) < chunk_size) {
                chunk_end = current_level.end();
            } else {
                std::advance(chunk_end, chunk_size);
            }
            hp_type* inode = new SortedArray<key_type, mapped_type>("SoA", partition_slot_size_inner, entry_slot_size_inner);
            for (; chunk_begin != chunk_end; ++chunk_begin) {
                hp_type* current_node = *chunk_begin;
                key_type key = current_node->max_key();
                inode->insert(key, current_node, false);
            }
            next_level.push_back(inode);
            chunk_begin = chunk_end;
        } while (std::distance(chunk_begin, current_level.end()) > 0);
        current_level.swap(next_level);
    }
    assert(current_level.size() == 1);
    return current_level[0];
}

/* works only for uni-dense datasets with an integer scale up factor */
template<typename Key, typename Value>
HorizontalPartition<Key, Value>* bulkload_existing_index(HorizontalPartition<Key, Value>* hp, unsigned int scale_factor) {
    using key_type = Key;
    using mapped_type = Value;
    using hp_type = HorizontalPartition<key_type, mapped_type>;

    /* Create empty, cloned index structure */
    auto clone = hp->clone();
    clone->clearEntries();
    
    /* Increase capacity and generate entries */
    std::vector<hp_type*> open_list;
    open_list.push_back(clone);
    while(not open_list.empty()) {
        auto node = open_list.back();
        open_list.pop_back();
        hp_type* original_node;
        if(node->getId() == clone->getId())
            original_node = hp;
        else {
            original_node = hp->lookupPartitionKey(node->max_key(), node->getId()).value();
        }
        node->get_partitions()->set_capacity(node->partition_capacity() * scale_factor);
        node->get_entries()->set_capacity(node->entry_capacity() * scale_factor);
        auto entry_keys = original_node->getEntryKeys(false);
        for (auto e : *entry_keys) {
            for(unsigned int i=0; i<scale_factor; i++) {
                auto k = e * scale_factor + i;
                auto v = original_node->pointQ(e).value() * scale_factor + i;
                node->insert(k, v);
            }
        }
        delete(entry_keys);
        /* fix partition keys */
        if(node->get_partitions()->size() > 0) {
            auto partition_keys = node->getPartitionKeys(false);
            std::sort(partition_keys->begin(), partition_keys->end());
            for (unsigned int x=0; x<partition_keys->size(); x++) {
                auto p = (*partition_keys)[partition_keys->size() - x - 1];
                auto child_partition = node->lookupPartitionKey(p).value();
                node->remove_partition(p);
                node->insert((p + 1) * scale_factor - 1, child_partition, false);
                open_list.push_back(child_partition);
            }
            delete(partition_keys);
        }
    }
    if (not clone->isValid()) {
        throw std::invalid_argument("Invalid Index after bulkload");
    }
    return clone;
}

/**
 * Scaling of an existing index up to a new dataset
 * Upper layers of the index remain the same
 * only exception: largest pivots might be changed to accomodate for keys larger than those in the original tree
 * leaves are replicated to store additional keys within a specific interval (defined by the pivot elements of the original tree)
 */
template<typename Key, typename Value, typename Iterator>
Iterator scale_existing_index_(HorizontalPartition<Key, Value>* hp, Iterator begin, Iterator end, Key pivot, bool replicate) {
    if (hp->partition_size() > 0) {
        /* inner node */
        auto p_keys = hp->getPartitionKeys(false);
        std::sort(p_keys->begin(), p_keys->end());
        for(auto k: *p_keys) {
            auto child = hp->lookupPartitionKey(k).value();
            if (child->get_partitions()->size() > 0) {
                /* child is inner node */
                begin = scale_existing_index_<Key, Value>(child, begin, end, k, replicate);
            } else {
                /* child is leaf node */
                begin = scale_existing_index_<Key, Value>(child, begin, end, k, replicate);
                if (begin != end and begin->first <= k) {
                    hp->remove_partition(k);
                    hp->insert(child->max_key(), child);
                }
                while(begin != end and begin->first <= k) {
                    auto new_child = child->clone();
                    begin = scale_existing_index_<Key, Value>(new_child, begin, end, k, replicate);
                    hp->set_partition_capacity(hp->partition_capacity() + 1);
                    if(begin != end and begin->first <= k)
                        hp->insert(new_child->max_key(), new_child);
                    else    
                        hp->insert(k, new_child);
                }
            }
        }
        delete(p_keys);
    } else {
        /* leaf node */
        unsigned int current_size = hp->entry_size();
        hp->clearEntries();
        for(; begin != end and begin->first <= pivot and (not replicate or hp->entry_size() < current_size); ++begin) {
            if(not replicate and hp->entry_size() >= hp->entry_capacity())
                hp->set_entry_capacity(hp->entry_capacity() + 1);
            hp->insert(begin->first, begin->second);
        }
    }
    return begin;
}

template<typename Key, typename Value, typename Iterator>
HorizontalPartition<Key, Value>* scale_existing_index(HorizontalPartition<Key, Value>* hp, Iterator begin, Iterator end, bool replicate=true) {
    /* Create empty, cloned index structure */
    auto clone = hp->clone();
    
    if(clone->get_partitions()->size() > 0) {
        /* multi-level index */
        begin = scale_existing_index_<Key, Value>(clone, begin, end, Key(), replicate);
        /* find rightmost leaf */
        auto current_node = clone;
        while(current_node->get_partitions()->size() > 0) {
            auto p_keys = current_node->getPartitionKeys(false);
            std::sort(p_keys->begin(), p_keys->end());
            current_node = current_node->lookupPartitionKey((*p_keys)[p_keys->size()-1]).value();
            delete(p_keys);
        }
        for (; begin != end; ++begin) {
            current_node->set_entry_capacity(current_node->entry_capacity() + 1);
            current_node->insert(begin->first, begin->second);
        }
    } else {
        /* single-level index */
        clone->clearEntries();
        clone->set_entry_capacity(std::distance(begin, end));
        for (; begin != end; ++begin) {
            clone->insert(begin->first, begin->second);
        }
    }
    return clone;
}
