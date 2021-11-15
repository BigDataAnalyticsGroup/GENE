#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <fstream>

#include "hp/IndexSearchMethods.hpp"
#include "hp/NodeStatistics.hpp"
#include "nodes/Node.hpp"
#include "util/enum.hpp"
#include "util/Exceptions.hpp"
#include "util/globals.hpp"

#define IMPORT(WHAT) using WHAT = typename Base::WHAT

template<typename key_type, typename mapped_type> struct SortedArray;
template<typename key_type, typename mapped_type> struct Tree;
template<typename key_type, typename mapped_type> struct Hash;

/*======================================================================================================================
 * HorizontalPartition
 *====================================================================================================================*/
template<typename Key, typename Value>
struct HorizontalPartition
{
    using size_type = std::size_t;
    using key_type = Key;
    using mapped_type = Value;
    using hp_type = HorizontalPartition<key_type, mapped_type>;
    using id_type = unsigned long;

    private:
    static id_type ID;
    id_type myId;

    protected:
    Node<key_type, hp_type*> *partitions_; //< maps keys to nested horizontal partitions
    Node<key_type, mapped_type> *entries_; //< maps keys to entries

    IndexSearchMethod<key_type, hp_type*> *partition_search_; //< search algorithm to find partitions
    IndexSearchMethod<key_type, mapped_type> *entry_search_; //< search algorithm to find entries

    std::unordered_map<std::string, std::unique_ptr<IndexSearchMethod<key_type, hp_type*>>> partition_search_methods_; //< set of partition search algorithms
    std::unordered_map<std::string, std::unique_ptr<IndexSearchMethod<key_type, mapped_type>>> entry_search_methods_; //< set of entry search algorithms

    key_type max_key_ = std::numeric_limits<key_type>::min();
    key_type min_key_ = std::numeric_limits<key_type>::max();

    PartitionType partitionType = PartitionType::UnknownPartition;

    NodeStatistics<Key, Value>* stats_;

    HorizontalPartition()
    : myId(ID++)
    {
        stats_ = new NodeStatistics<Key, Value>(PartitionType::UnknownPartition);
    }

    public:
    virtual ~HorizontalPartition() {
        for (auto it = partitions_->value_begin(), end = partitions_->value_end(); it != end; ++it) { delete *it; }
        delete partitions_;
        delete entries_;
        delete(stats_);
    }


    /*==================================================================================================================
     * Modifiers
     *================================================================================================================*/
    /** Try adding key `k` and value `v` locally to entries.  If maxiumum capacity is reached, insert recursively to
     * nested horizontal partitions, if any. */
    virtual std::pair<bool, std::optional<HorizontalPartition<key_type, mapped_type>*>> insert(key_type k, mapped_type v, bool split=true) {
        auto result = insert_(k, v, split);
        if (result.second) {
            auto new_root = new SortedArray<key_type, mapped_type>("SoA");
            new_root->insert(result.second.value()->max_key_, result.second.value(), false);
            new_root->insert(max_key_, this, false);
            return std::make_pair(result.first, new_root);
        } else {
            return std::make_pair(result.first, std::nullopt);
        }
    }
    
    /**
     * Insert a given key and value into the partition, descending recursively if necessary
     * If a leaf does not have enough space for the new entry, split the leaf to create more space
     * Returns a boolean indicator of success as well as an optional HP resulting from a split which has to be inserted in the parent
     */
    virtual std::pair<bool, std::optional<HorizontalPartition<key_type, mapped_type>*>> insert_(key_type k, mapped_type v, bool split=true) {
        this->stats_->add_execution(0, Workload_kind::Insert);
        if (entries_->insert(k, v)) { // may raise `DuplicateKeyException`
            /* Local insert successful. */
            try {
                auto lr_search = dynamic_cast<LinearRegressionSearch<key_type, mapped_type>*>(entry_search_methods_.at("LinearRegressionSearch").get());
                lr_search->untrain();
            } catch (std::out_of_range &e) { /* do nothing */ }

            if (k > max_key_) max_key_ = k;
            if (k < min_key_) min_key_ = k;
            return std::make_pair(true, std::nullopt);
        } else {
            /* Local insert failed, insert recursively. */
            std::pair<bool, std::optional<HorizontalPartition<key_type, mapped_type>*>> result;
            if (not partitions_->empty()) {
                /* Inner node */
                auto hp = choose_partition(k); // find corresponding partition
                if (not hp) { // no valid partition exists
                    auto partition_keys = getPartitionKeys(false);
                    std::sort(partition_keys->begin(), partition_keys->end());
                    key_type largest_key = (*partition_keys)[partition_keys->size()-1];
                    assert(largest_key < k);
                    auto largest_partition = lookupPartitionKey(largest_key).value();
                    remove_partition(largest_key);
                    insert(k, largest_partition, false);
                    delete(partition_keys);
                    result = largest_partition->insert_(k, v, split);
                } else {
                    result = (*hp)->insert_(k, v, split); // continue to traverse horizontal partitions
                }
                if (result.second) {
                    // child hp was split, insert newly created hp in this hp
                    if (partitions_->size() >= partitions_->capacity()) {
                        auto new_partition = split_partition();
                        if (result.second.value()->max_key_ > new_partition->max_key_)
                            insert(result.second.value()->max_key_, result.second.value(), false);
                        else
                            new_partition->insert(result.second.value()->max_key_, result.second.value(), false);
                        return std::make_pair(result.first, new_partition);
                    } else {
                        insert(result.second.value()->max_key_, result.second.value(), false);
                        return std::make_pair(result.first, std::nullopt);
                    }
                } else {
                    return result;
                }            
            } else {
                /* Leaf node */
                if (split) {
                    auto new_partition = split_partition();
                    if (k > new_partition->max_key()) {
                        insert_(k, v, false);
                    } else {
                        new_partition->insert_(k, v, false);
                    }
                    return std::make_pair(true, new_partition);
                } else {
                    stats_->add_penalty(Configuration::Get().penaltyFactorMissingKeys); // penalty
                    return std::make_pair(false, std::nullopt);    
                }
            }
        }
    }

    virtual HorizontalPartition<key_type, mapped_type>* split_partition() {
        /* Creating new neighbor partition with same configuration */
        this->stats_->add_split();
        HorizontalPartition<key_type, mapped_type>* neighbor;
        switch(partitionType) {
            case PartitionType::SoAPartition:
                neighbor = new SortedArray<key_type, mapped_type>("SoA", partitions_->capacity(), entries_->capacity());
                break;
            case PartitionType::TreePartition:
                neighbor = new Tree<key_type, mapped_type>(partitions_->capacity(), entries_->capacity());
                break;
            case PartitionType::HashPartition:
                neighbor = new Hash<key_type, mapped_type>(partitions_->capacity(), entries_->capacity());
                break;
            default:
                throw std::invalid_argument("not implemented yet!");
        }
        neighbor->set_partition_search_method(toString(partition_search_->getSearchType()));
        neighbor->set_entry_search_method(toString(entry_search_->getSearchType()));
        if (partitions_->capacity() == 0 and entries_->capacity() == 0)
            throw std::invalid_argument("Partition with zero capacities makes no sense");
        if (partitions_->size() == 0) {
            /* Splitting leaf node */
            std::vector<key_type>* entry_keys = getEntryKeys(false);
            std::sort(entry_keys->begin(), entry_keys->end());
            unsigned int end_index = (unsigned int)(0.5 * entries_->size());
            /* Moving the lower half of the entries to the neighbor */
            for (unsigned int i=0; i<end_index; i++) {
                key_type entry_to_move = (*entry_keys)[i];
                mapped_type value_to_move = pointQ(entry_to_move).value();
                remove_entry(entry_to_move);
                neighbor->insert_(entry_to_move, value_to_move, false);
            }
            delete(entry_keys);
            return neighbor;
        } else {
            /* Splitting inner node */
            std::vector<key_type>* partition_keys = getPartitionKeys(false);
            std::sort(partition_keys->begin(), partition_keys->end());
            std::vector<key_type>* entry_keys = getEntryKeys(false);
            std::sort(entry_keys->begin(), entry_keys->end());
            unsigned int end_index = (unsigned int)(0.5 * partitions_->size());
            key_type largest_key = std::numeric_limits<key_type>::min();
            /* Moving the lower half of the partitions to the neighbor */
            for (unsigned int i=0; i<end_index; i++) {
                key_type key_to_move = (*partition_keys)[i];
                if (key_to_move > largest_key)
                    largest_key = key_to_move;
                hp_type* partition_to_move = lookupPartitionKey(key_to_move).value();
                remove_partition(key_to_move);
                neighbor->insert(key_to_move, partition_to_move, false);
            }
            /* Moving all entries with lower or equal keys than `largest_key` to the neighbor */
            for (unsigned int i=0; i<entry_keys->size(); i++) {
                key_type entry_to_move = (*entry_keys)[i];
                if (entry_to_move > largest_key)
                    break;
                mapped_type value_to_move = pointQ(entry_to_move).value();
                remove_entry(entry_to_move);
                neighbor->insert_(entry_to_move, value_to_move, false);
            }
            delete(partition_keys);
            delete(entry_keys);
            return neighbor;
        }
    }

    /** Add key `k` and horizontal partition `hp` to `partitions`. */
    virtual bool insert(key_type, hp_type *, bool = true);

    virtual bool remove_entry(key_type k) {
        if (not entries_->remove(k)) {
            /* Local remove failed. */
            if (not partitions_->empty()) {
                auto hp = choose_partition(k); // find corresponding partition
                if (not hp) { // no valid partition exists
                    /* No key to delete. */
                    return false;
                } else return (*hp)->remove_entry(k); // recursively delete k
            } else return false;
        }
        return true;
    }
    virtual bool remove_partition(key_type k) {
        if (not partitions_->remove(k)) {
            /* Local remove failed. */
            if (not partitions_->empty()) {
                auto hp = choose_partition(k); // find corresponding partition
                if (not hp) { // no valid partition exists
                    /* No key to delete. */
                    return false;
                } else return (*hp)->remove_partition(k); // recursively delete k
            } else return false;
        }
        return true;
    }

    virtual bool remove_delete_partition(key_type k) {
        if (not partitions_->remove_delete(k)) {
            /* Local remove failed. */
            if (not partitions_->empty()) {
                auto hp = choose_partition(k); // find corresponding partition
                if (not hp) { // no valid partition exists
                    /* No key to delete. */
                    return false;
                } else return (*hp)->remove_delete_partition(k); // recursively delete k
            } else return false;
        }
        return true;
    }

    void clearEntries() {
        this->entries_->clear();
        std::vector<Key>* keys = this->partitions_->get_keys();
        for(key_type k: *(keys))
            this->partitions_->get(k).value()->clearEntries();
        delete(keys);
    }

    /**
     * Reset the NodeStatistics
     */
    void reset_statistics(bool recursive=true) {
        stats_->reset_statistics();
        if (recursive) {
            for(auto it=partitions_->key_begin(), end=partitions_->key_end(); it != end; ++it) {
                partitions_->get(*it).value()->reset_statistics(recursive);
            }
        }
    }

    void normalize_repetitions(uint64_t num_rep, bool recursive=true) {
        stats_->normalize_repetitions(num_rep);
        if (recursive) {
            for(auto it=partitions_->key_begin(), end=partitions_->key_end(); it != end; ++it) {
                partitions_->get(*it).value()->normalize_repetitions(num_rep, recursive);
            }
        }
    }


    /**
     * Map the NodeStatistics, i.e. insert it into the provided mapping from id_type to NodeStatistics*
     */
    void map_statistics(std::map<id_type, NodeStatistics<Key, Value>*>* mapping, bool recursive=true, std::set<id_type>* ids=nullptr, bool clone=false) {
        if (not ids or (ids and ids->find(myId) != ids->end())) {
            if (clone)
                mapping->insert(std::make_pair(myId, new NodeStatistics(stats_)));
            else    
                mapping->insert(std::make_pair(myId, stats_));
        }
        for (auto it=partitions_->value_begin(), end=partitions_->value_end(); it != end; ++it) {
            (*it)->map_statistics(mapping, recursive, ids, clone);
        }
    } 

    protected:
    void add_partition_search(std::string name, std::unique_ptr<IndexSearchMethod<key_type, hp_type*>> method){
        partition_search_methods_.insert({name, std::move(method)});
    };
    void add_entry_search(std::string name, std::unique_ptr<IndexSearchMethod<key_type, mapped_type>> method){
        entry_search_methods_.insert({name, std::move(method)});
    };


    /*==================================================================================================================
     * Capacity
     *================================================================================================================*/
    public:
    bool empty() {
        if (this->entries_->size() > 0) { return false; }
        std::vector<Key>* keys = this->partitions_->get_keys();
        for (key_type k: *(keys)) { if (not this->partitions_->get(k).value()->empty()) return false; }
        delete(keys);
        return true;
    }

    bool isEmptyPartition() const {
        if (this->entries_->size() > 0) return false;
        if (this->partitions_->size() > 0) return false;
        return true;
    }

    unsigned int countEmptyPartitions() {
        if (isEmptyPartition())
            return 1;
        unsigned int counter = 0;
        auto partition_it = partitions_->value_begin();
        auto partition_end = partitions_->value_end();
        while (partition_it != partition_end) {
            counter += (*partition_it)->countEmptyPartitions();
            ++partition_it;
        }
        return counter;
    }

    unsigned int countEntries() {
        unsigned int counter = entries_->size();
        auto partition_it = partitions_->value_begin();
        auto partition_end = partitions_->value_end();
        while(partition_it != partition_end) {
            counter += (*partition_it)->countEntries();
            ++partition_it;
        }
        return counter;
    }


    /*==================================================================================================================
     * Element access / Queries
     *================================================================================================================*/
    /** Returns the corresponding mapped_type value for the given key iff the key is present.  The lookup initially
     * searches the `entries_` node and if the search fails, recursively searches the nested horizontal partitions. */
    std::optional<mapped_type> pointQ(key_type k) {
        this->stats_->add_execution(0, Workload_kind::Point);
        if (not entries_->empty()) {
            if (auto v = entry_search_->lookup(entries_, k)) {
                return v;
            }
        }
        if (not partitions_->empty()) {
            if (auto hp = choose_partition(k)) {
                return (*hp)->pointQ(k);
            }
        }
        return std::nullopt;
    }

    void rangeQ(uint64_t &results, const key_type &lower, const key_type &upper) {
        this->stats_->add_execution(0, Workload_kind::Range);
        if (not entries_->empty()) {
            auto part_in_range = entry_search_->find_range(entries_, results, lower, upper); // returns an empty std::vector
            assert(part_in_range.empty());
        }
        if (not partitions_->empty()) {
            auto part_in_range = partition_search_->find_range(partitions_, results, lower, upper);
            for (auto part : part_in_range)
                part->rangeQ(results, lower, upper);
            return;
        }
    }

    /** Returns the corresponding partition for the given key.  The returned horizontal partition matches the insertion
     * semantics. */
    std::optional<hp_type*> choose_partition(key_type k) { return partition_search_->choose_subset(partitions_, k); }

    std::optional<hp_type*> lookupPartitionKey(key_type k) {
        auto v = partition_search_->lookup(partitions_, k);
        if (v.has_value()) return v;
        return std::nullopt;
    }

    std::optional<hp_type*> lookupPartitionKey(key_type k, id_type id) {
        if (myId == id) return this;
        auto v = partition_search_->lookup(partitions_, k);
        if (v.has_value() && (*v)->myId == id) return v;
        if (partitions_->empty()) return std::nullopt;
        if(auto hp = choose_partition(k)) {
            if (hp.has_value()) {
                if ((*hp)->myId == id) return hp.value();
                else return (*hp)->lookupPartitionKey(k, id);
            }
        }
        return std::nullopt;
    }

    std::optional<hp_type*> getPartitionContainingKey(key_type k, id_type id, bool containing=false) {
        auto v = partition_search_->lookup(partitions_, k);
        if (v.has_value() && (*v)->myId == id) {
            return this;
        }
        if (partitions_->empty()) {
            return std::nullopt;
        }
        if(auto hp = choose_partition(k)) {
            if (hp.has_value()) {
                if(containing and (*hp)->myId == id)
                    return this;
                return (*hp)->getPartitionContainingKey(k, id);  
            }
        }
        return std::nullopt;
    }


    /*==================================================================================================================
     * Getters / Setters
     *================================================================================================================*/
    /** Returns the number of entries in the horizontal partition. */
    size_type size() const { return entries_->size(); }
    size_type partition_size() const { return partitions_->size(); }
    size_type entry_size() const { return entries_->size(); }

    size_type partition_capacity() const { return partitions_->capacity(); }
    size_type entry_capacity() const { return entries_->capacity(); }
    void set_partition_capacity(size_type c) const { partitions_->set_capacity(c); }
    void set_entry_capacity(size_type c) const { entries_->set_capacity(c); }

    key_type max_key() const { return max_key_; }
    key_type min_key() const { return min_key_; }

    PartitionType getPartitionType() { return partitionType; }

    void set_partition_search_method(std::string name) { partition_search_ = partition_search_methods_.at(name).get(); }
    void set_entry_search_method(std::string name) { entry_search_ = entry_search_methods_.at(name).get(); }

    NodeStatistics<key_type, mapped_type> * get_stats() { return stats_; }

    void set_partition_search_method_rec(std::string name) {
        this->set_partition_search_method(name);
        auto part_it = partitions_->value_begin();
        auto part_end = partitions_->value_end();
        for (; part_it != part_end; ++part_it) (*part_it)->set_partition_search_method_rec(name);
    }
    void set_entry_search_method_rec(std::string name) {
        this->set_entry_search_method(name);
        auto part_it = partitions_->value_begin();
        auto part_end = partitions_->value_end();
        for (; part_it != part_end; ++part_it) (*part_it)->set_entry_search_method_rec(name);
    }

    Node<key_type, hp_type*> * get_partitions() const { return partitions_; }
    Node<key_type, mapped_type> * get_entries() const { return entries_; }

    IndexSearchMethod<key_type, hp_type*>* get_partition_search_method() { return partition_search_; }
    IndexSearchMethod<key_type, mapped_type>* get_entry_search_method() { return entry_search_; }

    std::vector<key_type>* getEntryKeys(bool recursive=true) {
        std::vector<key_type>* keys = new std::vector<key_type>();
        std::vector<key_type>* currentKeys = this->entries_->get_keys();
        for(Key k: *(currentKeys))
            keys->push_back(k);
        delete(currentKeys);
        if (!recursive)
            return keys;
        currentKeys = this->partitions_->get_keys();
        for(Key k: *(currentKeys)) {
            HorizontalPartition<key_type, mapped_type>* hp2 = this->partitions_->get(k).value();
            std::vector<key_type>* childKeys = hp2->getEntryKeys();
            for (Key k2 : (*childKeys))
                keys->push_back(k2);
            delete(childKeys);
        }
        delete(currentKeys);
        return keys;
    }

    std::vector<key_type>* getPartitionKeys(bool recursive=true) {
        std::vector<key_type>* keys = new std::vector<key_type>();
        std::vector<key_type>* currentKeys = partitions_->get_keys();
        for(Key k : *(currentKeys)) {
            keys->push_back(k);
            if (recursive) {
                std::optional<HorizontalPartition*> hp2 = partitions_->get(k);
                if (hp2) {
                    std::vector<Key>* childKeys = (*hp2)->getPartitionKeys();
                    for (Key k2 : *childKeys)
                        keys->push_back(k2);
                    delete(childKeys);
                }
            }
        }
        delete(currentKeys);
        return keys;
    }

    std::vector<std::string> * getPartitionSearchMethods() {
        std::vector<std::string> * methods = new std::vector<std::string>;
        for (auto it = partition_search_methods_.begin(); it != partition_search_methods_.end(); it++)
            methods->push_back(it->first);
        return methods;
    }

    std::vector<std::string>* getEntrySearchMethods() {
        std::vector<std::string>* methods = new std::vector<std::string>();
        for(auto it = entry_search_methods_.begin(); it != entry_search_methods_.end(); it++)
            methods->push_back(it->first);
        return methods;
    }

    id_type getId() { return myId; }
    void setId(id_type newId) { myId = newId; }

    void get_ids(std::set<id_type>* ids) {
        ids->insert(myId);
        for(auto it=partitions_->value_begin(), end=partitions_->value_end(); it != end; ++it)
            (*it)->get_ids(ids);
    }

    /*==================================================================================================================
     * Cloning
     *================================================================================================================*/
    /** clone this instance */
    virtual HorizontalPartition<key_type, mapped_type>* cloneGivenCopy(HorizontalPartition<key_type, mapped_type>* copy) {
        if (this->partitions_) {
            Node<key_type, hp_type*>* partitions = this->partitions_->cloneEmpty();
            auto key_it = this->partitions_->key_begin();
            auto value_it = this->partitions_->value_begin();
            auto key_end = this->partitions_->key_end();
            for( ; key_it != key_end; ++key_it, ++value_it)
                partitions->insert(*key_it, (*value_it)->clone());
            delete(copy->partitions_);
            copy->partitions_ = partitions;
        }
        if (this->entries_) {
            delete(copy->entries_);
            copy->entries_ = this->entries_->clone();
        }
        copy->partition_search_methods_.clear();
        for(auto const& [key, value] : this->partition_search_methods_) {
            std::string k = key;
            IndexSearchMethod<key_type, hp_type*>* val = value->clone();
            copy->partition_search_methods_.insert({k, std::move(std::unique_ptr<IndexSearchMethod<key_type, hp_type*>>(val))});
            if (value.get() == this->partition_search_) { copy->partition_search_ = val; }
        }
        copy->entry_search_methods_.clear();
        for(auto const& [key, value] : this->entry_search_methods_) {
            std::string k = key;
            IndexSearchMethod<key_type, mapped_type>* val = value->clone();
            copy->entry_search_methods_.insert({k, std::move(std::unique_ptr<IndexSearchMethod<key_type, mapped_type>>(val))});
            if (value.get() == this->entry_search_) { copy->entry_search_ = val; }
        }
        copy->partitionType = this->partitionType;
        copy->max_key_ = this->max_key_;
        if (this->stats_) {
            delete(copy->stats_);
            copy->stats_ = new NodeStatistics(this->stats_);
        }
        copy->myId = myId;
        return copy;
    }

    virtual HorizontalPartition<key_type, mapped_type>* clone() {
        HorizontalPartition<key_type, mapped_type>* copy = new HorizontalPartition<key_type, mapped_type>();
        return cloneGivenCopy(copy);
    }


    /*==================================================================================================================
     * Misc
     *================================================================================================================*/
    bool isValid(std::set<id_type> closed=std::set<id_type>()) {
        if (closed.find(myId) != closed.end()) {
            std::cerr << "Found invalid cycle" << std::endl;
            return false;
        }
        closed.insert(myId);
        std::vector<Key>* keysParent = partitions_->get_keys();
        if(keysParent->empty()) {
            delete(keysParent);
            return true;
        }
        std::sort(keysParent->begin(), keysParent->end());
        Key largestNonDefault = (*keysParent)[0];
        for(Key k: *(keysParent)) {
            if(k > largestNonDefault && k < std::numeric_limits<Key>::max())
                largestNonDefault = k;
            std::vector<Key>* keysChild = partitions_->get(k).value()->get_partitions()->get_keys();
            for(Key k2: *(keysChild)) {
                // found non-default children in default partition although parent is not full
                if(k > largestNonDefault && keysParent->size() < partitions_->capacity() && k2 < std::numeric_limits<Key>::max()) {
                    std::cerr << "Invalid default children under non-full parent combination: " << largestNonDefault << std::endl;
                    return false;
                }
                // found child key larger than parent key
                if(k2 > k && k2 < std::numeric_limits<Key>::max()) {
                    std::cerr << "Invalid parent-child-combination: " << k2 << " > " << k << std::endl;
                    return false;
                }
                // found child key in parent default partition smaller than largest non-default parent key
                if(k > largestNonDefault && k2 <= largestNonDefault) {
                    std::cerr << "Invalid default-child-combination: " << k2 << " < " << largestNonDefault << std::endl;
                    return false;
                }
            }
            delete(keysChild);
            if (not partitions_->get(k).value()->isValid(closed))
                return false;
        }
        delete(keysParent);
        return true;
    }

    void train_linreg_data() const {
        auto entry_linreg = static_cast<LinearRegressionSearch<key_type, mapped_type>*>(entry_search_methods_.at("LinearRegressionSearch").get());
        entry_linreg->train(entries_);
    }

    void train_linreg_partition() {
        auto partition_linreg = static_cast<LinearRegressionSearch<key_type, hp_type*>*>(partition_search_methods_.at("LinearRegressionSearch").get());
        partition_linreg->train(partitions_);
    }


    /*==================================================================================================================
     * Comparison
     *================================================================================================================*/
    bool operator==(const HorizontalPartition<key_type, mapped_type> &other) const {
        if (partitionType != other.partitionType) return false;

        auto other_entries = other.get_entries();
        auto other_partitions = other.get_partitions();

        /* Check entry node for equality. */
        auto e_key_it = entries_->key_begin();
        auto e_o_key_it = other_entries->key_begin();
        auto e_key_end = entries_->key_end();

        auto e_value_it = entries_->value_begin();
        auto e_o_value_it = other_entries->value_begin();

        if (entries_->size() != other_entries->size()) return false;
        if (entries_->capacity() != other_entries->capacity()) return false;
        for(; e_key_it != e_key_end; ++e_key_it, ++e_o_key_it, ++e_value_it, ++e_o_value_it)
            if ((*e_key_it != *e_o_key_it) or (*e_value_it != *e_o_value_it)) return false;

        /* Check partition node for equality. */
        auto p_key_it = partitions_->key_begin();
        auto p_o_key_it = other_partitions->key_begin();
        auto p_key_end = partitions_->key_end();

        auto p_value_it = partitions_->value_begin();
        auto p_o_value_it = other_partitions->value_begin();

        if (partitions_->size() != other_partitions->size()) return false;
        if (partitions_->capacity() != other_partitions->capacity()) return false;
        for(; p_key_it != p_key_end; ++p_key_it, ++p_o_key_it, ++p_value_it, ++p_o_value_it)
            if ((*p_key_it != *p_o_key_it) or (*(*p_value_it) != *(*p_o_value_it))) return false;

        if (partition_search_->getSearchType() != other.partition_search_->getSearchType()) return false;
        if (entry_search_->getSearchType() != other.entry_search_->getSearchType()) return false;

        return true;
    };
    bool operator!=(const HorizontalPartition<key_type, mapped_type> &other) const { return not (*this == other); }


    /*==================================================================================================================
     * Visualization
     *================================================================================================================*/
    void dot(std::ostream &out, bool include_partitions, bool include_entries, const std::string &zcode_prefix="Id") const {
        /* Configuration. */
        std::string color = "black";
        std::string type = "HorizontalPartition";

        /* Differentiate based on concrete HP type. */
        switch (this->partitionType) {
            case PartitionType::SoAPartition:
                    color = "orangered";
                    type  = "SortedArray (SoA)";
                    break;
            case PartitionType::TreePartition:
                    color = "deepskyblue";
                    type  = "Tree";
                    break;
            case PartitionType::HashPartition:
                    color = "palegreen";
                    type  = "Hash";
                    break;
            default:
                    break;
        }

        /* Differentiate based on concrete search method type. */
        std::string p_search = "SearchMethod";
        std::string p_color = "black";
        switch (partition_search_->getSearchType()) {
            case SearchType::DefaultSearchMethod:
                p_search = "Default Search";
                p_color = "lightpink";
                break;
            case SearchType::LinearSearchMethod:
                p_search = "Linear Search";
                p_color = "aquamarine";
                break;
            case SearchType::BinarySearchMethod:
                p_search = "Binary Search";
                p_color = "burlywood";
                break;
            case SearchType::InterpolationSearchMethod:
                p_search = "Interpolation Search";
                p_color = "green";
                break;
            case SearchType::ExponentialSearchMethod:
                p_search = "Exponential Search";
                p_color = "yellow";
                break;
            case SearchType::LinearRegressionSearchMethod:
                p_search = "Linear Regression Search";
                p_color = "indigo";
                break;
            default:
                break;
        }

        std::string e_search = "SearchMethod";
        std::string e_color = "black";
        switch (entry_search_->getSearchType()) {
            case SearchType::DefaultSearchMethod:
                e_search = "Default Search";
                e_color = "lightpink";
                break;
            case SearchType::LinearSearchMethod:
                e_search = "Linear Search";
                e_color = "aquamarine";
                break;
            case SearchType::BinarySearchMethod:
                e_search = "Binary Search";
                e_color = "burlywood";
                break;
            case SearchType::InterpolationSearchMethod:
                e_search = "Interpolation Search";
                e_color = "green";
                break;
            case SearchType::ExponentialSearchMethod:
                e_search = "Exponential Search";
                e_color = "yellow";
                break;
            case SearchType::LinearRegressionSearchMethod:
                e_search = "Linear Regression Search";
                e_color = "indigo";
                break;
            default:
                break;
        }


        auto p_capacity = partition_capacity();
        auto e_capacity = entry_capacity();

        auto include_p = include_partitions && p_capacity > 0;
        auto include_e = include_entries && e_capacity > 0;

        auto node_width = 0;
        node_width = include_p ? p_capacity : node_width;
        node_width = include_e ? node_width + e_capacity : node_width;
        node_width = (include_p and include_e) ? node_width + 1 : node_width;

        if (include_p or include_e) {
            out << "node" << zcode_prefix << " [label=<\n<table cellspacing=\"0\" border=\"0\" cellborder=\"1\" color=\""
                << color << "\">\n"
                << "\t<tr><td colspan=\"" << std::to_string(node_width) << "\">" << type 
                << " ID: " << myId
                << " Runtime: " << stats_->get_total_time()
                << "</td></tr>\n";

            /* Add search methods. */
            out << "\t<tr>";
            if (include_p)
                out << "<td colspan=\"" << std::to_string(p_capacity) << "\">" << p_search << "</td>";
            if (include_p and include_e)
                out << "<td BGCOLOR=\"" << color << "\"></td>";
            if (include_e)
                out << "<td colspan=\"" << std::to_string(e_capacity) << "\">" << e_search << "</td>";
            out << "</tr>\n";

            /* Keys. */
            {
                out << "\t<tr>";
                if (include_p) {
                    /* Partition keys. */
                    auto pkeys_it = partitions_->key_begin();
                    auto p_size = partitions_->size();
                    for (size_type i = 0; i < p_capacity; ++i) {
                        if (i < p_size) {
                            out << "<td>" << *pkeys_it << "</td>";
                            ++pkeys_it;
                        } else out << "<td></td>";
                    }
                }

                if (include_p and include_e)
                    out << "<td BGCOLOR=\"" << color << "\"></td>";

                if (include_e) {
                    /* Entry keys. */
                    auto ekeys_it = entries_->key_begin();
                    auto e_size = entries_->size();
                    for (size_type i = 0; i < e_capacity; ++i) {
                        if (i < e_size) {
                            out << "<td>" << *ekeys_it << "</td>";
                            ++ekeys_it;
                        } else out << "<td></td>";
                    }
                }
                out << "</tr>\n";
            }

            /* Values. */
            {
                out << "\t<tr>";
                if (include_p) {
                    /* Partition values. */
                    auto p_size = partitions_->size();
                    for (size_type i = 0; i < p_capacity; ++i) {
                        if (i < p_size)
                            out << "<td PORT=\"P" << i << "\">P" << i << "</td>";
                        else out << "<td></td>";
                    }
                }

                if (include_p and include_e)
                    out << "<td BGCOLOR=\"" << color << "\"></td>";

                if (include_e) {
                    /* Entry values. */
                    auto evalues_it = entries_->value_begin();
                    auto e_size = entries_->size();
                    for (size_type i = 0; i < e_capacity; ++i) {
                        if (i < e_size) {
                            out << "<td>" << *evalues_it << "</td>";
                            ++evalues_it;
                        } else out << "<td></td>";
                    }
                }
                out << "</tr>\n";
            }

            out << "</table>>];\n";
        } else {
            /* Minimal Visualization. */
            std::string label = "ID: " + std::to_string(myId) + \
                                " | Runtime: " + std::to_string(stats_->get_total_time()) + \
                                "<BR/>PQ: " + std::to_string(stats_->get_num_point()) + \
                                " | RQ: " + std::to_string(stats_->get_num_range()) + \
                                " | IN: " + std::to_string(stats_->get_num_insert()) + \
                                " | #p: " + std::to_string(partition_size()) + \
                                " | #e: " + std::to_string(entry_size());
            out << "node" << zcode_prefix << " [label=<\n<table cellspacing=\"0\" border=\"0\" cellborder=\"1\" color=\""
                << color << "\">\n"
                << "\t<tr><td bgcolor=\"" << p_color << "\"> </td><td bgcolor=\"" << color << "\">" << label << "</td><td bgcolor=\"" << e_color << "\"> </td></tr></table>>];\n";
        }


        uint32_t zcode_length = std::ceil(std::log10(p_capacity));
        auto p_size = partitions_->size();
        auto pvalues_it = partitions_->value_begin();
        for (size_type i = 0; i < p_size; ++i) {
            std::stringstream stream;
            stream << std::setw(zcode_length) << std::setfill('0') << i;
            std::string prefix(zcode_prefix + stream.str());
            (*pvalues_it)->dot(out, include_partitions, include_entries, prefix);
            if (include_p)
                out << "node" << zcode_prefix << ":P" << std::to_string(i) << " -> node" << prefix << ";\n";
            else
                out << "node" << zcode_prefix << " -> node" << prefix << ";\n";
            ++pvalues_it;
        }
    }

    public:
    void to_graphviz(std::ostream &out, bool include_partitions, bool include_entries) const {
        out << "digraph index {\n"
            << "node [shape=plain, height=.1];\n";
        dot(out, include_partitions, include_entries);
        out << "}\n";
    }
    void to_graphviz(std::ostream &out) const { to_graphviz(out, true, true); }

    friend std::string to_string(const hp_type &hp) {
        std::ostringstream os;
        os << hp;
        return os.str();
    }

    /** Write the `entries` of this horizontal partition to `out`. */
    friend std::ostream & operator<<(std::ostream &out, HorizontalPartition<key_type, mapped_type> &hp) {
        out << "partitions:\n";
        if (auto p = hp.get_partitions()) p->print(out);
        out << "entries:\n";
        if (auto e = hp.get_entries()) e->print(out);
        return out;
    }

    friend std::ostream & operator<<(std::ostream &out, const HorizontalPartition<key_type, mapped_type> &hp) {
        out << "partitions:\n";
        if (auto p = hp.get_partitions()) p->print(out);
        out << "entries:\n";
        if (auto e = hp.get_entries()) e->print(out);
        return out;
    }
};


/*======================================================================================================================
 * SortedArray
 *====================================================================================================================*/
template<typename Key, typename Value>
struct SortedArray : public HorizontalPartition<Key, Value>
{
    using Base = HorizontalPartition<Key, Value>;
    IMPORT(size_type);
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(hp_type);

    SortedArray(std::string container_type, size_type partition_capacity, size_type entry_capacity) {
        if (container_type == "SoA") {
            Base::partitions_ = new SoA<key_type, hp_type*>(partition_capacity);
            Base::entries_  = new SoA<key_type, mapped_type>(entry_capacity);
            Base::partitionType = PartitionType::SoAPartition;
            Base::stats_->set_partition_type(PartitionType::SoAPartition);
        } else throw std::invalid_argument("Not a valid data structure for a sorted array");

        /* Add valid search methods. */
        Base::add_partition_search("LinearSearch", std::move(std::make_unique<LinearSearch<key_type, hp_type*>>()));
        Base::add_entry_search("LinearSearch", std::move(std::make_unique<LinearSearch<key_type, mapped_type>>()));

        Base::add_partition_search("BinarySearch", std::move(std::make_unique<BinarySearch<key_type, hp_type*>>()));
        Base::add_entry_search("BinarySearch", std::move(std::make_unique<BinarySearch<key_type, mapped_type>>()));

        Base::add_partition_search("InterpolationSearch", std::move(std::make_unique<InterpolationSearch<key_type, hp_type*>>()));
        Base::add_entry_search("InterpolationSearch", std::move(std::make_unique<InterpolationSearch<key_type, mapped_type>>()));

        Base::add_partition_search("ExponentialSearch", std::move(std::make_unique<ExponentialSearch<key_type, hp_type*>>()));
        Base::add_entry_search("ExponentialSearch", std::move(std::make_unique<ExponentialSearch<key_type, mapped_type>>()));

        Base::add_partition_search("LinearRegressionSearch", std::move(std::make_unique<LinearRegressionSearch<key_type, hp_type*>>()));
        Base::add_entry_search("LinearRegressionSearch", std::move(std::make_unique<LinearRegressionSearch<key_type, mapped_type>>()));

        /* Set search methods. */
        Base::set_partition_search_method("BinarySearch");
        Base::set_entry_search_method("BinarySearch");
    }
    SortedArray(std::string container_type, size_type capacity) : SortedArray<key_type, mapped_type>(container_type, capacity, capacity) { }
    SortedArray(std::string container_type) : SortedArray<key_type, mapped_type>(container_type, 4, 4) { }

    std::pair<bool, std::optional<HorizontalPartition<key_type, mapped_type>*>> insert(key_type k, mapped_type v, bool split=true) override { return Base::insert(k, v, split); }
    
    /**
     * Insert a given key and value into the partition, descending recursively if necessary
     * If a leaf does not have enough space for the new entry, split the leaf to create more space
     * Returns a boolean indicator of success as well as an optional HP resulting from a split which has to be inserted in the parent
     */
    std::pair<bool, std::optional<HorizontalPartition<key_type, mapped_type>*>> insert_(key_type k, mapped_type v, bool split=true) override {
        Base::stats_->add_execution(0, Workload_kind::Insert);
        if (auto soa = dynamic_cast<SoA<key_type, mapped_type>*>(Base::get_entries())) {
            if (soa->insert_sorted(k, v)) { // may raise `DuplicateKeyException`
                /* Local insert successful. */
                try {
                    auto lr_search = dynamic_cast<LinearRegressionSearch<key_type, mapped_type>*>(Base::entry_search_methods_.at("LinearRegressionSearch").get());
                    lr_search->untrain();
                } catch (std::out_of_range &e) { /* do nothing */ }

                if (k > Base::max_key_) Base::max_key_ = k;
                if (k < Base::min_key_) Base::min_key_ = k;
                return std::make_pair(true, std::nullopt);
            } else {
                /* Local insert failed, insert recursively. */
                std::pair<bool, std::optional<HorizontalPartition<key_type, mapped_type>*>> result;
                if (not Base::partitions_->empty()) {
                    /* Inner node */
                    auto hp = Base::choose_partition(k); // find corresponding partition
                    if (not hp) { // no valid partition exists
                        auto partition_keys = Base::getPartitionKeys(false);
                        std::sort(partition_keys->begin(), partition_keys->end());
                        key_type largest_key = (*partition_keys)[partition_keys->size()-1];
                        assert(largest_key < k);
                        auto largest_partition = Base::lookupPartitionKey(largest_key).value();
                        Base::remove_partition(largest_key);
                        insert(k, largest_partition, false);
                        delete(partition_keys);
                        result = largest_partition->insert_(k, v, split);
                    } else {
                        result = (*hp)->insert_(k, v, split); // continue to traverse horizontal partitions
                    }
                    if (result.second) {
                        // child hp was split, insert newly created hp in this hp
                        if (Base::partitions_->size() >= Base::partitions_->capacity()) {
                            auto new_partition = Base::split_partition();
                            if (result.second.value()->max_key() > new_partition->max_key())
                                insert(result.second.value()->max_key(), result.second.value(), false);
                            else
                                new_partition->insert(result.second.value()->max_key(), result.second.value(), false);
                            return std::make_pair(result.first, new_partition);
                        } else {
                            insert(result.second.value()->max_key(), result.second.value(), false);
                            return std::make_pair(result.first, std::nullopt);
                        }
                    } else {
                        return result;
                    }            
                } else {
                    /* Leaf node */
                    if (split) {
                        auto new_partition = Base::split_partition();
                        if (k > new_partition->max_key()) {
                            insert_(k, v, false);
                        } else {
                            new_partition->insert_(k, v, false);
                        }
                        return std::make_pair(true, new_partition);
                    } else {
                        Base::stats_->add_penalty(Configuration::Get().penaltyFactorMissingKeys); // penalty
                        return std::make_pair(false, std::nullopt);    
                    }
                }
            }
        } 
        else throw std::invalid_argument("no valid node for sorted insert");
    }

    bool insert(key_type k, hp_type *hp, bool default_partition=true) override {
        if (auto soa = dynamic_cast<SoA<key_type, hp_type*>*>(Base::get_partitions())) {
            if (soa->insert_sorted(k, hp)) { // may raise `DuplicateKeyException`
                /* Local insert successful. */
                try {
                    auto lr_search = dynamic_cast<LinearRegressionSearch<key_type, hp_type*>*>(Base::partition_search_methods_.at("LinearRegressionSearch").get());
                    lr_search->untrain();
                } catch (std::out_of_range &e) { /* do nothing */ }

                hp->get_stats()->set_key(k);
                if (k > Base::max_key_) Base::max_key_ = k;
                if (k < Base::min_key_) Base::min_key_ = k;
                if(soa->size() == 1 and default_partition)
                    if (not soa->insert_sorted(std::numeric_limits<key_type>::max(), new SortedArray<key_type, mapped_type>("SoA", Base::partition_capacity(), Base::entry_capacity()))) {
                        return false;
                    }
                return true;
            } else {
                /* Local insert failed, insert recursively. */
                if (not Base::partitions_->empty()) {
                    auto rec_hp = Base::choose_partition(k); // find corresponding partition
                    if (rec_hp) // valid partition exists
                        return (*rec_hp)->insert(k, hp, default_partition); // continue to traverse horizontal partitions
                }
                return false;
            }
        }
        else throw std::invalid_argument("no valid node for sorted insert");
    }
    
    /** Write the `entries` of this horizontal partition to `out`. */
    friend std::ostream & operator<<(std::ostream &out, SortedArray<key_type, mapped_type> &hp) {
        out << "horizontal sorted array partition:\n";
        out << "partitions:\n";
        if (auto sorted_hp = hp.get_partitions()) sorted_hp->print(out);
        out << "entries:\n";
        if (auto sorted_entries = hp.get_entries()) sorted_entries->print(out);
        return out;
    }

    HorizontalPartition<key_type, mapped_type>* clone() override {
        if(this->partitionType == PartitionType::SoAPartition) {
            HorizontalPartition<key_type, mapped_type>* hpCopy = new SortedArray<key_type, mapped_type>("SoA", Base::partition_capacity(), Base::entry_capacity());
            return this->cloneGivenCopy(hpCopy);
        } else
            throw 666;
    }
};


/*======================================================================================================================
 * Hash
 *====================================================================================================================*/
template<typename Key, typename Value>
struct Hash : public HorizontalPartition<Key, Value>
{
    using Base = HorizontalPartition<Key, Value>;
    IMPORT(size_type);
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(hp_type);

    /* NOTE: Currently, finding a nested partition in a hash node is only supported by an equality check, i.e. given a
     * key-hp mapping in a hash node, the corresponding horizontal partition is only found if the search key matches the
     * key in the mapping. Therefore, currently it makes only sense to use hash partitions as leaf nodes.
     */
    public:
    Hash(size_type partition_capacity, size_type entry_capacity) {
        // assert(partition_capacity == 0);
        Base::partitions_ = new HashNode<Key, hp_type*>(partition_capacity);
        Base::entries_ = new HashNode<Key, mapped_type>(entry_capacity);

        /* Add valid search methods. */
        Base::add_partition_search("DefaultSearch", std::move(std::make_unique<DefaultSearch<key_type, hp_type*>>()));
        Base::add_entry_search("DefaultSearch", std::move(std::make_unique<DefaultSearch<key_type, mapped_type>>()));

        Base::add_partition_search("LinearSearch", std::move(std::make_unique<LinearSearch<key_type, hp_type*>>()));
        Base::add_entry_search("LinearSearch", std::move(std::make_unique<LinearSearch<key_type, mapped_type>>()));

        /* Set search methods. */
        Base::set_partition_search_method("LinearSearch");
        Base::set_entry_search_method("DefaultSearch");

        Base::partitionType = PartitionType::HashPartition;
        Base::stats_->set_partition_type(PartitionType::HashPartition);
    }

    Hash(size_type capacity) : Hash<key_type, mapped_type>(0, capacity) { }
    Hash() : Hash<key_type, mapped_type>(0, 4) { }

    std::pair<bool, std::optional<HorizontalPartition<key_type, mapped_type>*>> insert(key_type k, mapped_type v, bool split=true) override {
        return Base::insert(k, v, split);
    }

    bool insert(key_type, hp_type *, bool) override { throw InvalidInsertionException(); }

    /** Write the `entries` of this horizontal partition to `out`. */
    friend std::ostream & operator<<(std::ostream &out, Hash<Key, Value> &hp) {
        out << "horizontal hash partition:\n";
        out << "partitions:\n";
        if (auto partitions = hp.get_partitions()) partitions->print(out);
        out << "entries:\n";
        if (auto entries = hp.get_entries()) entries->print(out);

    return out;
}

    HorizontalPartition<key_type, mapped_type>* clone() override {
        HorizontalPartition<key_type, mapped_type>* hpCopy = new Hash<key_type, mapped_type>(Base::partition_capacity(), Base::entry_capacity());
        return this->cloneGivenCopy(hpCopy);
    }
};


/*======================================================================================================================
 * Tree
 *====================================================================================================================*/
template<typename Key, typename Value>
struct Tree : public HorizontalPartition<Key, Value>
{
    using Base = HorizontalPartition<Key, Value>;
    IMPORT(size_type);
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(hp_type);

    Tree(size_type partition_capacity, size_type entry_capacity) {
        Base::partitions_ = new TreeNode<Key, hp_type*>(partition_capacity);
        Base::entries_ = new TreeNode<Key, mapped_type>(entry_capacity);
        Base::partitionType = PartitionType::TreePartition;
        Base::stats_->set_partition_type(PartitionType::TreePartition);

        /* Add valid search methods. */
        Base::add_partition_search("LinearSearch", std::move(std::make_unique<LinearSearch<key_type, hp_type*>>()));
        Base::add_entry_search("LinearSearch", std::move(std::make_unique<LinearSearch<key_type, mapped_type>>()));

        Base::add_partition_search("BinarySearch", std::move(std::make_unique<BinarySearch<key_type, hp_type*>>()));
        Base::add_entry_search("BinarySearch", std::move(std::make_unique<BinarySearch<key_type, mapped_type>>()));

        /* Set search methods. */
        Base::set_partition_search_method("BinarySearch");
        Base::set_entry_search_method("BinarySearch");
    }

    Tree(size_type capacity) : Tree<key_type, mapped_type>(capacity, capacity) { }
    Tree() : Tree<key_type, mapped_type>(4, 4) { }

    /** Write the `entries` of this horizontal partition to `out`. */
    friend std::ostream & operator<<(std::ostream &out, Tree<Key, Value> &hp) {
        out << "horizontal tree partition:\n";
        out << "partitions:\n";
        if (auto partitions = hp.get_partitions()) partitions->print(out);
        out << "entries:\n";
        if (auto entries = hp.get_entries()) entries->print(out);

        return out;
    }

    HorizontalPartition<key_type, mapped_type>* clone() override {
        HorizontalPartition<key_type, mapped_type>* hpCopy = new Tree<key_type, mapped_type>(Base::partition_capacity(), Base::entry_capacity());
        return this->cloneGivenCopy(hpCopy);
    }
};

/* ======================================== IMPLEMENTATION ===========================================================*/
template<typename Key, typename Value>
bool HorizontalPartition<Key, Value>::insert(key_type k, hp_type *hp, bool default_partition) {
    if (partitions_->insert(k, hp)) { // may raise `DuplicateKeyException`
        /* Local insert successful. */
        try {
            auto lr_search = dynamic_cast<LinearRegressionSearch<key_type, hp_type*>*>(partition_search_methods_.at("LinearRegressionSearch").get());
            lr_search->untrain();
        } catch (std::out_of_range &e) { /* do nothing */ }

        hp->get_stats()->set_key(k);
        if (k > max_key_) max_key_ = k;
        if (k < min_key_) min_key_ = k;
        if (partitions_->size() == 1 and default_partition) {
            /* Insert default max. partition. */
            if (not partitions_->insert(std::numeric_limits<key_type>::max(), new SortedArray<key_type, mapped_type>("SoA", partition_capacity(), entry_capacity()))) {
                return false;
            }
        }
        return true;
    } else {
        /* Local insert failed, insert recursively. */
        if (not partitions_->empty()) {
            auto rec_hp = choose_partition(k); // find corresponding partition
            if (rec_hp) // valid partition exists
                return (*rec_hp)->insert(k, hp, default_partition); // continue to traverse horizontal partitions
        }
        return false;
    }
}

template<typename Key, typename Value>
unsigned long HorizontalPartition<Key, Value>::ID = 0;
#undef IMPORT
