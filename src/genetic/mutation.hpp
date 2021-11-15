#pragma once

#include "individual.hpp"
#include "util/globals.hpp"
#include "util/enum.hpp"
#include <algorithm>
#include <limits>

using id_type = unsigned long;

/*===================================================================================================
 * Mutation Base Class
 * ================================================================================================*/

/**
 * Abstract base class for all mutations
 */
template<typename Key, typename Value>
class Mutation {
    public:
    Mutation(double prob, const MutationType type)
        : probability(prob),
          type(type),
          executionCounter(0)
    {}

    virtual ~Mutation() {}

    /** return the probability associated with this mutation */
    double getProbability() { return probability; }

    /** perform this mutation based on the given individual */
    virtual const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) = 0;

    bool operator <(const Mutation<Key, Value>& other) const { (void) other; return false; }

    /** return the execution counter denoting how often this mutation has been performed so far */
    unsigned int getExecutionCounter() { return executionCounter; }

    /** return the name of this mutation */
    const std::string getName() { return toString(type); }

    /** return the type of this mutation */
    MutationType getType() { return type; }

    private:
    double probability;
    const MutationType type;

    protected:
    unsigned int executionCounter;
    Key last_key_mutated;
};

/*===================================================================================================
 * Mutation Helper Functions
 * ================================================================================================*/

/**
 * draw a random id according to the node statistics mapping
 * the method constructs a random distribution based on the total time spent in each node
 * it then draws a random id based on this distribution
 * if `idsToChooseFrom` is defined, only these ids are considered
 */
template<typename Key, typename Value>
id_type getRandomId(std::map<id_type, NodeStatistics<Key, Value>*> *mapping, 
        const std::set<id_type> *idsToChooseFrom=nullptr, bool use_time=false, bool use_calls=false) {
    std::mt19937 gen(std::random_device{}());
    std::vector<id_type> ids;
    std::vector<long double> weights;
    long double min_weight = std::numeric_limits<long double>::max();
    if(idsToChooseFrom && idsToChooseFrom->size() > 0) {
        if (not mapping) {
            throw std::runtime_error("Mapping is null");
        }
        for (auto it = idsToChooseFrom->begin(); it != idsToChooseFrom->end(); it++) {
            if (mapping->find(*it) != mapping->end()) {
                ids.push_back(*it);
                if (use_time)
                    weights.push_back((*mapping)[*it]->get_total_time());   
                else if (use_calls)
                    weights.push_back((*mapping)[*it]->get_num_calls());   
                else
                    weights.push_back(1);
                if (weights[weights.size() - 1] > 0 and weights[weights.size() - 1] < min_weight)
                    min_weight = weights[weights.size() - 1];
            } else {
                std::cout << "Could not find ID: " << *it << " in mapping" << std::endl;
                throw std::runtime_error("Key missing in mapping");
            }
        }
    } else {
        if (not mapping) {
            throw std::runtime_error("Mapping is null");
        }
        for (auto it = mapping->begin(); it != mapping->end(); it++) {
            ids.push_back(it->first);
            if (use_time)
                weights.push_back(it->second->get_total_time());   
            else if (use_calls)
                weights.push_back(it->second->get_num_calls());   
            else
                weights.push_back(1);   
            if (weights[weights.size() - 1] > 0 and weights[weights.size() - 1] < min_weight)
                min_weight = weights[weights.size() - 1];
        }
    }
    // use minimum, non-zero weight determined above
    // if minimum is still the initial value (i.e. no non-zero minimum was found), 
    // set all weights to 1 to get a valid distribution
    // else: set all zero weights to 10% of the minimum weight to consider them at least a few times
    for(size_t i=0; i<weights.size(); i++) {
        if(weights[i] == 0) {
            if(min_weight < std::numeric_limits<long double>::max()) {
                weights[i] = 0.1 * min_weight;
            } else {
                weights[i] = 1;
            }
        }
    }
    std::discrete_distribution<id_type> dist (weights.begin(), weights.end());
    auto result = ids[dist(gen)];
    return result;
} 

/** create a new partition with random data layout using the same capacities as the given HorizontalPartition */
template<typename Key, typename Value>
HorizontalPartition<Key, Value>* createPartition(HorizontalPartition<Key, Value>* hpCopy=nullptr, PartitionType preferred=PartitionType::UnknownPartition) {
    std::mt19937 gen(std::random_device{}());
    if (preferred == PartitionType::UnknownPartition) {
        std::vector<PartitionType> options{PartitionType::SoAPartition, PartitionType::TreePartition};
        std::discrete_distribution<int> dist {1, 1};
        preferred = options[dist(gen)];
    }
    switch(preferred) {
        case PartitionType::SoAPartition: { // SoA
            auto partition = new SortedArray<Key, Value>("SoA", hpCopy->partition_capacity(), hpCopy->entry_capacity());
            return partition; }
        case PartitionType::TreePartition: { // Tree
            auto partition = new Tree<Key, Value>(hpCopy->partition_capacity(), hpCopy->entry_capacity());
            return partition; }
        default: // Hash or other HPs not supported
            throw std::invalid_argument("Unsupported partition type.");
    }
}

/** 
 * Replace the given HP by creating a new HP with random type and search methods and copying the entries over 
 * If `newType` is not `PartitionType::UnknownPartition`, it strictly defined the new type
 * Otherwise the new type is randomly chosen based on the given `priorities` by building a corresponding distribution
 */
template<typename Key, typename Value>
HorizontalPartition<Key, Value>* replacePartition(HorizontalPartition<Key, Value>* hp, PartitionType newType, 
                                                std::map<PartitionType, double> priorities=std::map<PartitionType, double>(), 
                                                std::map<std::string, double> searchPriorities=std::map<std::string, double>(),
                                                bool allowSearchMethodChange=false) {
    HorizontalPartition<Key, Value>* hpCopy = nullptr;
    std::mt19937 gen(std::random_device{}());
    // create a new partition based on the given type
    // if newType is UnknownPartition, create a partition with random data layout
    switch(newType) {
        case PartitionType::SoAPartition:
            hpCopy = new SortedArray<Key, Value>("SoA", hp->partition_capacity(), hp->entry_capacity());
            break;
        case PartitionType::TreePartition:
            hpCopy = new Tree<Key, Value>(hp->partition_capacity(), hp->entry_capacity());
            break;
        case PartitionType::HashPartition:
            hpCopy = new Hash<Key, Value>(hp->partition_capacity(), hp->entry_capacity());
            break;
        // create a partition with random layout
        case PartitionType::UnknownPartition:
            // create probability distribution based on the possible types
            // either take default priority of 1 or the priority given in priorities
            std::vector<double> probabilities;
            std::vector<PartitionType> possibleTypes = {PartitionType::SoAPartition, PartitionType::TreePartition, PartitionType::HashPartition};
            for (auto t : possibleTypes) {
                if(hp->getPartitionType() == t || (t == PartitionType::HashPartition && hp->get_partitions()->size() > 0))
                    probabilities.push_back(0.0);
                else {
                    if (priorities.find(t) != priorities.end())
                        probabilities.push_back(priorities[t]);
                    else
                        probabilities.push_back(1.0);
                }
            }
            std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
            // draw data layout and create partition
            auto tNew = possibleTypes[distribution(gen)];
            switch(tNew) {
                case PartitionType::SoAPartition:
                    hpCopy = new SortedArray<Key, Value>("SoA", hp->partition_capacity(), hp->entry_capacity());
                    break;
                case PartitionType::TreePartition:
                    hpCopy = new Tree<Key, Value>(hp->partition_capacity(), hp->entry_capacity());
                    break;
                case PartitionType::HashPartition:
                    hpCopy = new Hash<Key, Value>(hp->partition_capacity(), hp->entry_capacity());
                    break;
                default:
                    throw std::invalid_argument("Unexpected random choice");
            }
            break;
    }

    // copy over partitions
    auto p_key_it = hp->get_partitions()->key_begin();
    auto p_key_end = hp->get_partitions()->key_end();
    auto p_value_it = hp->get_partitions()->value_begin();
    for( ; p_key_it != p_key_end; ++p_key_it, ++p_value_it) {
        hpCopy->remove_delete_partition(*p_key_it);
        try {
            //hpCopy->insert(*p_key_it, (*p_value_it)->clone(), false);
            hpCopy->insert(*p_key_it, *p_value_it, false);
        }
        catch(std::bad_optional_access& e) {
            std::cerr << "ReplacePartition: Bad partitions lookup: " << *p_key_it << std::endl;
            throw e;
        }
        catch (DuplicateKeyException& e) {
            std::cerr << "DuplicateKeyException in " << __FUNCTION__ << " case 1 with key " << *p_key_it << std::endl;
            throw e;
        }
    }
    auto p_keys = hp->get_partitions()->get_keys();
    for(auto k = p_keys->begin(); k != p_keys->end(); ++k) {
        hp->remove_partition(*k);
    }
    delete(p_keys);
        
    // copy over entries
    auto e_key_it = hp->get_entries()->key_begin();
    auto e_key_end = hp->get_entries()->key_end();
    auto e_value_it = hp->get_entries()->value_begin();
    for( ; e_key_it != e_key_end; ++e_key_it, ++e_value_it) {
        try {
            hpCopy->insert(*e_key_it, *e_value_it);
        } 
        catch(std::bad_optional_access& e) {
            std::cerr << "ReplacePartition: Bad entries lookup: " << *e_key_it << std::endl;
            throw e;
        }
        catch (DuplicateKeyException& e) {
            std::cerr << "DuplicateKeyException in " << __FUNCTION__ << " case 2" << std::endl;
            throw e;
        }
    }
    auto e_keys = hp->get_entries()->get_keys();
    for(auto k = e_keys->begin(); k != e_keys->end(); ++k) {
        hp->remove_entry(*k);
    }
    delete(e_keys);
    
    // copy over id
    hpCopy->setId(hp->getId());

    // determine search method for new partition
    std::uniform_int_distribution boolDist (0, 1);
    std::vector<std::string>* methods = hpCopy->getPartitionSearchMethods();
    if (allowSearchMethodChange) {
        /* Alternative 1:
         * Choose search methods according to following method:
         * with 50% chance: choose random search method
         * otherwise: try to set same search method as in old node type
         * if not available: keep default search method set in constructor */
        if (boolDist(gen)) {
            std::vector<double> probabilities;
            for (auto method : *methods) {
                if (searchPriorities.find(method) != searchPriorities.end())
                    probabilities.push_back(searchPriorities[method]);
                else
                    probabilities.push_back(1.0);
            }
            std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
            std::mt19937 gen(std::random_device{}());
            std::string newMethod = (*methods)[distribution(gen)];
            if (newMethod != toString(hpCopy->get_partition_search_method()->getSearchType()))
                hpCopy->set_partition_search_method(newMethod);
        } else {
            std::string newMethod = toString(hp->get_partition_search_method()->getSearchType());
            if ((not newMethod.empty()) && (std::find(methods->begin(), methods->end(), newMethod) != methods->end()) && (newMethod != toString(hpCopy->get_partition_search_method()->getSearchType())))
                hpCopy->set_partition_search_method(newMethod);
        }
    } else {
        /* Alternative 2:
         * Always try to keep previous search method */
        std::string newMethod = toString(hp->get_partition_search_method()->getSearchType());
        if ((not newMethod.empty()) && (std::find(methods->begin(), methods->end(), newMethod) != methods->end()) && (newMethod != toString(hpCopy->get_partition_search_method()->getSearchType())))
            hpCopy->set_partition_search_method(newMethod);
    }
    delete(methods);
    
    methods = hpCopy->getEntrySearchMethods();
    if (allowSearchMethodChange) {
        /* Alternative 1 */
        if (boolDist(gen)) {
            std::vector<double> probabilities;
            for (auto method : *methods) {
                if (searchPriorities.find(method) != searchPriorities.end())
                    probabilities.push_back(searchPriorities[method]);
                else
                    probabilities.push_back(1.0);
            }
            std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
            std::mt19937 gen(std::random_device{}());
            std::string newMethod = (*methods)[distribution(gen)];
            if (newMethod != toString(hpCopy->get_entry_search_method()->getSearchType()))
                hpCopy->set_entry_search_method(newMethod);
        } else {
            std::string newMethod = toString(hp->get_entry_search_method()->getSearchType());
            if ((not newMethod.empty()) && (std::find(methods->begin(), methods->end(), newMethod) != methods->end()) && (newMethod != toString(hpCopy->get_entry_search_method()->getSearchType())))
                hpCopy->set_entry_search_method(newMethod);
        }
    } else {
        /* Alternative 2 */
        std::string newMethod = toString(hp->get_entry_search_method()->getSearchType());
        if ((not newMethod.empty()) && (std::find(methods->begin(), methods->end(), newMethod) != methods->end()) && (newMethod != toString(hpCopy->get_entry_search_method()->getSearchType())))
            hpCopy->set_entry_search_method(newMethod);
    }
    delete(methods);

    return hpCopy;
}

/*===================================================================================================
 * Mutation Functions On Horizontal Partitions
 * ================================================================================================*/

/** Change type of HP */
template<typename Key, typename Value>
std::optional<HorizontalPartition<Key, Value>*> changePartitionType(HorizontalPartition<Key, Value>* hp, Key key, id_type id, bool clone=true, 
                                                                    PartitionType newType=PartitionType::UnknownPartition, 
                                                                    std::map<PartitionType, double> priorities=std::map<PartitionType, double>(),
                                                                    std::map<std::string, double> searchPriorities=std::map<std::string, double>()) {
    HorizontalPartition<Key, Value>* hpCopy = hp;
    if (clone)
        hpCopy = hp->clone();
    // if key is not present, change type of root partition
    auto lookupResult = hpCopy->lookupPartitionKey(key, id);
    if(not lookupResult || lookupResult.value() == hpCopy) {
        HorizontalPartition<Key, Value>* replaced = replacePartition(hpCopy, newType, priorities, searchPriorities);
        delete(hpCopy);
        if(!replaced->isValid()) {
            std::cerr << "Invalid tree structure after changePartitionType" << std::endl;
            throw 666;
        }
        return std::optional<HorizontalPartition<Key, Value>*>{replaced};
    } 
    // change type of existing partition
    else {
        HorizontalPartition<Key, Value>* partitionToReplace;
        HorizontalPartition<Key, Value>* parentPartition;
        try {
            partitionToReplace = lookupResult.value();
            parentPartition = hpCopy->getPartitionContainingKey(key, id).value();
        }
        catch (std::bad_optional_access& e) {
            std::cerr << "changePartitionType: Bad partitions lookup: " << key << std::endl;
            throw e;
        }
        HorizontalPartition<Key, Value>* replaced = replacePartition(partitionToReplace, newType, priorities, searchPriorities);
        parentPartition->remove_partition(key);
        try {
            parentPartition->insert(key, replaced, false);
        } 
        catch (DuplicateKeyException& e) {
            std::cerr << "DuplicateKeyException in " << __FUNCTION__  << " with Key " << key << std::endl;
            throw e;
        }
        delete(partitionToReplace);
        if(!hpCopy->isValid()) {
            std::cerr << "Invalid tree structure after changePartitionType" << std::endl;
            throw 666;
        }
        return std::optional<HorizontalPartition<Key, Value>*>{hpCopy};
    }
}

/** Change search method for partition node */
template<typename Key, typename Value>
std::optional<HorizontalPartition<Key, Value>*> changePartitionSearchMethod(HorizontalPartition<Key, Value>* hp, Key key, id_type id, bool clone=true, std::string newMethod="", 
                                                                            std::map<std::string, double> priorities=std::map<std::string, double>()) {
    HorizontalPartition<Key, Value>* hpCopy = hp;
    HorizontalPartition<Key, Value>* partitionToMutate = hpCopy;
    std::optional<HorizontalPartition<Key, Value>*> lookupResult = hpCopy->lookupPartitionKey(key, id);
    if(lookupResult) {
        partitionToMutate = lookupResult.value();
    }
    // determine search method probability distribution and draw newMethod
    std::vector<std::string>* methods = partitionToMutate->getPartitionSearchMethods();
    if (newMethod.empty() || std::find(methods->begin(), methods->end(), newMethod) == methods->end()) {
        std::vector<double> probabilities;
        for (auto method : *methods) {
            if (priorities.find(method) != priorities.end())
                probabilities.push_back(priorities[method]);
            else
                probabilities.push_back(1.0);
        }
        std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
        std::mt19937 gen(std::random_device{}());
        newMethod = (*methods)[distribution(gen)];
    }
    // abort if method change was not succesful, i.e. the partition still uses the same method
    IndexSearchMethod<Key, HorizontalPartition<Key, Value>*>* formerMethod = partitionToMutate->get_partition_search_method();
    if(newMethod == toString(formerMethod->getSearchType())) {
        delete(methods);
        return std::nullopt;
    }
    // clone index if necessary
    if (clone) {
        hpCopy = hp->clone();
        partitionToMutate = hpCopy;
        lookupResult = hpCopy->lookupPartitionKey(key, id);
        if(lookupResult) {
            partitionToMutate = lookupResult.value();
        }
    }
    // replace formerMethod with newMethod
    partitionToMutate->set_partition_search_method(newMethod);
    delete(methods);
    if(!hpCopy->isValid()) {
        std::cerr << "Invalid tree structure after changePartitionSearchMethod" << std::endl;
        throw 666;
    }
    return std::optional<HorizontalPartition<Key, Value>*>{hpCopy};
}

/** Change search method for entry node */
template<typename Key, typename Value>
std::optional<HorizontalPartition<Key, Value>*> changeEntrySearchMethod(HorizontalPartition<Key, Value>* hp, Key key, id_type id, bool clone=true, std::string newMethod="",
                                                                        std::map<std::string, double> priorities=std::map<std::string, double>()) {
    HorizontalPartition<Key, Value>* hpCopy = hp;
    HorizontalPartition<Key, Value>* partitionToMutate = hpCopy;
    std::optional<HorizontalPartition<Key, Value>*> lookupResult = hpCopy->lookupPartitionKey(key, id);
    if(lookupResult) {
        partitionToMutate = lookupResult.value();
    }
    // determine search method probability distribution and draw newMethod
    std::vector<std::string>* methods = partitionToMutate->getEntrySearchMethods();
    if (newMethod.empty() || std::find(methods->begin(), methods->end(), newMethod) == methods->end()) {
        std::vector<double> probabilities;
        for (auto method : *methods) {
            if (priorities.find(method) != priorities.end())
                probabilities.push_back(priorities[method]);
            else
                probabilities.push_back(1.0);
        }
        std::discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
        std::mt19937 gen(std::random_device{}());
        newMethod = (*methods)[distribution(gen)];
    }
    // replace formerMethod with newMethod
    IndexSearchMethod<Key, Value>* formerMethod = partitionToMutate->get_entry_search_method();
    // abort if method change was not succesful, i.e. the partition still uses the same method
    if(newMethod == toString(formerMethod->getSearchType())) {
        delete(methods);
        return std::nullopt;
    }
    if (clone) {
        hpCopy = hp->clone();
        partitionToMutate = hpCopy;
        lookupResult = hpCopy->lookupPartitionKey(key, id);
        if(lookupResult) {
            partitionToMutate = lookupResult.value();
        }
    }
    partitionToMutate->set_entry_search_method(newMethod);
    delete(methods);
    if(!hpCopy->isValid()) {
        std::cerr << "Invalid tree structure after changeEntrySearchMethod" << std::endl;
        throw 666;
    }
    return std::optional<HorizontalPartition<Key, Value>*>{hpCopy};
}

/** 
 * Split partition into k child partitions 
 * The partition boundaries are determined by the given `quantiles`
 */
template<typename Key, typename Value>
std::optional<HorizontalPartition<Key, Value>*> splitPartitionHorizontally(HorizontalPartition<Key, Value>* hp, Key key, id_type id, bool clone=true, std::vector<double> quantiles = {0.5}) {
    auto hpCopy = hp;
    /* If key is present: split corresponding partition
     * If not: split root node */
    auto lookupResult = hpCopy->lookupPartitionKey(key, id);
    auto lookupResultPreviousParent = hpCopy->getPartitionContainingKey(key, id);
    HorizontalPartition<Key, Value>* partitionToSplit;
    HorizontalPartition<Key, Value>* parent;
    if (lookupResult && lookupResult.value() != hpCopy) {
        partitionToSplit = lookupResult.value();
        parent = lookupResultPreviousParent.value();
    } else {
        partitionToSplit = hpCopy;
        parent = nullptr;
    }
    if (parent) {
        unsigned int freeSlots = parent->partition_capacity() - parent->get_partitions()->size();
        if (freeSlots < quantiles.size()) {
            return std::nullopt;
        }
    }
    
    if (clone) {
        hpCopy = hp->clone();
        lookupResult = hpCopy->lookupPartitionKey(key, id);
        lookupResultPreviousParent = hpCopy->getPartitionContainingKey(key, id);
        if (lookupResult && lookupResult.value() != hpCopy) {
            partitionToSplit = lookupResult.value();
            parent = lookupResultPreviousParent.value();
        } else {
            partitionToSplit = hpCopy;
            parent = nullptr;
        }
    }

    /* Retrieve and sort partitionToSplit's keys */
    std::vector<Key>* partitionKeys = partitionToSplit->getPartitionKeys(false);
    std::sort(partitionKeys->begin(), partitionKeys->end());
    std::vector<Key>* entryKeys = partitionToSplit->getEntryKeys(false);
    std::sort(entryKeys->begin(), entryKeys->end());

    /* Create new parent node in case the root node is split (which has no parent) */
    if (!parent) {
        parent = createPartition(partitionToSplit);
        if (parent->partition_capacity() == 0)
            parent->set_partition_capacity(std::max(quantiles.size() + 1, 4ul));
        parent->insert(partitionToSplit->max_key(), partitionToSplit, false);
        hpCopy = parent;
    }

    /* Create a new partition per quantile and insert it into parent
     * redistribute partitionToSplit's child partitions within the current quantile
     * redistribute partitionToSplit's entries which belong to this new node */
    unsigned int startIndex = 0;
    bool split_performed = false;
    if (partitionToSplit->get_partitions()->size() == 0) {
        /* Splitting leaf node. */
        for (double q : quantiles) {
            unsigned int endIndex = (unsigned int)(q * partitionToSplit->get_entries()->size());
            if (endIndex == startIndex) {
                continue;
            }

            auto neighborPartition = createPartition(partitionToSplit);
            Key largestKey = (*entryKeys)[startIndex];
            for (unsigned int i=startIndex; i < endIndex; i++) {
                Key entry_to_move = (*entryKeys)[i];
                Value value_to_move = *partitionToSplit->pointQ(entry_to_move);
                partitionToSplit->remove_entry(entry_to_move);
                neighborPartition->insert(entry_to_move, value_to_move);
                if (entry_to_move > largestKey)
                    largestKey = entry_to_move;
            }
            parent->insert(largestKey, neighborPartition, false);
            startIndex = endIndex;
            split_performed = true;
        }
    } else {
        /* Splitting internal node. */
        for (double q : quantiles) {
            unsigned int endIndex = (unsigned int)(q * partitionToSplit->get_partitions()->size());
            /* if endIndex == startIndex: nothing to do, continue */
            if (endIndex == startIndex) {
                continue;
            }

            auto neighborPartition = createPartition(partitionToSplit);
            Key largestKey = (*partitionKeys)[startIndex];
            for (unsigned int i=startIndex; i < endIndex; i++) {
                auto partitionToMove = partitionToSplit->lookupPartitionKey((*partitionKeys)[i]).value();
                partitionToSplit->remove_partition((*partitionKeys)[i]);
                neighborPartition->insert((*partitionKeys)[i], partitionToMove, false);
                if ((*partitionKeys)[i] > largestKey)
                    largestKey = (*partitionKeys)[i];
            }
            for (auto it=entryKeys->begin(); it != entryKeys->end(); ) {
                if (*it > largestKey)
                    break;
                Value v = partitionToSplit->pointQ(*it).value();
                partitionToSplit->remove_entry(*it);
                neighborPartition->insert(*it, v);
                it = entryKeys->erase(it);
            }
            parent->insert(largestKey, neighborPartition, false);
            startIndex = endIndex;
            split_performed = true;
        }
    }
    delete(partitionKeys);
    delete(entryKeys);

    if (not split_performed) {
        if (clone)
            delete(hpCopy);
        return std::nullopt;
    }
    if(!hpCopy->isValid()) {
        std::cerr << "Invalid tree structure after splitPartitionHorizontally" << std::endl;
        throw 666;
    }
    return std::optional<HorizontalPartition<Key, Value>*>{hpCopy};
}

/** Insert additional partition directly above */
template<typename Key, typename Value>
std::optional<HorizontalPartition<Key, Value>*> splitPartitionVertically(HorizontalPartition<Key, Value>* hp, Key key, id_type id, bool clone=true) {
    auto hpCopy = hp;
    if (clone)
        hpCopy = hp->clone();

    auto lookupResult = hpCopy->lookupPartitionKey(key, id);
    auto lookupResultParent = hpCopy->getPartitionContainingKey(key, id);
    HorizontalPartition<Key, Value>* partitionToSplit = nullptr;
    HorizontalPartition<Key, Value>* parentPartition = nullptr;
    /* if key does not exist: split root node */
    if (!lookupResult || lookupResult.value() == hpCopy) {
        partitionToSplit = hpCopy;
        parentPartition = nullptr;
    } else {
        partitionToSplit = lookupResult.value();
        parentPartition = lookupResultParent.value();
    }
    HorizontalPartition<Key, Value>* newPartition;
    /* If parentPartition exists: Insert new partition between parentPartition and partitionToSplit */
    if (parentPartition) {
        newPartition = createPartition(parentPartition);
        parentPartition->remove_partition(key);
        parentPartition->insert(key, newPartition, false);
        newPartition->insert(key, partitionToSplit, false);
    } 
    /* Else: Insert new partition above partitionToSplit as new root */
    else {
        newPartition = createPartition(partitionToSplit);
        if (newPartition->partition_capacity() == 0)
            newPartition->set_partition_capacity(4);
        newPartition->insert(partitionToSplit->max_key(), partitionToSplit, false);
        hpCopy = newPartition;
    }
    
    if(!hpCopy->isValid()) {
        std::cerr << "Invalid tree structure after splitPartitionVertically" << std::endl;
        throw 666;
    }
    return std::optional<HorizontalPartition<Key, Value>*>{hpCopy};
}

/** merge adjacent partitions, try left neighbor first */
template<typename Key, typename Value>
std::optional<HorizontalPartition<Key, Value>*> mergePartitionHorizontally(HorizontalPartition<Key, Value>* hp, Key key, id_type id, bool clone=true) {
    auto hpCopy = hp;

    /* If key is not present or parent node is not found: abort */
    auto lookupResultMerge = hpCopy->lookupPartitionKey(key, id);
    auto lookupResultParent = hpCopy->getPartitionContainingKey(key, id);
    if (!lookupResultMerge || !lookupResultParent) {
        return std::nullopt;
    }
    auto partitionToMerge = lookupResultMerge.value();
    auto parentPartition = lookupResultParent.value();
    auto keys = parentPartition->getPartitionKeys(false);
    std::sort(keys->begin(), keys->end()); // TODO: why sort? partition keys are always sorted in our setting?
    auto itMerge = std::find(keys->begin(), keys->end(), key);
    /* If key is not found in parent node: abort
     * Should not happen as previous lookup was successful at this point */
    if (itMerge == keys->end()) {
        delete(keys);
        return std::nullopt;
    }
    HorizontalPartition<Key, Value>* neighborPartition = nullptr;
    bool right = false;
    Key neighborKey = *itMerge;
    /* If key is not the first key in parent: Merge with left neighbor */
    if (itMerge != keys->begin()) {
        neighborKey = *(itMerge - 1);
        auto lookupResultNeighbor = parentPartition->lookupPartitionKey(neighborKey);
        if (lookupResultNeighbor) {
            neighborPartition = lookupResultNeighbor.value();
            /* If partitionToMerge does not have enough capacity to additionally store all partitions and entries of neighborPartition: abort */
            unsigned int combinedPartitionSize = neighborPartition->get_partitions()->size() + partitionToMerge->get_partitions()->size();
            unsigned int combinedEntrySize = neighborPartition->get_entries()->size() + partitionToMerge->get_entries()->size();
            if (combinedPartitionSize > partitionToMerge->partition_capacity() || combinedEntrySize > partitionToMerge->entry_capacity()) {
                neighborPartition = nullptr;
            }
        }
    }
    /* If key is the first and not the last key in parent: Merge with right neighbor */
    if (!neighborPartition && itMerge != (keys->end() - 1)) {
        neighborKey = *(itMerge + 1);
        auto lookupResultNeighbor = parentPartition->lookupPartitionKey(neighborKey);
        if (lookupResultNeighbor) {
            neighborPartition = lookupResultNeighbor.value();
            right = true;
            /* If partitionToMerge does not have enough capacity to additionally store all partitions and entries of neighborPartition: abort */
            unsigned int combinedPartitionSize = neighborPartition->get_partitions()->size() + partitionToMerge->get_partitions()->size();
            unsigned int combinedEntrySize = neighborPartition->get_entries()->size() + partitionToMerge->get_entries()->size();
            if (combinedPartitionSize > partitionToMerge->partition_capacity() || combinedEntrySize > partitionToMerge->entry_capacity()) {
                neighborPartition = nullptr;
            }
        }
    }
    delete(keys);
    /* If key is the only key in parentPartition or the neighborPartition was not found: abort 
     * If partitionToMerge is a Hash partition and neighbor contains partitions: abort */
    if (!neighborPartition || (partitionToMerge->getPartitionType() == PartitionType::HashPartition && neighborPartition->get_partitions()->size() > 0)) {
        return std::nullopt;
    }
    
    if (clone) {
        hpCopy = hp->clone();
        lookupResultMerge = hpCopy->lookupPartitionKey(key, id);
        lookupResultParent = hpCopy->getPartitionContainingKey(key, id);
        partitionToMerge = lookupResultMerge.value();
        parentPartition = lookupResultParent.value();
        auto lookupResultNeighbor = parentPartition->lookupPartitionKey(neighborKey);
        neighborPartition = lookupResultNeighbor.value();
    }
    
    /* Copy over partitions */
    auto keysP = neighborPartition->getPartitionKeys(false);
    for (Key k: *keysP) {
        partitionToMerge->insert(k, neighborPartition->lookupPartitionKey(k).value(), false);
        neighborPartition->remove_partition(k);
    }
    delete(keysP);

    /* Copy over entries */
    auto keysE = neighborPartition->getEntryKeys(false);
    for (Key k: *keysE) {
        partitionToMerge->insert(k, neighborPartition->pointQ(k).value());
        neighborPartition->remove_entry(k);
    }
    delete(keysE);

    /* Remove neighborPartition from parentPartition
     * if merged with right neighbor: reinsert partitionToMerge in parentPartition with neighborKey */
    if (right) {
        parentPartition->remove_partition(key);
        parentPartition->remove_delete_partition(neighborKey);
        try {
            parentPartition->insert(neighborKey, partitionToMerge, false);
        }
        catch (DuplicateKeyException& e) {
            std::cerr << "DuplicateKeyException in " << __FUNCTION__ << " case 2: " << " with Key " << neighborKey << std::endl;
            throw e;
        }
    } else {
        parentPartition->remove_delete_partition(neighborKey);
    }
    
    if(!hpCopy->isValid()) {
        std::cerr << "Invalid tree structure after mergePartitionHorizontally" << std::endl;
        throw 666;
    }
    return std::optional<HorizontalPartition<Key, Value>*>{hpCopy};
}

/** merge (vertically) adjacent partitions */
template<typename Key, typename Value>
std::optional<HorizontalPartition<Key, Value>*> mergePartitionVertically(HorizontalPartition<Key, Value>* hp, Key key, id_type id, bool clone=true) {
    auto hpCopy = hp;

    /* If key is not present or parent node is not found: abort */
    auto lookupResultMerge = hpCopy->lookupPartitionKey(key, id);
    auto lookupResultParent = hpCopy->getPartitionContainingKey(key, id);
    if (!lookupResultMerge || !lookupResultParent) {
        return std::nullopt;
    }
    auto partitionToMerge = lookupResultMerge.value();
    auto parentPartition = lookupResultParent.value();

    auto lookupResultGrandParent = hpCopy->getPartitionContainingKey(key, parentPartition->getId(), true);
    /* For the merge to be possible, the following condition must hold:
     * 1) The parent must have exactly one child, namely partitionToMerge
     * 2) The child must have enough entry slots to store all entries of parentPartition 
     * If this condition does not hold: abort */
    unsigned int freeEntrySlotsChild = partitionToMerge->get_entries()->capacity() - partitionToMerge->get_entries()->size();
    if (parentPartition->get_partitions()->size() > 1 || freeEntrySlotsChild < parentPartition->get_entries()->size()) {
        return std::nullopt;
    }
    
    if (clone) {
        hpCopy = hp->clone();
        lookupResultMerge = hpCopy->lookupPartitionKey(key, id);
        lookupResultParent = hpCopy->getPartitionContainingKey(key, id);
        partitionToMerge = lookupResultMerge.value();
        parentPartition = lookupResultParent.value();
        lookupResultGrandParent = hpCopy->getPartitionContainingKey(key, parentPartition->getId(), true);
    }
    
    /* remove partitionToMerge from parentPartition */
    parentPartition->remove_partition(key);

    /* copy parentPartition's entries to child node */
    auto keysE = parentPartition->getEntryKeys(false);
    for (Key k: *keysE) {
        partitionToMerge->insert(k, parentPartition->pointQ(k).value());
        parentPartition->remove_entry(k);
    }
    delete(keysE);

    /* if grandParentPartition exists: remove parentPartition from grandParentPartition
     * and insert partitionToMerge at its place */
    if (lookupResultGrandParent) {
        auto grandParentPartition = lookupResultGrandParent.value();
        auto partitionKeys = grandParentPartition->getPartitionKeys(false);
        for (Key p: *partitionKeys) {
            if (grandParentPartition->lookupPartitionKey(p, parentPartition->getId()).has_value()) {
                grandParentPartition->remove_partition(p);
                grandParentPartition->insert(p, partitionToMerge, false);
                break;
            }
        }
        delete(partitionKeys);
    }
    /* otherwise: partitionToMerge becomes the new root */
    else {
        hpCopy = partitionToMerge;
    }

    /* delete now superfluous parentPartition */
    delete(parentPartition);

    if(!hpCopy->isValid()) {
        std::cerr << "Invalid tree structure after mergePartitionVertically" << std::endl;
        throw 666;
    }
    return std::optional<HorizontalPartition<Key, Value>*>{hpCopy};
}


/*===================================================================================================
 * Basic Mutations
 * ================================================================================================*/

/**
 * Perform multiple mutations in a row
 * Each mutation is chosen randomly
 * The node to perform the first mutation on is chosen randomly
 * Afterwards, each mutation adds the node it operated on as well as its child nodes to a set of valid nodes
 * future mutations are then performed based on a node randomly chosen from this set
 */
template<typename Key, typename Value>
class LocalChainMutation : public Mutation<Key, Value> {
    public:
        LocalChainMutation(double probability, Key lowerBound, Key upperBound, std::set<Mutation<Key, Value>*> mutations, bool maxFitness)
            : Mutation<Key, Value>(probability, MutationType::LocalChainMutationType),
            lowerBound(lowerBound),
            upperBound(upperBound),
            chainLength(Configuration::Get().chainMutationLength),
            mutations(mutations),
            partitionPriorities(Configuration::Get().partitionPriorities),
            searchPriorities(Configuration::Get().searchPriorities),
            gen(std::random_device{}()), //TODO: why no seed?
            maxFitness(maxFitness)
        {
            for (auto mut : mutations)
                probabilities.push_back(mut->getProbability());
            dist = std::discrete_distribution<int>(probabilities.begin(), probabilities.end());
        }

        const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) override {
            this->executionCounter++;
            auto indCasted = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
            auto mutant = indCasted->getHorizontalPartition()->clone();
            auto bestIndividual = ind;
            try {
                std::set<id_type> idsToChooseFrom;
                auto mapping = indCasted->get_stats();
                auto firstMapping = mapping;
                auto oldMapping = mapping;
                for(auto it=mapping->begin(); it!=mapping->end(); it++)
                    idsToChooseFrom.insert(it->first);
                std::set<id_type> usedIds;
                for (unsigned int i=0; i<chainLength; i++) {
                    int randomChoice = dist(gen); 
                    auto mutation = mutations.begin();
                    std::advance(mutation, randomChoice);
                    /* in the first iteration, use all available keys
                     * in all following iterations, use only keys of nodes already selected
                     * or keys of direct children of these nodes */
                    if (i > 0 && usedIds.size() > 0) {
                        oldMapping = mapping;
                        mapping = new std::map<id_type, NodeStatistics<Key, Value>*>();
                        mutant->map_statistics(mapping, true);
                        bool firstMatch = true;
                        /* include previously used ids as well as their child nodes */
                        for (auto it=usedIds.begin(); it != usedIds.end(); it++) {
                            if (mapping->find(*it) == mapping->end())
                                continue;
                            auto currentPartition = mutant->lookupPartitionKey((*mapping)[*it]->get_key(), *it);
                            if (currentPartition) {
                                if (firstMatch) {
                                    idsToChooseFrom.clear();
                                    firstMatch = false;
                                }
                                idsToChooseFrom.insert(*it);
                                for(auto it2=currentPartition.value()->get_partitions()->value_begin(), it2_end=currentPartition.value()->get_partitions()->value_end(); it2 != it2_end; ++it2) {
                                    idsToChooseFrom.insert((*it2)->getId());
                                }
                            }
                        }
                        /* include ids introduced by previous mutation */
                        for (auto it=mapping->begin(), end=mapping->end(); it != end; ++it) {
                            if (oldMapping->find(it->first) == oldMapping->end()) {
                                if (firstMatch) {
                                    idsToChooseFrom.clear();
                                    firstMatch = false;
                                }
                                idsToChooseFrom.insert(it->first);
                            }
                        }
                        /* include all IDs if no matches so far */
                        if (firstMatch) {
                            idsToChooseFrom.clear();
                            for(auto it=mapping->begin(); it!=mapping->end(); it++)
                                idsToChooseFrom.insert(it->first);
                        }
                    }
                    id_type id = getRandomId(mapping, &idsToChooseFrom);
                    Key key = (*(mapping))[id]->get_key();
                    if (i>0 && usedIds.size() > 0)
                        if (oldMapping != firstMapping)
                            delete(oldMapping);
                    std::optional<HorizontalPartition<Key, Value>*> mutationResult;
                    switch((*mutation)->getType()) {
                        case MutationType::PartitionTypeChangeMutationType:
                            mutationResult = changePartitionType<Key, Value>(mutant, key, id, false, PartitionType::UnknownPartition, partitionPriorities, searchPriorities);
                            break;
                        case MutationType::PartitionSearchChangeMutationType:
                            mutationResult = changePartitionSearchMethod<Key, Value>(mutant, key, id, false, "", searchPriorities);
                            break;
                        case MutationType::EntrySearchChangeMutationType:
                            mutationResult = changeEntrySearchMethod<Key, Value>(mutant, key, id, false, "", searchPriorities);
                            break;
                        case MutationType::HorizontalMergeMutationType:
                            mutationResult = mergePartitionHorizontally<Key, Value>(mutant, key, id, false);
                            break;
                        case MutationType::VerticalMergeMutationType:
                            mutationResult = mergePartitionVertically<Key, Value>(mutant, key, id, false);
                            break;
                        case MutationType::HorizontalSplitMutationType:
                            mutationResult = splitPartitionHorizontally<Key, Value>(mutant, key, id, false);
                            break;
                        case MutationType::VerticalSplitMutationType:
                            mutationResult = splitPartitionVertically<Key, Value>(mutant, key, id, false);
                            break;
                        default:
                            std::cerr << "Mutation: " << (*mutation)->getName() << " not yet implemented in LocalChainMutation" << std::endl;
                            throw 666;
                    }
                    if (mutationResult) {
                        mutant = mutationResult.value();
                        usedIds.insert(id);
                        auto candidate = new HorizontalPartitionIndividual<Key, Value>(mutant->clone(), indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                             indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                             indCasted->getID(), indCasted->getMissingEntries());
                        if (maxFitness && candidate->getFitness() > bestIndividual->getFitness()) {
                            if(bestIndividual != ind)
                                delete(bestIndividual);
                            bestIndividual = candidate;
                        } else if (!maxFitness && candidate->getFitness() < bestIndividual->getFitness()) {
                            if (bestIndividual != ind)
                                delete(bestIndividual);
                            bestIndividual = candidate;
                        } else {
                            delete(candidate);
                        }
                    }
                }
                if (mapping != firstMapping)
                    delete(mapping);
                delete(mutant);
                return bestIndividual;
            } catch(std::exception& e) {
                std::cerr << "Exception in LocalChainMutation" << std::endl;
                auto candidate = new HorizontalPartitionIndividual<Key, Value>(mutant, indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                     indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                     indCasted->getID(), indCasted->getMissingEntries());
                candidate->toDotFile("ExceptionBefore.dot", true, true);
                throw e;
            }
        }

    private:
        Key lowerBound, upperBound;
        unsigned int chainLength;
        std::set<Mutation<Key, Value>*> mutations;
        std::map<PartitionType, double> partitionPriorities;
        std::map<std::string, double> searchPriorities;
        std::mt19937 gen;
        bool maxFitness;
        std::vector<double> probabilities;
        std::discrete_distribution<int> dist;
};

/**
 * change the type of a randomly chosen node 
 */
template<typename Key, typename Value>
class PartitionTypeChangeMutation : public Mutation<Key, Value> {
    public:
        PartitionTypeChangeMutation(double probability, Key lowerBound, Key upperBound)
            : Mutation<Key, Value>(probability, MutationType::PartitionTypeChangeMutationType),
            lowerBound(lowerBound),
            upperBound(upperBound),
            priorities(Configuration::Get().partitionPriorities),
            searchPriorities(Configuration::Get().searchPriorities)
        {}

        const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) override {
            Key key;
            id_type id;
            try {
                this->executionCounter++;
                auto indCasted = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
                id = getRandomId(indCasted->get_stats());
                key = (*(indCasted->get_stats()))[id]->get_key();
                auto mutationResult = changePartitionType<Key, Value>(indCasted->getHorizontalPartition(), key, id, true, PartitionType::UnknownPartition, priorities, searchPriorities);
                if (mutationResult)
                    return new HorizontalPartitionIndividual<Key, Value>(mutationResult.value(), indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                             indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                             indCasted->getID(), indCasted->getMissingEntries());
                else
                    return ind;
            } catch(std::exception& e) {
                std::cerr << "Exception in PartitionTypeChangeMutation with key: " << key << " and Id: " << id << std::endl;
                ind->toDotFile("ExceptionBefore.dot", true, true);
                throw e;
            }
        }

    private:
        Key lowerBound, upperBound;
        std::map<PartitionType, double> priorities;
        std::map<std::string, double> searchPriorities;
};

/**
 * change the search method for the partitions in a randomly chosen node 
 */
template<typename Key, typename Value>
class PartitionSearchChangeMutation : public Mutation<Key, Value> {
    public:
        PartitionSearchChangeMutation(double probability, Key lowerBound, Key upperBound)
            : Mutation<Key, Value>(probability, MutationType::PartitionSearchChangeMutationType),
            lowerBound(lowerBound),
            upperBound(upperBound),
            priorities(Configuration::Get().searchPriorities)
        {}

        const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) override {
            Key key;
            id_type id;
            try {
                this->executionCounter++;
                auto indCasted = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
                id = getRandomId(indCasted->get_stats());
                key = (*(indCasted->get_stats()))[id]->get_key();
                auto mutationResult = changePartitionSearchMethod<Key, Value>(indCasted->getHorizontalPartition(), key, id, true, "", priorities);
                if (mutationResult)
                    return new HorizontalPartitionIndividual<Key, Value>(mutationResult.value(), indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                             indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                             indCasted->getID(), indCasted->getMissingEntries());
                else
                    return ind;
            } catch(std::exception& e) {
                std::cerr << "Exception in PartitionSearchChangeMutation with key: " << key << " and Id: " << id << std::endl;
                ind->toDotFile("ExceptionBefore.dot", true, true);
                throw e;
            }
        }

    private:
        Key lowerBound, upperBound;
        std::map<std::string, double> priorities;
};

/**
 * change the search method for the entries in a randomly chosen node 
 */
template<typename Key, typename Value>
class EntrySearchChangeMutation : public Mutation<Key, Value> {
    public:
        EntrySearchChangeMutation(double probability, Key lowerBound, Key upperBound)
            : Mutation<Key, Value>(probability, MutationType::EntrySearchChangeMutationType),
            lowerBound(lowerBound),
            upperBound(upperBound),
            priorities(Configuration::Get().searchPriorities)
        {}

        const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) override {
            Key key;
            id_type id;
            try {
                this->executionCounter++;
                auto indCasted = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
                id = getRandomId(indCasted->get_stats());
                key = (*(indCasted->get_stats()))[id]->get_key();
                auto mutationResult = changeEntrySearchMethod<Key, Value>(indCasted->getHorizontalPartition(), key, id, true, "", priorities);
                if (mutationResult)
                    return new HorizontalPartitionIndividual<Key, Value>(mutationResult.value(), indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                             indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                             indCasted->getID(), indCasted->getMissingEntries());
                else
                    return ind;
            } catch(std::exception& e) {
                std::cerr << "Exception in EntrySearchChangeMutation with key: " << key << " and Id: " << id << std::endl;
                ind->toDotFile("ExceptionBefore.dot", true, true);
                throw e;
            }
        }

    private:
        Key lowerBound, upperBound;
        std::map<std::string, double> priorities;
};

/**
 * split a randomly chosen node horizontally,
 * i.e. the its child partitions are first redistributed according to the given quantiles
 * in a second step, the entries of the node are then also redistributed to ensure a valid tree structure
 */
template<typename Key, typename Value>
class HorizontalSplitMutation : public Mutation<Key, Value> {
    public:
        HorizontalSplitMutation(double probability, Key lowerBound, Key upperBound) : Mutation<Key, Value>(probability, MutationType::HorizontalSplitMutationType), lowerBound(lowerBound), upperBound(upperBound) {}

        const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) override {
            Key key;
            id_type id;
            try {
                this->executionCounter++;
                auto indCasted = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
                id = getRandomId(indCasted->get_stats());
                key = (*(indCasted->get_stats()))[id]->get_key();
                auto mutationResult = splitPartitionHorizontally<Key, Value>(indCasted->getHorizontalPartition(), key, id);
                if (mutationResult)
                    return new HorizontalPartitionIndividual<Key, Value>(mutationResult.value(), indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                             indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                             indCasted->getID(), indCasted->getMissingEntries());
                else
                    return ind;
            } catch(std::exception& e) {
                std::cerr << "Exception in HorizontalSplitMutation with key: " << key << " and Id: " << id << std::endl;
                ind->toDotFile("ExceptionBefore.dot", true, true);
                throw e;
            }
        }

    private:
        Key lowerBound, upperBound;
};

/**
 * perform a vertical split on a randonly chosen node,
 * i.e. insert a new parent node above the chosen node, possibly creating a new root
 */
template<typename Key, typename Value>
class VerticalSplitMutation : public Mutation<Key, Value> {
    public:
        VerticalSplitMutation(double probability, Key lowerBound, Key upperBound) : Mutation<Key, Value>(probability, MutationType::VerticalSplitMutationType), lowerBound(lowerBound), upperBound(upperBound) {}

        const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) override {
            Key key;
            id_type id;
            try {
                this->executionCounter++;
                auto indCasted = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
                id = getRandomId(indCasted->get_stats());
                key = (*(indCasted->get_stats()))[id]->get_key();
                auto mutationResult = splitPartitionVertically<Key, Value>(indCasted->getHorizontalPartition(), key, id);
                if (mutationResult)
                    return new HorizontalPartitionIndividual<Key, Value>(mutationResult.value(), indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                             indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                             indCasted->getID(), indCasted->getMissingEntries());
                else
                    return ind;
            } catch(std::exception& e) {
                std::cerr << "Exception in VerticalSplitMutation with key: " << key << " and Id: " << id << std::endl;
                ind->toDotFile("ExceptionBefore.dot", true, true);
                throw e;
            }
        }

    private:
        Key lowerBound, upperBound;
};

/**
 * perform a horizontal merge of a given node,
 * i.e. try to shift all child partitions as well as entries of a neighbor node to the chosen node
 * the mutation prefers the left neighbor node
 * if there is no direct, left neighbor, a merge with the right neighbor is performed
 */
template<typename Key, typename Value>
class HorizontalMergeMutation : public Mutation<Key, Value> {
    public:
        HorizontalMergeMutation(double probability, Key upperBound) : Mutation<Key, Value>(probability, MutationType::HorizontalMergeMutationType), upperBound(upperBound) {}

        const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) override {
            Key key;
            id_type id;
            try {
                this->executionCounter++;
                auto indCasted = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
                id = getRandomId(indCasted->get_stats());
                key =(*(indCasted->get_stats()))[id]->get_key();
                auto mutationResult = mergePartitionHorizontally<Key, Value>(indCasted->getHorizontalPartition(), key, id);
                if (mutationResult)
                    return new HorizontalPartitionIndividual<Key, Value>(mutationResult.value(), indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                             indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                             indCasted->getID(), indCasted->getMissingEntries());
                else
                    return ind;
            } catch(std::exception& e) {
                std::cerr << "Exception in HorizontalMergeMutation with key: " << key << " and Id: " << id << std::endl;
                ind->toDotFile("ExceptionBefore.dot", true, true);
                throw e;
            }
        }

    private:
        Key upperBound;
};

/**
 * perform a vertical merge on a randomly chosen node,
 * i.e. try to shift all values of the chosen node (child partitions as well as entries) to the parent node
 */
template<typename Key, typename Value>
class VerticalMergeMutation : public Mutation<Key, Value> {
    public:
        VerticalMergeMutation(double probability, Key upperBound) : Mutation<Key, Value>(probability, MutationType::VerticalMergeMutationType), upperBound(upperBound) {}

        const Individual<Key, Value>* mutate(const Individual<Key, Value>* ind) override {
            Key key;
            id_type id;
            try {
                this->executionCounter++;
                auto indCasted = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
                id = getRandomId(indCasted->get_stats());
                key = (*(indCasted->get_stats()))[id]->get_key();
                auto mutationResult = mergePartitionVertically<Key, Value>(indCasted->getHorizontalPartition(), key, id);
                if (mutationResult)
                    return new HorizontalPartitionIndividual<Key, Value>(mutationResult.value(), indCasted->getWorkload(), indCasted->getPartitionCapacity(), indCasted->getEntryCapacity(),
                                                             indCasted->getPenaltyFactorMissingKeys(), indCasted->getPenaltyFactorEmptyPartitions(), indCasted->getRepetitions(), 
                                                             indCasted->getID(), indCasted->getMissingEntries());
                else
                    return ind;
            } catch(std::exception& e) {
                std::cerr << "Exception in VerticalMergeMutation with key: " << key << " and Id: " << id << std::endl;
                ind->toDotFile("ExceptionBefore.dot", true, true);
                throw e;
            }
        }

    private:
        Key upperBound;
};
