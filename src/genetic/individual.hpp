#pragma once

#include <set>
#include <unordered_set>
#include "hp/Partition.hpp"
#include "util/hash.hpp"
#include "util/Dataset.hpp"
#include "util/Workload.hpp"
#include "util/enum.hpp"
#include <ostream>
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <limits>

#define IMPORT(WHAT) using WHAT = typename Base::WHAT

template<typename Key, typename Value>
struct IndividualIdComparator;

/**
 * abstract base class for individuals
 */
template<typename Key, typename Value>
class Individual {
    public:
    using key_type = Key;
    using mapped_type = Value;

    private:
        static unsigned int ID;

    public:
        Individual() : fitness(-1.0), currentID(ID++) {}

        virtual ~Individual() {}

        virtual double getFitness() const = 0;

        virtual bool operator== (const Individual<Key, Value>& other) const = 0;
        virtual bool operator!= (const Individual<Key, Value>& other) const { return not (*this == other); }

        virtual void print(std::ostream& stream) const = 0;
        std::string to_string(Individual<Key, Value>& ind) const {
            std::ostringstream oss;
            oss << ind;
            return oss.str();
        }

        virtual void to_graphviz(std::ostream &out, const std::string& name, bool includePartitions=false, bool includeEntries=false, const std::string &annotation="", bool after_image=false) const = 0;
        virtual void toDotFile(const std::string& filename="", bool includePartitions=false, bool includeEntries=false, const std::string &annotation="", bool after_image=false) const = 0;

        virtual HorizontalPartition<Key, Value>* getPartition() const = 0;

        unsigned int getID() const { return currentID; }

        double get_stddev_fitness() const { return stddev_fitness; }

        virtual std::string describe() const = 0;
        virtual std::string describe(std::ofstream& csv, char separator=',') const = 0;

    protected:
        double stddev_fitness;
        double fitness;
        unsigned int currentID;
};

/**
 * implementation of a fitness result storing all information gathered during the computation of the fitness value
 */
template<typename Key, typename Value>
class FitnessResult {
    
    /* =================================================================================================================
     * Constructors 
     * ===============================================================================================================*/
    
    public:
    FitnessResult(double fitness, double stddev_fitness, std::unordered_set<std::pair<Key, Value>, pair_hash> erase, 
            std::unordered_set<std::pair<Key, Value>, pair_hash> include, HorizontalPartition<Key, Value>* after)
        : fitness(fitness),
        stddev_fitness(stddev_fitness),
        erase(erase),
        include(include),
        after(after)
    {}

    /* =================================================================================================================
     * Class members 
     * ===============================================================================================================*/
    double fitness;
    double stddev_fitness;
    std::unordered_set<std::pair<Key, Value>, pair_hash> erase;
    std::unordered_set<std::pair<Key, Value>, pair_hash> include;
    HorizontalPartition<Key, Value>* after;
    
    /* =================================================================================================================
     * Getter & Setter 
     * ===============================================================================================================*/
    double get_fitness() const { return fitness; }
    double get_stddev_fitness() const { return stddev_fitness; }
    std::unordered_set<std::pair<Key, Value>, pair_hash> get_entries_to_erase() { return erase; }
    std::unordered_set<std::pair<Key, Value>, pair_hash> get_entries_to_include() { return include; }
    HorizontalPartition<Key, Value>* get_after_image() { return after; }
};

/**
 * concrete implementation for an individual storing a HorizontalPartition as tree structure
 */
template<typename Key, typename Value>
class HorizontalPartitionIndividual : public Individual<Key, Value> {

    using Base = Individual<Key, Value>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    using id_type = unsigned long;

    public:
    /** Use if individual has already been initialized and only needs to maintain a set of missing entries */
    HorizontalPartitionIndividual(HorizontalPartition<Key, Value>* partition,
                                  Workload<Key, Value>* workload,
                                  size_t partition_capacity=4, size_t entry_capacity=4, double penaltyFactorMissingKeys=1000.0,
                                  double penaltyFactorEmptyPartitions=0.0,
                                  unsigned int repetitions=1, unsigned int predecessorID=0, 
                                  std::unordered_set<std::pair<Key, Value>, pair_hash> missingEntries=std::unordered_set<std::pair<Key, Value>, pair_hash>())
        : Individual<Key, Value>(),
          hp(partition),
          hp_after(nullptr),
          dataset(0),
          workload(workload),
          partition_capacity(partition_capacity),
          entry_capacity(entry_capacity),
          penaltyFactorMissingKeys(penaltyFactorMissingKeys),
          penaltyFactorEmptyPartitions(penaltyFactorEmptyPartitions),
          repetitions(repetitions),
          predecessorID(predecessorID),
          missingEntries(missingEntries),
          id_to_stats(new std::map<id_type, NodeStatistics<Key, Value>*>())
    {
        auto eval = calculateFitness(workload);
        this->fitness = eval.get_fitness();
        this->stddev_fitness = eval.get_stddev_fitness();
        for (auto e : eval.get_entries_to_erase())
            missingEntries.erase(e);
        for (auto e : eval.get_entries_to_include())
            missingEntries.insert(e);
        hp_after = eval.get_after_image();
    }

    /** 
     * Use if individual has not yet been initialized
     * All entries contained in dataset will be inserted
     * In case of CapacityException, entry will be appended to set of missing entries
     */
    HorizontalPartitionIndividual(HorizontalPartition<Key, Value>* partition,
                                  Dataset<Key, Value>* dataset, Workload<Key, Value>* workload,
                                  size_t partition_capacity=4, size_t entry_capacity=4, double penaltyFactorMissingKeys=1000.0,
                                  double penaltyFactorEmptyPartitions=0.0,
                                  unsigned int repetitions=1, unsigned int predecessorID=0) 
        : Individual<Key, Value>(),
          hp(partition),
          hp_after(nullptr),
          dataset(dataset),
          workload(workload),
          partition_capacity(partition_capacity),
          entry_capacity(entry_capacity),
          penaltyFactorMissingKeys(penaltyFactorMissingKeys),
          penaltyFactorEmptyPartitions(penaltyFactorEmptyPartitions),
          repetitions(repetitions),
          predecessorID(predecessorID),
          missingEntries(std::unordered_set<std::pair<Key, Value>, pair_hash>()),
          id_to_stats(new std::map<id_type, NodeStatistics<Key, Value>*>())
    {
        for(auto it = dataset->begin(); it != dataset->end(); ++it) {
            try {
                auto result = hp->insert(it->first, it->second);
                if(result.second.has_value())
                    hp = result.second.value();
            } catch(CapacityException& e) {
                missingEntries.insert(*it);
            }
        }
        auto eval = calculateFitness(workload);
        this->fitness = eval.get_fitness();
        this->stddev_fitness = eval.get_stddev_fitness();
        for (auto e : eval.get_entries_to_erase())
            missingEntries.erase(e);
        for (auto e : eval.get_entries_to_include())
            missingEntries.insert(e);
        hp_after = eval.get_after_image();
    }

    ~HorizontalPartitionIndividual() {
        if (hp_after and hp_after != hp) {
            delete(hp_after);
        }
        delete(hp);
        clear_id_stats(false);
        delete(id_to_stats);
    }

    /** return saved fitness value */
    double getFitness() const override { return this->fitness; }

    /** test for equality, i.e. equal horizontal partition */
    bool operator== (const Individual<Key, Value>& other) const override {
        auto casted = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(&other);
        if (not casted)
            return false;
        return *hp == *casted->hp;
    };

    private:

    /* =================================================================================================================
     * Class members 
     * ===============================================================================================================*/

    HorizontalPartition<Key, Value>* hp;
    HorizontalPartition<Key, Value>* hp_after;
    Dataset<Key, Value>* dataset;
    Workload<Key, Value>* workload;
    size_t partition_capacity;
    size_t entry_capacity;
    double penaltyFactorMissingKeys, penaltyFactorEmptyPartitions;
    unsigned int repetitions;
    const unsigned int predecessorID;
    std::unordered_set<std::pair<Key, Value>, pair_hash> missingEntries;
    std::map<id_type, NodeStatistics<Key, Value>*>* id_to_stats;
    
    /* =================================================================================================================
     * Node statistics management 
     * ===============================================================================================================*/

    /** insert or update the node statistics for a given id */
    void insert_update_id_stats(id_type id, NodeStatistics<Key, Value>* ns) {
        auto it = id_to_stats->find(id);
        if (it == id_to_stats->end())
            id_to_stats->insert(std::make_pair(id, ns));
        else {
            delete(it->second);
            it->second = ns;
        }
    }

    /** clear all entries in the id_to_stats mapping */
    void clear_id_stats(bool clone=false) const {
        if (clone) {
            for (auto it=id_to_stats->begin(), end=id_to_stats->end(); it != end; ++it) {
                delete(it->second);
            }
        }
        id_to_stats->clear();
    }


    /* =================================================================================================================
     * Fitness Function 
     * ===============================================================================================================*/

    public:
    /** 
     * calculate the fitness and return a FitnessResult
     */
    FitnessResult<Key, Value> calculateFitness(Workload<Key, Value> *workload) const {
        // initialize variables
        double fit = 0.0;
        double penalty = 0.0;
        std::unordered_set<std::pair<Key, Value>, pair_hash> entriesToErase;
        std::unordered_set<std::pair<Key, Value>, pair_hash> entriesToInsert;
        // try to insert yet missing entries
        for (auto it = missingEntries.begin(); it != missingEntries.end(); ) {
            try {
                hp->insert(it->first, it->second);
                entriesToErase.insert(*it);
            }
            catch(CapacityException& e) { it++; }
        }
        // count empty partitions and apply penalty
        unsigned int emptyPartitions = hp->countEmptyPartitions();
        fit = emptyPartitions * penaltyFactorEmptyPartitions + missingEntries.size() * penaltyFactorMissingKeys;

        // reset the NodeStatistics of the index recursively
        hp->reset_statistics(true);

        HorizontalPartition<Key, Value>* hp_after = nullptr;
        HorizontalPartition<Key, Value>* hp_copy = nullptr;
        std::set<id_type>* ids = new std::set<id_type>();
        hp->get_ids(ids);
        // perform measurement repetitions to fill a vector of durations
        // median value of these measurements is taken as fitness
        std::vector<long double> durations;
        double sum_durations = 0.0;
        for(unsigned int r=0; r < repetitions; r++) {
            uint64_t result = 0;
            if (workload->has_inserts())
                hp_copy = hp->clone();
            else
                hp_copy = hp;
            auto start = std::chrono::high_resolution_clock::now();
            // iterate over workload and execute each query
            penalty = 0.0;
            for (auto it=workload->begin(); it!=workload->end(); ++it) {
                switch((*it)->workload_kind()) {
                    case Workload_kind::Point:
                        {
                            auto p = static_cast<PointQuery<key_type, mapped_type>*>(*it);
                            result += *(hp_copy->pointQ(p->key()));
                            break;
                        }
                    case Workload_kind::Range:
                        {
                            auto r = static_cast<RangeQuery<key_type, mapped_type>*>(*it);
                            hp_copy->rangeQ(result, r->lower(), r->upper());
                            break;
                        }
                    case Workload_kind::Insert:
                        {
                            auto i = static_cast<Insertion<key_type, mapped_type>*>(*it);
                            auto res = hp_copy->insert(i->key(), i->value());
                            if (res.second)
                                hp_copy = res.second.value();
                            if (not res.first)
                                penalty += penaltyFactorMissingKeys;
                            break;
                        }
                    default:
                        throw std::invalid_argument("Unknown Workload_type");
                    }
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = (std::chrono::duration<long double, std::milli>(end - start)).count();

            sum_durations += duration;
            durations.push_back(duration);
            
            // collect the NodeStatistics generated during workload execution
            clear_id_stats();
            hp_copy->map_statistics(id_to_stats, true, ids);

            // delete all but last cloned HPs
            if (hp_after and workload->has_inserts()) {
                delete(hp_after);
            }
            hp_after = hp_copy;
        }
        delete(ids);
        if(repetitions > 1)
            hp_copy->normalize_repetitions(repetitions);

        std::sort(durations.begin(), durations.end());
        long double medianDuration = repetitions % 2 == 0 ? (durations[repetitions / 2 - 1] + durations[repetitions / 2]) / 2 : durations[repetitions / 2];
        // Compute stddev of fitness value
        double mean = sum_durations / durations.size();
        double sum_nominator = 0.0;
        for (auto d : durations) sum_nominator += (d - mean) * (d - mean);
        double stddev = sqrt((sum_nominator) / durations.size());
        fit += medianDuration + penalty;

        return FitnessResult<Key, Value>(fit, stddev, entriesToErase, entriesToInsert, hp_after);
    }

    /* =================================================================================================================
     * IO Functions 
     * ===============================================================================================================*/

    void print(std::ostream& stream) const override {
        stream << "HorizontalPartitionIndividual:" << std::endl << describe() << std::endl;
    }

    /** generate a graphviz dot string and append it to the given output stream */
    void to_graphviz(std::ostream &out,const std::string& name="index", bool includePartitions=false, bool includeEntries=false, const std::string &annotation="", bool after_image=false) const override {
        out << "digraph " << name << " {\n";
        out << "node [shape=plain, height=.1];\n";
        if (after_image)
            hp_after->dot(out, includePartitions, includeEntries);
        else    
            hp->dot(out, includePartitions, includeEntries);
        /* Annotations. */
        out << "labelloc=\"t\";";
        out << "label=\"" << annotation << "\";\n";
        out << "}\n";
    }

    /** generate a graphviz dot string and write it to the given file */
    void toDotFile(const std::string& filename="", bool includePartitions=false, bool includeEntries=false, const std::string &annotation="", bool after_image=false) const override {
        std::string fname;
        if(filename.empty())
            fname = "Individual"+std::to_string(Base::currentID)+".dot";
        else
            fname = filename;
        std::ofstream out(fname);
        if (not out.is_open()) {
            std::cerr << "COULD NOT OPEN DOTFILE FILE" << std::endl;
        }
        to_graphviz(out, "index", includePartitions, includeEntries, annotation, after_image);
        out.close();
    }
    
    /** 
     * describe the HorizontalPartition of this individual and return it as string
     * description contains information per layer including node types, search methods etc.
     */
    std::string describe() const override {
        auto nodesPerLevel = new std::vector<std::map<PartitionType, unsigned int>>();
        auto partitionSearchesPerLevel = new std::vector<std::map<SearchType, unsigned int>>();
        auto entrySearchesPerLevel = new std::vector<std::map<SearchType, unsigned int>>();
        auto emptyPartitionsPerLevel = new std::vector<unsigned int>();
        stats(hp, 0, nodesPerLevel, partitionSearchesPerLevel, entrySearchesPerLevel, emptyPartitionsPerLevel);
        return describeString(nodesPerLevel, partitionSearchesPerLevel, entrySearchesPerLevel, emptyPartitionsPerLevel, true);
    }

    /** describe the HorizontalPartition of this individual, return it as string and write it to the given output stream */
    std::string describe(std::ofstream& csv, char separator=',') const override {
        auto nodesPerLevel = new std::vector<std::map<PartitionType, unsigned int>>();
        auto partitionSearchesPerLevel = new std::vector<std::map<SearchType, unsigned int>>();
        auto entrySearchesPerLevel = new std::vector<std::map<SearchType, unsigned int>>();
        auto emptyPartitionsPerLevel = new std::vector<unsigned int>();
        stats(hp, 0, nodesPerLevel, partitionSearchesPerLevel, entrySearchesPerLevel, emptyPartitionsPerLevel);
        describeCSV(nodesPerLevel, partitionSearchesPerLevel, entrySearchesPerLevel, emptyPartitionsPerLevel, csv, false, separator);
        return describeString(nodesPerLevel, partitionSearchesPerLevel, entrySearchesPerLevel, emptyPartitionsPerLevel, true);
    }

    /** describe the HorizontalPartition of this individual based on the given maps and vectors and return it as string
     * each entry in the vectors represents one layer of the horizontal partitionSearchesPerLevel
     * each map then contains a mapping from node type / search method to number of occurences
     */
    std::string describeString(std::vector<std::map<PartitionType, unsigned int>>* nodesPerLevel,
                 std::vector<std::map<SearchType, unsigned int>>* partitionSearchesPerLevel,
                 std::vector<std::map<SearchType, unsigned int>>* entrySearchesPerLevel,
                 std::vector<unsigned int>* emptyPartitionsPerLevel, bool deleteVectors=false) const {
        std::stringstream ss;
        ss << "Individual ID: " << Base::currentID << std::endl;
        ss << "Fitness: " << Base::fitness << std::endl;
        ss << "Missing Entries: " << missingEntries.size() << std::endl;
        for (unsigned int l=0; l < nodesPerLevel->size(); l++) {
            ss << "Level: " << l << std::endl;
            ss << "\tPartition Types:" << std::endl;
            ss << "\t\tSorted Array Partitions: " << (*nodesPerLevel)[l][PartitionType::SoAPartition] << std::endl;
            ss << "\t\tTree Partitions: " << (*nodesPerLevel)[l][PartitionType::TreePartition] << std::endl;
            ss << "\t\tHash Partitions: " << (*nodesPerLevel)[l][PartitionType::HashPartition] << std::endl;
            ss << "\tEmpty Partitions: " << (*emptyPartitionsPerLevel)[l] << std::endl;
            ss << "\tPartition Search Methods:" << std::endl;
            ss << "\t\tDefault Search: " << (*partitionSearchesPerLevel)[l][SearchType::DefaultSearchMethod] << std::endl;
            ss << "\t\tLinear Search: " << (*partitionSearchesPerLevel)[l][SearchType::LinearSearchMethod] << std::endl;
            ss << "\t\tBinary Search: " << (*partitionSearchesPerLevel)[l][SearchType::BinarySearchMethod] << std::endl;
            ss << "\t\tInterpolation Search: " << (*partitionSearchesPerLevel)[l][SearchType::InterpolationSearchMethod] << std::endl;
            ss << "\t\tExponential Search: " << (*partitionSearchesPerLevel)[l][SearchType::ExponentialSearchMethod] << std::endl;
            ss << "\t\tLinear Regression Search: " << (*partitionSearchesPerLevel)[l][SearchType::LinearRegressionSearchMethod] << std::endl;
            ss << "\tEntry Search Methods:" << std::endl;
            ss << "\t\tDefault Search: " << (*entrySearchesPerLevel)[l][SearchType::DefaultSearchMethod] << std::endl;
            ss << "\t\tLinear Search: " << (*entrySearchesPerLevel)[l][SearchType::LinearSearchMethod] << std::endl;
            ss << "\t\tBinary Search: " << (*entrySearchesPerLevel)[l][SearchType::BinarySearchMethod] << std::endl;
            ss << "\t\tInterpolation Search: " << (*entrySearchesPerLevel)[l][SearchType::InterpolationSearchMethod] << std::endl;
            ss << "\t\tExponential Search: " << (*entrySearchesPerLevel)[l][SearchType::ExponentialSearchMethod] << std::endl;
            ss << "\t\tLinear Regression Search: " << (*entrySearchesPerLevel)[l][SearchType::LinearRegressionSearchMethod] << std::endl;
        }
        if (deleteVectors) {
            delete(nodesPerLevel);
            delete(partitionSearchesPerLevel);
            delete(entrySearchesPerLevel);
            delete(emptyPartitionsPerLevel);
        }
        return ss.str();
    }

    /** describe the HorizontalPartition of this individual based on the given maps and vectors and write it to the given output stream in csv format */
    void describeCSV(std::vector<std::map<PartitionType, unsigned int>>* nodesPerLevel,
        std::vector<std::map<SearchType, unsigned int>>* partitionSearchesPerLevel,
        std::vector<std::map<SearchType, unsigned int>>* entrySearchesPerLevel,
        std::vector<unsigned int>* emptyPartitionsPerLevel, 
        std::ofstream& csv, bool deleteVectors=false, char separator=',') const {

        csv << Base::currentID << separator << 0 << separator << "Fitness" << separator << Base::fitness << std::endl;
        csv << Base::currentID << separator << 0 << separator << "MissingEntries" << separator << missingEntries.size() << std::endl;
        for (unsigned int l=0; l < nodesPerLevel->size(); l++) {
            csv << Base::currentID << separator << l << separator << "SoAPartitions" << separator << (*nodesPerLevel)[l][PartitionType::SoAPartition] << std::endl;
            csv << Base::currentID << separator << l << separator << "TreePartitions" << separator << (*nodesPerLevel)[l][PartitionType::TreePartition] << std::endl;
            csv << Base::currentID << separator << l << separator << "HashPartitions" << separator << (*nodesPerLevel)[l][PartitionType::HashPartition] << std::endl;
            csv << Base::currentID << separator << l << separator << "EmptyPartitions" << separator << (*emptyPartitionsPerLevel)[l] << std::endl;
            csv << Base::currentID << separator << l << separator << "Partition::DefaultSearch" << separator << (*partitionSearchesPerLevel)[l][SearchType::DefaultSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Partition::LinearSearch" << separator << (*partitionSearchesPerLevel)[l][SearchType::LinearSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Partition::BinarySearch" << separator <<(*partitionSearchesPerLevel)[l][SearchType::BinarySearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Partition::InterpolationSearch" << separator << (*partitionSearchesPerLevel)[l][SearchType::InterpolationSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Partition::ExponentialSearch" << separator << (*partitionSearchesPerLevel)[l][SearchType::ExponentialSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Partition::LinearRegressionSearch" << separator << (*partitionSearchesPerLevel)[l][SearchType::LinearRegressionSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Entry::DefaultSearch" << separator << (*entrySearchesPerLevel)[l][SearchType::DefaultSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Entry::LinearSearch" << separator << (*entrySearchesPerLevel)[l][SearchType::LinearSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Entry::BinarySearch" << separator << (*entrySearchesPerLevel)[l][SearchType::BinarySearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Entry::InterpolationSearch" << separator << (*entrySearchesPerLevel)[l][SearchType::InterpolationSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Entry::ExponentialSearch" << separator << (*entrySearchesPerLevel)[l][SearchType::ExponentialSearchMethod] << std::endl;
            csv << Base::currentID << separator << l << separator << "Entry::LinearRegressionSearch" << separator << (*entrySearchesPerLevel)[l][SearchType::LinearRegressionSearchMethod] << std::endl;
        }
        if (deleteVectors) {
            delete(nodesPerLevel);
            delete(partitionSearchesPerLevel);
            delete(entrySearchesPerLevel);
            delete(emptyPartitionsPerLevel);
        }
    }

    /** compute statistical properties of the HorizontalPartition stored in this individual and write it to the given maps and vectors */
    void stats(HorizontalPartition<Key, Value>* hp, unsigned int level, std::vector<std::map<PartitionType, unsigned int>>* nodesPerLevel,
        std::vector<std::map<SearchType, unsigned int>>* partitionSearchesPerLevel,
        std::vector<std::map<SearchType, unsigned int>>* entrySearchesPerLevel,
        std::vector<unsigned int>* emptyPartitionsPerLevel) const {
        if (level >= nodesPerLevel->size()) {
            nodesPerLevel->push_back(std::map<PartitionType, unsigned int>());
            partitionSearchesPerLevel->push_back(std::map<SearchType, unsigned int>());
            entrySearchesPerLevel->push_back(std::map<SearchType, unsigned int>());
            emptyPartitionsPerLevel->push_back(0);
        }
        (*nodesPerLevel)[level][hp->getPartitionType()] = (*nodesPerLevel)[level][hp->getPartitionType()] + 1;
        (*partitionSearchesPerLevel)[level][hp->get_partition_search_method()->getSearchType()] = (*partitionSearchesPerLevel)[level][hp->get_partition_search_method()->getSearchType()] + 1;
        (*entrySearchesPerLevel)[level][hp->get_entry_search_method()->getSearchType()] = (*entrySearchesPerLevel)[level][hp->get_entry_search_method()->getSearchType()] + 1;
        if(hp->isEmptyPartition())
            (*emptyPartitionsPerLevel)[level] = (*emptyPartitionsPerLevel)[level] + 1;
        for(auto it = hp->get_partitions()->value_begin(), it_end = hp->get_partitions()->value_end(); it != it_end; ++it) {
            stats(*it, level+1, nodesPerLevel, partitionSearchesPerLevel, entrySearchesPerLevel, emptyPartitionsPerLevel);
        }
    }

    /* =================================================================================================================
     * Getter & Setter 
     * ===============================================================================================================*/

    HorizontalPartition<Key, Value>* getPartition() const override { return hp; }

    std::unordered_set<std::pair<Key, Value>, pair_hash> getMissingEntries() const { return missingEntries; }

    bool isValid() const { return hp->isValid(); }

    HorizontalPartition<Key, Value>* getHorizontalPartition() const { return hp; }
    HorizontalPartition<Key, Value>* get_hp_after() const { return hp_after; }

    Workload<Key, Value>* getWorkload() const { return workload; }
    size_t getPartitionCapacity() const { return partition_capacity; }
    size_t getEntryCapacity() const { return entry_capacity; }
    double getPenaltyFactorMissingKeys() const { return penaltyFactorMissingKeys; }
    double getPenaltyFactorEmptyPartitions() const { return penaltyFactorEmptyPartitions; }
    unsigned int getRepetitions() const { return repetitions; }
    std::map<id_type, NodeStatistics<Key, Value>*>* get_stats() const { return id_to_stats; }
};

/**
 * comparator based on the id of the individuals
 */
template<typename Key, typename Value>
struct IndividualIdComparator
{
    bool operator()(const Individual<Key, Value>* lhs, const Individual<Key, Value>* rhs) const { return lhs->getID() < rhs->getID(); }
};

/**
 * comparator based on the fitness of the individuals
 * if the underlying HorizontalPartitions are equal, none of the two is considered smaller
 * otherwise return the comparison of the fitness values if they differ
 * otherwise return the comparison of the IDs
 */
template<typename Key, typename Value>
struct IndividualFitnessComparator
{
    bool operator()(const Individual<Key, Value>* lhs, const Individual<Key, Value>* rhs) const { 
        auto lhsCasted = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(lhs);
        auto rhsCasted = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(rhs);
        if (lhsCasted && rhsCasted && (*(lhsCasted->getHorizontalPartition())) == (*(rhsCasted->getHorizontalPartition()))) {
            return false;
        }
        double lFitness = lhs->getFitness();
        double rFitness = rhs->getFitness();
        if(lFitness != rFitness)
            return lFitness < rFitness;
        return lhs->getID() < rhs->getID();
    }
};

template<typename Key, typename Value>
std::ostream& operator << (std::ostream& stream, const Individual<Key, Value>& ind) {
    ind.print(stream);
    return stream;
}

/**
 * initialize the ID of the individuals
 */
template<typename Key, typename Value>
unsigned int Individual<Key, Value>::ID = 0;

#undef IMPORT
