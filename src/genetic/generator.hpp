#pragma once

#include "util/Dataset.hpp"
#include "util/Workload.hpp"
#include "genetic/individual.hpp"
#include "hp/Index.hpp"
#include "util/globals.hpp"
#include <set>
#include <ostream>
#include <sstream>
#include <fstream>

const bool DEBUG_GENERATOR = false;

/**
 * abstract base class for the Generator
 * uses its generateIndivuals function to generate initial individuals for the population
 */
template<typename Key, typename Value>
class Generator {
    public:
        virtual ~Generator() {};
        /**
         * abstract method generating a set of `num` individuals.
         */
        virtual std::set<const Individual<Key, Value>*, IndividualFitnessComparator<Key, Value>> generateIndividuals(unsigned int num) = 0;
};

/**
 * print a progress bar for the initialization process
 */
void printProgress(double progress, int barWidth=100) {
    std::cout << "[";
    int pos = barWidth * progress;
    for (int j = 0; j < barWidth; ++j) {
        if (j < pos) std::cout << "=";
        else if (j == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

/**
 * concrete Generator implementation
 * uses B-tree like bulkloading to generate individuals
 */
template<typename Key, typename Value>
class FullHorizontalPartitionGenerator : public Generator<Key, Value> {
    public:
        ~FullHorizontalPartitionGenerator() {}
        
        FullHorizontalPartitionGenerator(Dataset<Key, Value>* dataset, Workload<Key, Value>* workload)
            : seed(Configuration::Get().seed), // seed for random number generator
              numSamples(Configuration::Get().numSamples), // number of keys
              dataset(dataset), // dataset to bulkload from
              workload(workload), // workload to evaluate the fitness on
              partition_capacity_inner(Configuration::Get().partition_capacity_inner),  // capacity for child partitions in inner nodes
              partition_capacity_leaf(Configuration::Get().partition_capacity_leaf), // capacity for child partitions in leaf nodes
              entry_capacity_inner(Configuration::Get().entry_capacity_inner), // capacity for entries in inner nodes
              entry_capacity_leaf(Configuration::Get().entry_capacity_leaf), // capacity for entries in leaf nodes
              penaltyFactorMissingKeys(Configuration::Get().penaltyFactorMissingKeys), // fitness penalty for missing keys
              penaltyFactorEmptyPartitions(Configuration::Get().penaltyFactorEmptyPartitions), // fitness penalty for empty partitions
              repetitions(Configuration::Get().repetitions), // measurement repetitions during fitness evaluation
              loadFactorEntriesLeaves(Configuration::Get().load_factor_entries_leaves), // load factor for the entries contained in the leaves
              loadFactorPartitionsInner(Configuration::Get().load_factor_partitions_inner) // load factor the child partitions contained in the inner nodes
        {}

        /**
         * generate a set containing `num` individuals
         * uses a B-tree like bulkloading algorithm based on random data layouts
         */
        std::set<const Individual<Key, Value>*, IndividualFitnessComparator<Key, Value>> generateIndividuals(unsigned int num) override {
            std::set<const Individual<Key, Value>*, IndividualFitnessComparator<Key, Value>> individuals;
            printProgress(0);
            // loop until all individuals have been created
            for(unsigned int i=0; i<num; i++) {
                // bulkload a HorizontalPartition
                HorizontalPartition<Key, Value>* hp;
                if (Configuration::Get().random_layouts) {
                    // as B-tree like structure with random data layout and random search methods
                    hp = bulkload_random_tree<Key, Value>(dataset->begin(), dataset->end(), partition_capacity_leaf, entry_capacity_leaf, partition_capacity_inner, entry_capacity_inner, loadFactorEntriesLeaves, loadFactorPartitionsInner);
                } else {
                    // as B-tree like structure with SoA data layout and binary search
                    hp = bulkload_soa_tree<Key, Value>(dataset->begin(), dataset->end(), partition_capacity_leaf, entry_capacity_leaf, partition_capacity_inner, entry_capacity_inner, loadFactorEntriesLeaves, loadFactorPartitionsInner);
                }

                // initialize an individual and measure its fitness
                HorizontalPartitionIndividual<Key, Value>* ind = new HorizontalPartitionIndividual<Key, Value>(hp,
                                                                        this->workload,
                                                                        this->partition_capacity_inner, this->entry_capacity_leaf,
                                                                        this->penaltyFactorMissingKeys, this->penaltyFactorEmptyPartitions,
                                                                        this->repetitions, 0);
                // save individual in result set
                if(not individuals.insert(ind).second) {
                    if (DEBUG_GENERATOR)
                        std::cout << "Generator deleting individual in generateIndividuals: " << ind->getID() << std::endl;
                    delete(ind);
                }
                printProgress(double(i+1) / num);
                }
            std::cout << std::endl;
            // return generated individuals
            return individuals;
        }

    private:
        unsigned int seed, numSamples;
        Dataset<Key, Value>* dataset;
        Workload<Key, Value>* workload;
        size_t partition_capacity_inner;
        size_t partition_capacity_leaf;
        size_t entry_capacity_inner;
        size_t entry_capacity_leaf;
        double penaltyFactorMissingKeys, penaltyFactorEmptyPartitions;
        unsigned int repetitions;
        double loadFactorEntriesLeaves, loadFactorPartitionsInner;
};
