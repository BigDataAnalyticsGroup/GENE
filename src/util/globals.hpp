#pragma once

#include <cstdint>
#include <map>


/** Singleton class representing the configuration parameters for our generic framework and the genetic algorithm. */
struct Configuration
{
    unsigned int numGenerations;
    unsigned int numMutantsPerGeneration;
    unsigned int populationSize;
    unsigned int tournamentSelectionSize;
    unsigned int initialPopulationSize;
    unsigned int chainMutationLength;
    unsigned int progressBarWidth;
    unsigned int seed;
    unsigned int numSamples;
    unsigned int repetitions;
    unsigned int cores;

    bool maxFitness;
    bool saveDotFiles;
    bool printOutput;
    bool tracing;
    bool structuralMutations;
    bool includePartitions;
    bool includeEntries;
    bool parallel;
    bool sparse;
    bool normal;
    bool printFinalPopulation;
    bool train_test_split;

    size_t partition_capacity_inner;
    size_t partition_capacity_leaf;
    size_t entry_capacity_inner;
    size_t entry_capacity_leaf;

    double quantile;
    double penaltyFactorMissingKeys;
    double penaltyFactorEmptyPartitions;
    double load_factor_entries_leaves;
    double load_factor_partitions_inner;
    double early_abort_fitness;

    std::string resultFile;
    std::string datasetFile;
    std::string workloadFile;
    std::string outputPrefix;

    std::map<PartitionType, double> partitionPriorities;
    std::map<std::string, double> searchPriorities;

    std::vector<std::string> scalingDatasets;
    std::vector<std::string> scalingWorkloads;
    unsigned int scalingStepsize;

    private:
    Configuration() = default;

    public:
    Configuration(Configuration const &) = delete;
    void operator=(Configuration const &) = delete;

    /** Return a reference to the single `Configuration` instance. */
    static Configuration & Get() {
        static Configuration the_configuration;
        return the_configuration;
    }
};
