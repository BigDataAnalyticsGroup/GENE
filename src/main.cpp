#include "cxxopts.hpp"
#include "genetic/generator.hpp"
#include "genetic/mutation.hpp"
#include "genetic/population.hpp"
#include "pool.hpp"
#include "util/globals.hpp"
#include "util/hash.hpp"
#include "util/util.hpp"
#include "util/Dataset.hpp"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <fstream>
#include <future>
#include <limits>
#include <ostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <sys/time.h>

const bool DEBUG = false;

/**
 * print a progress bar containing information about
 * - the current progress measured in generations
 * - the current average fitness (regarding the tournament selection set)
 * - the current best fitness
 * - the approximated memory consumption of the process
 */
void printProgressBar(unsigned int currentGeneration, unsigned int numGenerations, unsigned int barWidth, double quantile, double fitness, double bestFitness, double memory) {
    std::cout << "\r[";
    double progress = double(currentGeneration) / numGenerations;
    unsigned int pos = progress * barWidth;
    for(unsigned int i=0; i<barWidth; i++) {
        if(i<pos) std::cout << "=";
        else if(i==pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100) << "% (" << currentGeneration << " / " << numGenerations << ")  " << (100 * quantile) << "% quantile: " << fitness << " Best: " << bestFitness << " Memory: " << memory << " kB";
    std::cout.flush();
}

/**
 * extract the ratio of insert queries from a given workload file
 */
std::string extract_insert_ratio(const std::string &workload_file) {
    size_type pos_underscore = util::nth_occurrence(workload_file, "_", 5) + 1;
    size_type pos_dot = util::nth_occurrence(workload_file, ".", 2);
    return workload_file.substr(pos_underscore, pos_dot - pos_underscore);
}

/**
 * format a std::chrono::duration in a more human readable form,
 * i.e. print it in days, hours, minutes, seconds and milliseconds
 */
template<typename T>
void printDuration(std::ostream& stream, std::chrono::duration<T> duration) {
    auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    auto const milliseconds = diff % 1000;
    diff /= 1000;
    auto const seconds = diff % 60;
    diff /= 60;
    auto const minutes = diff % 60;
    diff /= 60;
    auto const hours = diff % 24;
    diff /= 24;
    auto const days = diff;

    bool startedPrinting = false;
    if(days >= 1) {
        startedPrinting = true;
        stream << days << (1 != days ? " days" : " day") << " ";
    }
    if(startedPrinting || hours >= 1) {
        startedPrinting = true;
        stream << hours << (1 != hours ? " hours" : " hour") << " ";
    }
    if(startedPrinting || minutes >= 1) {
        startedPrinting = true;
        stream << minutes << (1 != minutes ? " minutes" : " minute") << " ";
    }
    if(startedPrinting || seconds >= 1) {
        startedPrinting = true;
        stream << seconds << (1 != seconds ? " seconds" : " second") << " ";
    }
    stream << milliseconds << (1 != milliseconds ? " milliseconds" : " millisecond") << " ";
}

/**
 * print all settings defined for this execution,
 * useful for logging the execution
 */
void printSettings()
{
    std::cout << "Settings:" << std::endl << std::endl;
    std::cout << "Number of generations: " << Configuration::Get().numGenerations << std::endl;
    std::cout << "Number of mutations per generation: " << Configuration::Get().numMutantsPerGeneration << std::endl;
    std::cout << "Maximum size of population: " << Configuration::Get().populationSize << std::endl;
    std::cout << "Number of individuals considered in tournament selection: " << Configuration::Get().tournamentSelectionSize << std::endl;
    std::cout << "Size of initial population: " << Configuration::Get().initialPopulationSize << std::endl;
    std::cout << "Maximize Fitness: " << Configuration::Get().maxFitness << std::endl;
    std::cout << "Random seed: " << Configuration::Get().seed << std::endl;
    std::cout << "Number of samples to insert in index: " << Configuration::Get().numSamples << std::endl;
    std::cout << "Penalty Factor for samples not inserted: " << Configuration::Get().penaltyFactorMissingKeys << std::endl;
    std::cout << "Penalty Factor for empty partitions: " << Configuration::Get().penaltyFactorEmptyPartitions << std::endl;
    std::cout << "Number of measurement repetitions in fitness function (median taken): " << Configuration::Get().repetitions << std::endl;
    std::cout << "Partition capacity_inner: " << Configuration::Get().partition_capacity_inner << std::endl;
    std::cout << "Partition capacity_leaf: " << Configuration::Get().partition_capacity_leaf << std::endl;
    std::cout << "Entry capacity_inner: " << Configuration::Get().entry_capacity_inner << std::endl;
    std::cout << "Entry capacity_leaf: " << Configuration::Get().entry_capacity_leaf << std::endl;
    std::cout << "Quantile to be reached for population insertion: " << Configuration::Get().quantile << std::endl;
    std::cout << "Result file: " << Configuration::Get().resultFile << std::endl;
    std::cout << "Include partitions in dot files: " << Configuration::Get().includePartitions << std::endl;
    std::cout << "Include entries in dot files: " << Configuration::Get().includeEntries << std::endl;
    std::cout << "Multi-Threaded execution: " << Configuration::Get().parallel << std::endl;
    std::cout << "Number of cores used: " << Configuration::Get().cores << std::endl;
    std::cout << "Use sparse uniform data: " << Configuration::Get().sparse << std::endl;
    std::cout << "Use normal distributed data: " << Configuration::Get().normal << std::endl;
    std::cout << "Print final population: " << Configuration::Get().printFinalPopulation << std::endl;
    std::cout << "File containing dataset: " << Configuration::Get().datasetFile << std::endl;
    std::cout << "File containing workload: " << Configuration::Get().workloadFile << std::endl;
    std::cout << "Load factor of the leaves: " << Configuration::Get().load_factor_entries_leaves << std::endl;
    std::cout << "Load factor of the inner nodes: " << Configuration::Get().load_factor_partitions_inner << std::endl;
    std::cout << "Tracing: " << Configuration::Get().tracing << std::endl;
    std::cout << "Dot file prefix: " << Configuration::Get().outputPrefix << std::endl;
    std::cout << "Structural Mutations: " << Configuration::Get().structuralMutations << std::endl;
    std::cout << "Chain Mutation Length (0 means disabled): " << Configuration::Get().chainMutationLength << std::endl;
    for (auto it = Configuration::Get().partitionPriorities.begin(); it != Configuration::Get().partitionPriorities.end(); it++)
        std::cout << "Priority " << toString(it->first) << ": " << it->second << std::endl;
    for (auto it = Configuration::Get().searchPriorities.begin(); it != Configuration::Get().searchPriorities.end(); it++)
        std::cout << "Priority " << it->first << ": " << it->second << std::endl;
    std::cout << "Fitness to reach for early abort (0 means disabled): " << Configuration::Get().early_abort_fitness << std::endl;
    std::ostringstream oss;
    std::copy(Configuration::Get().scalingDatasets.begin(), Configuration::Get().scalingDatasets.end(), std::ostream_iterator<std::string>(oss, ","));
    std::cout << "Scaling datasets: " << oss.str() << std::endl;
    oss = std::ostringstream();
    std::copy(Configuration::Get().scalingWorkloads.begin(), Configuration::Get().scalingWorkloads.end(), std::ostream_iterator<std::string>(oss, ","));
    std::cout << "Scaling workloads: " << oss.str() << std::endl;
    std::cout << "Scaling stepsize: " << Configuration::Get().scalingStepsize << std::endl;
    std::cout << std::endl;
}

/**
 * main function and entry point for the genetic algorithm
 */
int main(int argc, char* argv[]) {

    // define command line options
    cxxopts::Options options("main_genetic", "genetic algorithm to construct index structures");

    options.add_options()
        ("g,generations", "Number of generations", cxxopts::value<unsigned int>()->default_value("10"))
        ("m,mutants", "Number of mutants per generation", cxxopts::value<unsigned int>()->default_value("10"))
        ("p,population", "Size of the population", cxxopts::value<unsigned int>()->default_value("10"))
        ("t,tournament", "Size of the subset used for tournament selection", cxxopts::value<unsigned int>()->default_value("10"))
        ("i,initial", "Size of the initial population", cxxopts::value<unsigned int>()->default_value("10"))
        ("k,keys", "Number of key-value-pairs to insert", cxxopts::value<unsigned int>()->default_value("10"))
        ("maximize", "Maximize the value of the fitness function", cxxopts::value<bool>()->default_value("false"))
        ("s,seed", "Seed for the random number generator", cxxopts::value<unsigned int>()->default_value("1"))
        ("penaltyMissingKeys", "Penalty used for missing keys", cxxopts::value<double>()->default_value("0")) // superfluous
        ("penaltyEmptyPartitions", "Penalty used for empty partitions", cxxopts::value<double>()->default_value("0")) // superfluous
        ("r,repetitions", "Number of repetitions of the fitness benchmark (median is used as result)", cxxopts::value<unsigned int>()->default_value("3"))
        ("partition_capacity_inner", "Partition capacity of single inner nodes", cxxopts::value<size_t>()->default_value("4"))
        ("partition_capacity_leaf", "Partition capacity of single leaf nodes", cxxopts::value<size_t>()->default_value("4"))
        ("entry_capacity_inner", "Entry capacity of single inner nodes", cxxopts::value<size_t>()->default_value("4"))
        ("entry_capacity_leaf", "Entry capacity of single leaf nodes", cxxopts::value<size_t>()->default_value("4"))
        ("q,quantile", "Quantile to be reached for population insertion", cxxopts::value<double>()->default_value("0.5"))
        ("resultFile", "CSV file containing results of the program run", cxxopts::value<std::string>()->default_value("results.csv"))
        ("includePartitions", "Include partitions in dot files", cxxopts::value<bool>()->default_value("false"))
        ("includeEntries", "Include entries in dot files", cxxopts::value<bool>()->default_value("false"))
        ("parallel", "Multi-threaded execution", cxxopts::value<bool>()->default_value("false"))
        ("c,cores", "Number of cores used for multi-threading", cxxopts::value<unsigned int>()->default_value("4"))
        ("sparse", "Use sparse uniform data", cxxopts::value<bool>()->default_value("false"))
        ("normal", "Use normal distribution for data", cxxopts::value<bool>()->default_value("false"))
        ("printFinalPopulation", "print the final population after finishing all mutations", cxxopts::value<bool>()->default_value("false"))
        ("datasetFile", "file to load dataset from, uses dataset with distribution defined by `sparse` and `normal` parameters and domain [1, maxKeyFactor * keys] as default if not given", cxxopts::value<std::string>()->default_value(""))
        ("workloadFile", "file to load workload from, uses point queries on left half of domain and range queries on right half as default if not given", cxxopts::value<std::string>()->default_value(""))
        ("loadFactorEntriesLeaves", "load factor to use when bulkloading the leaves", cxxopts::value<double>()->default_value("1.0"))
        ("loadFactorPartitionsInner", "load factor use when bulkloading the inner nodes", cxxopts::value<double>()->default_value("1.0"))
        ("prioritySoA", "Priority for SoA when creating new partitions", cxxopts::value<double>()->default_value("1.0"))
        ("priorityTree", "Priority for Tree when creating new partitions", cxxopts::value<double>()->default_value("1.0"))
        ("priorityHash", "Priority for Hash when creating new partitions", cxxopts::value<double>()->default_value("1.0"))
        ("priorityLinear", "Priority for LinearSearch", cxxopts::value<double>()->default_value("0.0"))
        ("priorityBinary", "Priority for BinarySearch", cxxopts::value<double>()->default_value("1.0"))
        ("priorityDefault", "Priority for DefaultSearch", cxxopts::value<double>()->default_value("1.0"))
        ("priorityInterpolation", "Priority for InterpolationSearch", cxxopts::value<double>()->default_value("1.0"))
        ("priorityExponential", "Priority for ExponentialSearch", cxxopts::value<double>()->default_value("1.0"))
        ("priorityLinearRegression", "Priority for LinearRegressionSearch", cxxopts::value<double>()->default_value("0.0"))
        ("noDot", "Do not save best individual of each generation as dot file", cxxopts::value<bool>()->default_value("false"))
        ("minimalOutput", "Print only minimum output to follow progress", cxxopts::value<bool>()->default_value("false"))
        ("tracing", "add annotations to plots", cxxopts::value<bool>()->default_value("false"))
        ("outputPrefix", "prefix for the output files generated during execution", cxxopts::value<std::string>()->default_value(""))
        ("disableStructuralMutations", "disable the use of mutations which affect the structure of the tree", cxxopts::value<bool>()->default_value("false"))
        ("chainMutationLength", "specify the length of chain mutations, zero to disable", cxxopts::value<unsigned int>()->default_value("0"))
        ("earlyAbortFitness", "fitness to reach for early abort, zero to disable", cxxopts::value<double>()->default_value("0.0"))
        ("disableTrainTestSplit", "disable the use of train-test-splits, using all the workload for training", cxxopts::value<bool>()->default_value("false"))
        ("scalingDatasets", "datasets used for upscaling", cxxopts::value<std::vector<std::string>>())
        ("scalingWorkloads", "workloads used for upscaling", cxxopts::value<std::vector<std::string>>())
        ("scalingStepsize", "apply upscaling every scalingStepsize generations", cxxopts::value<unsigned int>()->default_value("0"))
        ("h,help", "Usage")
    ;

    auto result = options.parse(argc, argv);

    // if the help flag is found, simply print the available arguments and their description, then exit
    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    // define settings and store them in the Configuration,
    // mainly based on command line options and their defaults
    // Configuration is a singleton that can easily be shared across all affected classes
    using Key = uint64_t;
    using Value = uint64_t;

    Configuration::Get().numGenerations = result["generations"].as<unsigned int>();
    Configuration::Get().numMutantsPerGeneration = result["mutants"].as<unsigned int>();
    Configuration::Get().populationSize = result["population"].as<unsigned int>();
    Configuration::Get().tournamentSelectionSize = result["tournament"].as<unsigned int>();
    Configuration::Get().initialPopulationSize = result["initial"].as<unsigned int>();
    Configuration::Get().seed = result["seed"].as<unsigned int>();
    Configuration::Get().numSamples = result["keys"].as<unsigned int>();
    Configuration::Get().repetitions = result["repetitions"].as<unsigned int>();
    Configuration::Get().progressBarWidth = 50;
    Configuration::Get().chainMutationLength = result["chainMutationLength"].as<unsigned int>();

    Configuration::Get().maxFitness = result["maximize"].as<bool>();
    Configuration::Get().includePartitions = result["includePartitions"].as<bool>();
    Configuration::Get().includeEntries = result["includeEntries"].as<bool>();

    Configuration::Get().parallel = result["parallel"].as<bool>();
    Configuration::Get().cores = Configuration::Get().parallel ? result["cores"].as<unsigned int>() : 1;

    Configuration::Get().sparse = result["sparse"].as<bool>();
    Configuration::Get().normal = result["normal"].as<bool>();
    Configuration::Get().printFinalPopulation = result["printFinalPopulation"].as<bool>();
    Configuration::Get().saveDotFiles = !(result["noDot"].as<bool>());
    Configuration::Get().printOutput = !(result["minimalOutput"].as<bool>());
    Configuration::Get().tracing = result["tracing"].as<bool>();
    Configuration::Get().structuralMutations = not(result["disableStructuralMutations"].as<bool>());
    Configuration::Get().train_test_split = not(result["disableTrainTestSplit"].as<bool>());

    Configuration::Get().quantile = result["quantile"].as<double>();
    Configuration::Get().penaltyFactorMissingKeys = result["penaltyMissingKeys"].as<double>();
    Configuration::Get().penaltyFactorEmptyPartitions = result["penaltyEmptyPartitions"].as<double>();
    Configuration::Get().load_factor_entries_leaves = result["loadFactorEntriesLeaves"].as<double>();
    Configuration::Get().load_factor_partitions_inner = result["loadFactorPartitionsInner"].as<double>();
    Configuration::Get().early_abort_fitness = result["earlyAbortFitness"].as<double>();

    Configuration::Get().partition_capacity_inner = result["partition_capacity_inner"].as<size_t>();
    Configuration::Get().partition_capacity_leaf = result["partition_capacity_leaf"].as<size_t>();
    Configuration::Get().entry_capacity_inner = result["entry_capacity_inner"].as<size_t>();
    Configuration::Get().entry_capacity_leaf = result["entry_capacity_leaf"].as<size_t>();

    Configuration::Get().resultFile = result["resultFile"].as<std::string>();
    Configuration::Get().datasetFile = result["datasetFile"].as<std::string>();
    Configuration::Get().workloadFile = result["workloadFile"].as<std::string>();
    Configuration::Get().outputPrefix = result["outputPrefix"].as<std::string>();

    Configuration::Get().partitionPriorities[PartitionType::SoAPartition] = result["prioritySoA"].as<double>();
    Configuration::Get().partitionPriorities[PartitionType::TreePartition] = result["priorityTree"].as<double>();
    Configuration::Get().partitionPriorities[PartitionType::HashPartition] = result["priorityHash"].as<double>();

    Configuration::Get().searchPriorities["LinearSearch"] = result["priorityLinear"].as<double>();
    Configuration::Get().searchPriorities["BinarySearch"] = result["priorityBinary"].as<double>();
    Configuration::Get().searchPriorities["DefaultSearch"] = result["priorityDefault"].as<double>();
    Configuration::Get().searchPriorities["InterpolationSearch"] = result["priorityInterpolation"].as<double>();
    Configuration::Get().searchPriorities["ExponentialSearch"] = result["priorityExponential"].as<double>();
    Configuration::Get().searchPriorities["LinearRegressionSearch"] = result["priorityLinearRegression"].as<double>();

    if (result.count("scalingDatasets")) {
        Configuration::Get().scalingDatasets = result["scalingDatasets"].as<std::vector<std::string>>();
    } else {
        Configuration::Get().scalingDatasets = std::vector<std::string>();
    }
    if (result.count("scalingWorkloads")) {
        Configuration::Get().scalingWorkloads = result["scalingWorkloads"].as<std::vector<std::string>>();
    } else {
        Configuration::Get().scalingWorkloads = std::vector<std::string>();
    }
    Configuration::Get().scalingStepsize = result["scalingStepsize"].as<unsigned int>();
    if(Configuration::Get().scalingWorkloads.size() > 0 and Configuration::Get().scalingWorkloads.size() != Configuration::Get().scalingDatasets.size())
        throw std::invalid_argument("scalingWorkloads must either be empty or have the same length as scalingDatasets");

    Key lowerBound = Key(1);
    Key upperBound = std::numeric_limits<Key>::max()-1;

    // a distribution can either be sparse or normal, not both at the same time
    if (Configuration::Get().sparse && Configuration::Get().normal)
    throw std::invalid_argument("Sparse and normal option can not be combined");

    std::cout << Configuration::Get().scalingDatasets.size() << std::endl;

    // print the settings for this execution
    printSettings();

    // define workload and dataset,
    // prefer files given via command line
    // only use default values if no file is given
    Dataset<Key, Value>* dataset;
    Workload<Key, Value>* workload;

    //Distribution distribution = Distribution::uniform_dense;
    //if (Configuration::Get().normal)
    //    distribution = Distribution::normal;
    if (not Configuration::Get().datasetFile.empty()) {
        dataset = Dataset<Key, Value>::fromFile(Configuration::Get().datasetFile);
    } else {
        if (Configuration::Get().sparse) {
            dataset = Dataset<Key, Value>::Create_Uniform_Sparse(Configuration::Get().numSamples, {lowerBound, upperBound}, Configuration::Get().seed);
        } else {
            dataset = Dataset<Key, Value>::Create_Uniform_Dense(Configuration::Get().numSamples, lowerBound);
        }
    }
    if (not Configuration::Get().workloadFile.empty()) {
        workload = Workload<Key, Value>::from_file(dataset, Configuration::Get().workloadFile);
    } else {
        workload = new Workload<Key, Value>(dataset, Configuration::Get().seed);
        workload->add_point_uniform_index_based(dataset->size(), {0, dataset->size()-1});
    }

    Workload<Key, Value>* train;
    Workload<Key, Value>* test;
    if (Configuration::Get().train_test_split) {
        /* Split workload in train/test. */
        auto split = workload->split_train_test(1, true);
        train = split.first;
        test = split.second;
    } else {
        train = workload;
        test = workload;
    }

    // load datasets and workloads used for upscaling experiments
    std::vector<Dataset<Key, Value>*> scalingDatasets;
    std::vector<Workload<Key, Value>*> scalingWorkloads;
    for (unsigned long i=0; i<Configuration::Get().scalingDatasets.size(); i++) {
        auto ds = Dataset<Key, Value>::fromFile(Configuration::Get().scalingDatasets[i]);
        scalingDatasets.push_back(ds);
        if (Configuration::Get().scalingWorkloads.size() > 0)
            scalingWorkloads.push_back(Workload<Key, Value>::from_file(ds, Configuration::Get().scalingWorkloads[i]));
        else
            scalingWorkloads.push_back(workload);
    }

    // define a Generator which initializes the population
    Generator<Key, Value>* gen = new FullHorizontalPartitionGenerator<Key, Value>(dataset, train);

    // initialize the thread pool and the csv files used to store the results
    thread_pool tp(Configuration::Get().cores);
    std::ofstream csv;
    csv.open(Configuration::Get().resultFile);
    csv << "Generation,Id,Fitness,MedianFitness,Runtime" << std::endl;
    std::ofstream csvPerLevel;
    std::string resultFilePerLevel = Configuration::Get().resultFile;
    std::string::size_type fileExtensionPos = resultFilePerLevel.rfind('.', resultFilePerLevel.length());
    std::string newExtension = "PerLevel.csv";
    if (fileExtensionPos != std::string::npos)
        resultFilePerLevel.replace(fileExtensionPos, newExtension.length(), newExtension);
    csvPerLevel.open(resultFilePerLevel);
    csvPerLevel << "ID" << "," << "Level" << "," << "ContentType" << "," << "Content" << std::endl;
    std::ofstream csvScaled;
    std::string resultFileScaled = Configuration::Get().resultFile;
    fileExtensionPos = resultFileScaled.rfind('.', resultFileScaled.length());
    newExtension = "_Scaled.csv";
    if (fileExtensionPos != std::string::npos)
        resultFileScaled.replace(fileExtensionPos, newExtension.length(), newExtension);
    csvScaled.open(resultFileScaled);
    csvScaled << "Generation,Dataset,Size,Id,Fitness" << std::endl;

    // initialize the population and retrieve the best individual
    auto startPop = std::chrono::system_clock::now();
    std::time_t tPop = std::chrono::system_clock::to_time_t(startPop);
    std::cout << "Starting population creation at: " << std::ctime(&tPop) << std::endl;
    //Population<Key, Value>* pop = new Population<Key, Value>(gen, populationSize, initialPopulationSize, maxFitness);
    Population<Key, Value>* pop = new Population<Key, Value>(gen);
    std::pair<const Individual<Key, Value>*, double> p = pop->tournamentSelection(pop->size(), 0.5);
    const HorizontalPartitionIndividual<Key, Value>* bestIndividual = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(p.first);
    double averageFitness = p.second;
    double stddev_fitness = pop->stddev_fitness();
    double min_fitness = pop->best_individual()->getFitness();
    double max_fitness = pop->worst_individual()->getFitness();

    // save initial results
    if (Configuration::Get().saveDotFiles) {
        if (Configuration::Get().tracing) {
            /* Annotations. */
            std::stringstream ss;
            ss << "Generation: 0\n";
            ss << "Population:"
                    << " size=" << pop->size()
                    << " | average_fitness=" << averageFitness
                    << " | stddev_fitness=" << stddev_fitness
                    << " | min_fitness=" << min_fitness
                    << " | max_fitness=" << max_fitness
                    << "\n";
            ss << "Best Individual:"
                    << " ID=" << bestIndividual->getID()
                    << " | fitness=" << bestIndividual->getFitness()
                    << "\n";
            bestIndividual->toDotFile(Configuration::Get().outputPrefix + "bestIndividualGeneration" + std::to_string(0) + ".dot", Configuration::Get().includePartitions, Configuration::Get().includeEntries, ss.str());
        } else
            bestIndividual->toDotFile(Configuration::Get().outputPrefix + "bestIndividualGeneration" + std::to_string(0) + ".dot", Configuration::Get().includePartitions, Configuration::Get().includeEntries);
    }

    if (Configuration::Get().printOutput) {
        std::cout << "Initial population created" << std::endl;
        std::cout << "Best individual of the initial population: " << std::endl << std::endl;
        std::cout << bestIndividual->describe(csvPerLevel) << std::endl << std::endl;
    } else {
        bestIndividual->describe(csvPerLevel);
    std::cout << "Best initial fitness: " << bestIndividual->getFitness() << std::endl;
    }
    csv << "0" << "," << bestIndividual->getID() << "," << bestIndividual->getFitness() << "," << averageFitness << "," << 0 << std::endl;

    // compute results for scaling experiments
    unsigned int bestIndividualScaledId = bestIndividual->getID();
    auto bestIndividualScaledFitness = std::vector<double>(Configuration::Get().scalingDatasets.size());
    if (Configuration::Get().scalingStepsize > 0) {
        for (unsigned int d=0; d<Configuration::Get().scalingDatasets.size(); d++) {
            auto scaled_hp = scale_existing_index<Key, Value>(bestIndividual->getHorizontalPartition(), scalingDatasets[d]->begin(), scalingDatasets[d]->end(), false);
            auto scaled_ind = new HorizontalPartitionIndividual<Key, Value>(scaled_hp, scalingWorkloads[d], bestIndividual->getPartitionCapacity(), bestIndividual->getEntryCapacity(),
                                         bestIndividual->getPenaltyFactorMissingKeys(), bestIndividual->getPenaltyFactorEmptyPartitions(), bestIndividual->getRepetitions(), 
                                         bestIndividual->getID(), bestIndividual->getMissingEntries());
            scaled_ind->toDotFile(Configuration::Get().outputPrefix + "bestIndividualGeneration" + std::to_string(0) + "_Scaled_" + std::to_string(d) + ".dot", Configuration::Get().includePartitions, Configuration::Get().includeEntries);
            bestIndividualScaledFitness[d] = scaled_ind->getFitness();
            delete(scaled_ind);
            csvScaled << 0 << "," << d << "," << scalingDatasets[d]->size() << "," << bestIndividualScaledId << "," << bestIndividualScaledFitness[d] << std::endl;
        }
    }

    // initialize the set of mutations to use
    std::set<Mutation<Key, Value>*> mutations;

    /* Basic mutations */
    mutations.insert(new PartitionTypeChangeMutation<Key, Value>(2.0, lowerBound, upperBound));
    mutations.insert(new PartitionSearchChangeMutation<Key, Value>(1.0, lowerBound, upperBound));
    mutations.insert(new EntrySearchChangeMutation<Key, Value>(1.0, lowerBound, upperBound));
    if (Configuration::Get().structuralMutations) {
        mutations.insert(new HorizontalMergeMutation<Key, Value>(4.0, upperBound));
        mutations.insert(new VerticalMergeMutation<Key, Value>(4.0, upperBound));
        mutations.insert(new HorizontalSplitMutation<Key, Value>(0.25, lowerBound, upperBound));
        mutations.insert(new VerticalSplitMutation<Key, Value>(0.25, lowerBound, upperBound));
    }

    if (Configuration::Get().chainMutationLength > 0) {
        /* Chain mutation */
        auto chainMutations = mutations;
        mutations.insert(new LocalChainMutation<Key, Value>(1.0, lowerBound, upperBound, chainMutations, Configuration::Get().maxFitness));
    }

    // define probability distribution for the mutations
    std::vector<double> probabilities;
    for(auto mut : mutations)
        probabilities.push_back(mut -> getProbability());
    std::mt19937 randomGen(std::random_device{}());
    std::discrete_distribution<int> dist(probabilities.begin(), probabilities.end());

    // initialize variables for best individual and memory usage,
    // print initial progress bar
    unsigned int bestIndividualId = bestIndividual->getID();
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    auto start = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(start);
    std::cout << "Starting mutations at: " << std::ctime(&t) << std::endl;
    printProgressBar(0, Configuration::Get().numGenerations, Configuration::Get().progressBarWidth, Configuration::Get().quantile, averageFitness, bestIndividual->getFitness(), r_usage.ru_maxrss);

    // initialize stringstream for potential debug annotations
    std::stringstream ss;
    // loop until numGenerations have been executed
    for(unsigned int g=0; g<Configuration::Get().numGenerations; g++) {
        // perform a tournament selection to determine the best individual to mutate
        std::stringstream ts_ss;
        p = pop->tournamentSelection(Configuration::Get().tournamentSelectionSize, Configuration::Get().quantile, ts_ss);
        auto bestIndividual_before = pop->best_individual();
        const Individual<Key, Value>* sample = p.first;
        averageFitness = p.second;
        
        // store further information in string stream if more detailed tracing is desired
        if (Configuration::Get().tracing) {
            /* Annotations. */
            ss << "Generation: " <<  g+1 << "\n";
            ss << "Population:"
                    << " size=" << pop->size()
                    << " | median_fitness=" << averageFitness
                    << " | stddev_fitness=" << stddev_fitness
                    << " | min_fitness=" << min_fitness
                    << " | max_fitness=" << max_fitness
                    << "\n";
            ss << "Sample (TS): " << ts_ss.str() << "\n";
            ss << "Best Individual global before:"
                    << " ID=" << bestIndividual_before->getID()
                    << " | fitness=" << bestIndividual_before->getFitness()
                    << " | stddev_fitness=" << bestIndividual_before->get_stddev_fitness()
                    << "\n";
            ss << "Best Individual sample:"
                    << " ID=" << sample->getID()
                    << " | fitness=" << sample->getFitness()
                    << " | stddev_fitness=" << sample->get_stddev_fitness()
                    << "\n";
        }

        /* Multi-threaded case */
        if (Configuration::Get().parallel) {
            // define futures, i.e. tasks for the thread pool
            std::vector<std::future<const Individual<Key, Value>*>> futures(Configuration::Get().numMutantsPerGeneration);
            std::vector<Mutation<Key, Value>*> applied_mutations;
            for(unsigned int m=0; m<Configuration::Get().numMutantsPerGeneration; m++) {
                int randomChoice = dist(randomGen);
                auto it = mutations.begin();
                std::advance(it, randomChoice);
                applied_mutations.push_back(*it);
                futures[m] = tp.enqueue_task(&Mutation<Key, Value>::mutate, *it, sample);
            }
            // loop over futures to retrieve the results
            // futures' get function is blocking until result is available
            std::vector<const Individual<Key, Value>*> indToAdd, indToDelete;
            ss << "Mutants(ID, mutation, fitness, stddev): ";
            for(unsigned int m=0; m<Configuration::Get().numMutantsPerGeneration; m++) {
                auto mutant = futures[m].get();
                if (Configuration::Get().tracing) {
                    if (m != 0) ss << " | ";
                    ss << "("
                        << mutant->getID() << ","
                        << applied_mutations[m]->getName() << ","
                        << mutant->getFitness() << ","
                        << mutant->get_stddev_fitness()
                        << ")";
                }
                if (DEBUG)
                    std::cout << mutant->describe() << std::endl;
                if (mutant == sample) {
                    continue;
                }
                /* Check if mutant has to be added to population or deleted. */
                auto mutant_fitness = mutant->getFitness();
                if ((Configuration::Get().maxFitness ? mutant_fitness >= averageFitness : mutant_fitness <= averageFitness) || Configuration::Get().quantile == 0.0) {
                    indToAdd.push_back(mutant);
                } else {
                    indToDelete.push_back(mutant);
                }
            }
            if (Configuration::Get().tracing)
                ss << "\n";

            // add individuals to population
            for (auto it = indToAdd.begin(); it != indToAdd.end(); ) {
                pop->add(*it);
                it = indToAdd.erase(it);
            }
            assert(indToAdd.empty());
            // delete individuals which are no longer necessary
            for (auto it = indToDelete.begin(); it != indToDelete.end(); ) {
                delete(*it);
                it = indToDelete.erase(it);
            }
            assert(indToDelete.empty());
        } 
        /* Single-threaded case */
        else {
            std::vector<const Individual<Key, Value>*> mutants(Configuration::Get().numMutantsPerGeneration);
            std::vector<Mutation<Key, Value>*> applied_mutations;
            for(unsigned int m=0; m<Configuration::Get().numMutantsPerGeneration; m++) {
                int randomChoice = dist(randomGen);
                auto it = mutations.begin();
                std::advance(it, randomChoice);
                applied_mutations.push_back(*it);
                mutants[m] = (*it)->mutate(sample);
            }
            // loop over mutants retrieved in previous step
            std::vector<const Individual<Key, Value>*> indToAdd, indToDelete;
            ss << "Mutants(ID, mutation, fitness, stddev): ";
            for (unsigned int  m=0; m < Configuration::Get().numMutantsPerGeneration; m++) {
                auto mutant = mutants[m];
                if (Configuration::Get().tracing) {
                    if (m > 0) ss << " | ";
                    ss << "("
                        << mutant->getID() << ","
                        << applied_mutations[m]->getName() << ","
                        << mutant->getFitness() << ","
                        << mutant->get_stddev_fitness()
                        << ")";
                }
                if (DEBUG)
                    std::cout << mutant->describe() << std::endl;

                // ignore failed mutations returning previous best individual
                if (mutant == sample) {
                    continue;
                }
                // Check if mutant has to be added to population
                if ((Configuration::Get().maxFitness  ? mutant->getFitness() >= averageFitness
                                        : mutant->getFitness() <= averageFitness) || Configuration::Get().quantile == 0.0) {
                    indToAdd.push_back(mutant);
                } 
                // delete individuals no longer necessary
                else {
                    indToDelete.push_back(mutant);
                }
            }
            if (Configuration::Get().tracing)
                ss << "\n";
            // add individuals to population
            for (auto it = indToAdd.begin(); it != indToAdd.end(); ) {
                pop->add(*it);
                it = indToAdd.erase(it);
            }
            assert(indToAdd.empty());
            // delete individuals which are no longer necessary
            for (auto it = indToDelete.begin(); it != indToDelete.end(); ) {
                delete(*it);
                it = indToDelete.erase(it);
            }
            assert(indToDelete.empty());
        }

        // retrieve best individual by applying a tournament selection over the whole population
        p = pop -> tournamentSelection(pop -> size(), 0.5);
        const HorizontalPartitionIndividual<Key, Value>* bestIndividual = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(p.first);
        ss << "Best Individual after (shown):"
                << " ID=" << bestIndividual->getID()
                << " | fitness=" << bestIndividual->getFitness()
                << " | stddev_fitness=" << bestIndividual->get_stddev_fitness()
                << "\n";

        /* store best Individual in graphviz dot file if desired. */
        if (Configuration::Get().tracing) {
            // With `tracing` flag set, dot best individual every generation.
            bestIndividual->toDotFile(Configuration::Get().outputPrefix + "bestIndividualGeneration" + std::to_string(g+1) + ".dot", Configuration::Get().includePartitions, Configuration::Get().includeEntries, ss.str());
            bestIndividualId = bestIndividual->getID();
            if (Configuration::Get().printOutput) {
                std::cout << std::endl << "Best individual of generation: " << g+1 << std::endl << std::endl;
                std::cout << bestIndividual->describe(csvPerLevel) << std::endl << std::endl;
            } else
                bestIndividual->describe(csvPerLevel);
        } else if (bestIndividualId != bestIndividual->getID()) {
            // Otherwise only dot best individual if there is a better one
            if (Configuration::Get().saveDotFiles)
                bestIndividual->toDotFile(Configuration::Get().outputPrefix + "bestIndividualGeneration" + std::to_string(g+1) + ".dot", Configuration::Get().includePartitions, Configuration::Get().includeEntries);
            bestIndividualId = bestIndividual->getID();
            if (Configuration::Get().printOutput) {
                std::cout << std::endl << "Best individual of generation: " << g+1 << std::endl << std::endl;
                std::cout << bestIndividual->describe(csvPerLevel) << std::endl << std::endl;
            } else
                bestIndividual->describe(csvPerLevel);
        }

        // print progress bar
        getrusage(RUSAGE_SELF, &r_usage);
        printProgressBar(g+1, Configuration::Get().numGenerations, Configuration::Get().progressBarWidth, Configuration::Get().quantile, averageFitness, bestIndividual->getFitness(), r_usage.ru_maxrss);

        // save results of this iteration
        std::chrono::duration<double> duration = std::chrono::system_clock::now() - start;
        auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
        csv << g+1 << "," << bestIndividual->getID() << "," << bestIndividual->getFitness() << "," << p.second << "," << diff << std::endl;

        // evaluate scaling experiments
        if (Configuration::Get().scalingStepsize > 0 and ((g+1) % Configuration::Get().scalingStepsize) == 0) {
            for (unsigned int d=0; d<Configuration::Get().scalingDatasets.size(); d++) {
                if(bestIndividualScaledId != bestIndividual->getID()) {
                    auto scaled_hp = scale_existing_index<Key, Value>(bestIndividual->getHorizontalPartition(), scalingDatasets[d]->begin(), scalingDatasets[d]->end(), false);
                    auto scaled_ind = new HorizontalPartitionIndividual<Key, Value>(scaled_hp, scalingWorkloads[d], bestIndividual->getPartitionCapacity(), bestIndividual->getEntryCapacity(),
                                                 bestIndividual->getPenaltyFactorMissingKeys(), bestIndividual->getPenaltyFactorEmptyPartitions(), bestIndividual->getRepetitions(), 
                                                 bestIndividual->getID(), bestIndividual->getMissingEntries());
                    scaled_ind->toDotFile(Configuration::Get().outputPrefix + "bestIndividualGeneration" + std::to_string(g+1) + "_Scaled_" + std::to_string(d) + ".dot", Configuration::Get().includePartitions, Configuration::Get().includeEntries);
                    bestIndividualScaledFitness[d] = scaled_ind->getFitness();
                    delete(scaled_ind);
                }
                csvScaled << g+1 << "," << d << "," << scalingDatasets[d]->size() << "," << bestIndividualScaledId << "," << bestIndividualScaledFitness[d] << std::endl;
            }
            bestIndividualScaledId = bestIndividual->getID();
        }
    
        // check for early abort
        if (Configuration::Get().early_abort_fitness > 0) {
            if (Configuration::Get().maxFitness and bestIndividual->getFitness() > Configuration::Get().early_abort_fitness)
                break;
            else if (not Configuration::Get().maxFitness and bestIndividual->getFitness() < Configuration::Get().early_abort_fitness)
                break;
        }
    }

    // print total runtime
    std::cout << std::endl << std::endl;
    std::chrono::duration<double> duration = std::chrono::system_clock::now() - start;
    std::cout << "Finished in: ";
    printDuration(std::cout, duration);
    std::cout << std::endl << std::endl;

    // retrieve best individual of final population
    p = pop->tournamentSelection(pop -> size(), Configuration::Get().quantile);
    const HorizontalPartitionIndividual<Key, Value>* bestIndividualFinalPopulation = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(p.first);
    averageFitness = p.second;
    std::cout << "Best fitness value stored: " << bestIndividualFinalPopulation->getFitness() << std::endl;
    if (Configuration::Get().printOutput) {
        std::cout << "Best individual of the final population:" << std::endl << std::endl;
        std::cout << bestIndividualFinalPopulation->describe() << std::endl << std::endl;
        if (Configuration::Get().printFinalPopulation)
            std::cout << "Final population: " << *pop << std::endl;

        std::cout << "Mutation execution counts:" << std::endl;
        for(Mutation<Key, Value>* mut : mutations) {
            std::cout << mut->getName() << " : " << mut->getExecutionCounter() << std::endl;
        }
        std::cout << std::endl;
    }
    std::stringstream annotation_train;
    annotation_train << "Fitness: " << bestIndividualFinalPopulation->getFitness() << std::endl;

    std::stringstream output_file_train;
    annotation_train << "Prefix: " << Configuration::Get().outputPrefix;
    output_file_train << Configuration::Get().outputPrefix << "bestTrain" << ".dot";
    bestIndividualFinalPopulation->toDotFile(output_file_train.str(), Configuration::Get().includePartitions, Configuration::Get().includeEntries, annotation_train.str());
    std::stringstream output_file_train_after;
    output_file_train_after << Configuration::Get().outputPrefix << "bestTrainAfter" << ".dot";
    bestIndividualFinalPopulation->toDotFile(output_file_train_after.str(), Configuration::Get().includePartitions, Configuration::Get().includeEntries, annotation_train.str(), true);

    /* Evaluate on test data. */
    std::pair<const HorizontalPartitionIndividual<Key, Value>*, double> best;
    if (Configuration::Get().maxFitness) { best.second = std::numeric_limits<double>::min(); }
    else { best.second = std::numeric_limits<double>::max(); }

    for (auto ind : *pop) {
        auto h_ind = static_cast<const HorizontalPartitionIndividual<Key, Value>*>(ind);
        auto fitness = (h_ind->calculateFitness(test)).get_fitness();
        if (Configuration::Get().maxFitness) {
            if (fitness > best.second) {
                best.first = h_ind;
                best.second = fitness;
            }
        } else {
            if (fitness < best.second) {
                best.first = h_ind;
                best.second = fitness;
            }
        }
    }
    if (not best.first) throw std::runtime_error("Best Individual evaluated on test data is null.");

    std::stringstream annotation_test;
    annotation_test << "Fitness: " << best.second << std::endl;

    std::stringstream output_file_test;
    output_file_test << Configuration::Get().outputPrefix << "bestTest" << ".dot";
    annotation_test << "Prefix: " << Configuration::Get().outputPrefix;
    best.first->toDotFile(output_file_test.str(), Configuration::Get().includePartitions, Configuration::Get().includeEntries, annotation_test.str());


    // clean up
    csv.close();
    csvPerLevel.close();
    csvScaled.close();
    delete(gen);
    delete(pop);
    while(!mutations.empty()) {
        auto it = mutations.begin();
        mutations.erase(it);
        delete(*it);
    }
    for (auto d: scalingDatasets) {
        if (d != dataset)
            delete(d);
    }
    delete(dataset);

    for (auto w: scalingWorkloads) {
        if (w != workload)
            delete(w);
    }
    delete(workload);
    if (Configuration::Get().train_test_split) {
        delete(train);
        delete(test);
    }

    return 0;
}
