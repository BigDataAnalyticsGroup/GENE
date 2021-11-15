#include "genetic/population.hpp"
#include "genetic/generator.hpp"
#include "genetic/mutation.hpp"
#include "util/globals.hpp"
#include "util/hash.hpp"
#include "cxxopts.hpp"
#include <ostream>
#include <set>
#include <chrono>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <future>
#include <sys/time.h>
#include <sys/resource.h>
#include <algorithm>
#include <limits>
#include "pool.hpp"

const bool DEBUG = false;

void printProgressBar(unsigned int currentGeneration, unsigned int numGenerations, unsigned int barWidth, double quantile, double fitness, double memory) {
    std::cout << "\r[";
    double progress = double(currentGeneration) / numGenerations;
    unsigned int pos = progress * barWidth;
    for(unsigned int i=0; i<barWidth; i++) {
        if(i<pos) std::cout << "=";
        else if(i==pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100) << "% (" << currentGeneration << " / " << numGenerations << ")  " << (100 * quantile) << "% quantile: " << fitness << " Memory: " << memory << " kB";
    std::cout.flush();
}

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
    std::cout << "Partition capacity: " << Configuration::Get().partition_capacity_inner << std::endl;
    std::cout << "Entry capacity: " << Configuration::Get().entry_capacity_leaf << std::endl;
    std::cout << "Quantile to be reached for population insertion: " << Configuration::Get().quantile << std::endl;
    std::cout << "Result file: " << Configuration::Get().resultFile << std::endl;
    std::cout << "Include partitions in dot files: " << Configuration::Get().includePartitions << std::endl;
    std::cout << "Include entries in dot files: " << Configuration::Get().includeEntries << std::endl;
    std::cout << "Multi-Threaded execution: " << Configuration::Get().parallel << std::endl;
    std::cout << "Number of cores used: " << Configuration::Get().cores << std::endl;
    std::cout << "Use sparse uniform data: " << Configuration::Get().sparse << std::endl;
    std::cout << "Print final population: " << Configuration::Get().printFinalPopulation << std::endl;
    std::cout << "File containing dataset: " << Configuration::Get().datasetFile << std::endl;
    std::cout << "File containing workload: " << Configuration::Get().workloadFile << std::endl;
    std::cout << "Load factor of the leaves: " << Configuration::Get().load_factor_entries_leaves << std::endl;
    std::cout << "Load factor of the inner nodes: " << Configuration::Get().load_factor_partitions_inner << std::endl;
    for (auto it = Configuration::Get().partitionPriorities.begin(); it != Configuration::Get().partitionPriorities.end(); it++)
        std::cout << "Priority " << toString(it->first) << ": " << it->second << std::endl;
    for (auto it = Configuration::Get().searchPriorities.begin(); it != Configuration::Get().searchPriorities.end(); it++)
        std::cout << "Priority " << it->first << ": " << it->second << std::endl;
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {

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
        ("penaltyMissingKeys", "Penalty used for missing keys", cxxopts::value<double>()->default_value("1000"))
        ("penaltyEmptyPartitions", "Penalty used for empty partitions", cxxopts::value<double>()->default_value("10"))
        ("r,repetitions", "Number of repetitions of the fitness benchmark (median is used as result)", cxxopts::value<unsigned int>()->default_value("3"))
        ("partition_capacity", "Partition capacity of single nodes", cxxopts::value<size_t>()->default_value("4"))
        ("entry_capacity", "Entry capacity of single nodes", cxxopts::value<size_t>()->default_value("4"))
        ("q,quantile", "Quantile to be reached for population insertion", cxxopts::value<double>()->default_value("0.5"))
        ("resultFile", "CSV file containing results of the program run", cxxopts::value<std::string>()->default_value("results.csv"))
        ("includePartitions", "Include partitions in dot files", cxxopts::value<bool>()->default_value("false"))
        ("includeEntries", "Include entries in dot files", cxxopts::value<bool>()->default_value("false"))
	    ("parallel", "Multi-threaded execution", cxxopts::value<bool>()->default_value("false"))
	    ("c,cores", "Number of cores used for multi-threading", cxxopts::value<unsigned int>()->default_value("4"))
	    ("sparse", "Use sparse uniform data", cxxopts::value<bool>()->default_value("false"))
        ("printFinalPopulation", "print the final population after finishing all mutations", cxxopts::value<bool>()->default_value("false"))
	    ("datasetFile", "file to load dataset from, uses dataset with distribution defined by `sparse` parameter and the largest domain possible given the datatype as default", cxxopts::value<std::string>()->default_value(""))
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
        ("h,help", "Usage")
    ;

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    using Key = uint32_t;
    using Value = uint32_t;

    Configuration::Get().numGenerations = result["generations"].as<unsigned int>();
    Configuration::Get().numMutantsPerGeneration = result["mutants"].as<unsigned int>();
    Configuration::Get().populationSize = result["population"].as<unsigned int>();
    Configuration::Get().tournamentSelectionSize = result["tournament"].as<unsigned int>();
    Configuration::Get().initialPopulationSize = result["initial"].as<unsigned int>();
    Configuration::Get().maxFitness = result["maximize"].as<bool>();
    Configuration::Get().seed = result["seed"].as<unsigned int>();
    Configuration::Get().numSamples = result["keys"].as<unsigned int>();
    Configuration::Get().penaltyFactorMissingKeys = result["penaltyMissingKeys"].as<double>();
    Configuration::Get().penaltyFactorEmptyPartitions = result["penaltyEmptyPartitions"].as<double>();
    Configuration::Get().repetitions = result["repetitions"].as<unsigned int>();
    Configuration::Get().partition_capacity_inner = result["partition_capacity"].as<size_t>();
    Configuration::Get().entry_capacity_leaf = result["entry_capacity"].as<size_t>();
    Configuration::Get().quantile = result["quantile"].as<double>();
    Configuration::Get().resultFile = result["resultFile"].as<std::string>();
    Configuration::Get().includePartitions = result["includePartitions"].as<bool>();
    Configuration::Get().includeEntries = result["includeEntries"].as<bool>();
    Configuration::Get().parallel = result["parallel"].as<bool>();
    Configuration::Get().cores = Configuration::Get().parallel ? result["cores"].as<unsigned int>() : 1;
    Configuration::Get().sparse = result["sparse"].as<bool>();
    Configuration::Get().printFinalPopulation = result["printFinalPopulation"].as<bool>();
    Configuration::Get().datasetFile = result["datasetFile"].as<std::string>();
    Configuration::Get().workloadFile = result["workloadFile"].as<std::string>();
    Configuration::Get().load_factor_entries_leaves = result["loadFactorEntriesLeaves"].as<double>();
    Configuration::Get().load_factor_partitions_inner = result["loadFactorPartitionsInner"].as<double>();

    Configuration::Get().partitionPriorities[PartitionType::SoAPartition] = result["prioritySoA"].as<double>();
    Configuration::Get().partitionPriorities[PartitionType::TreePartition] = result["priorityTree"].as<double>();
    Configuration::Get().partitionPriorities[PartitionType::HashPartition] = result["priorityHash"].as<double>();

    Configuration::Get().searchPriorities["LinearSearch"] = result["priorityLinear"].as<double>();
    Configuration::Get().searchPriorities["BinarySearch"] = result["priorityBinary"].as<double>();
    Configuration::Get().searchPriorities["DefaultSearch"] = result["priorityDefault"].as<double>();
    Configuration::Get().searchPriorities["InterpolationSearch"] = result["priorityInterpolation"].as<double>();
    Configuration::Get().searchPriorities["ExponentialSearch"] = result["priorityExponential"].as<double>();
    Configuration::Get().searchPriorities["LinearRegressionSearch"] = result["priorityLinearRegression"].as<double>();

    Configuration::Get().saveDotFiles = !(result["noDot"].as<bool>());
    Configuration::Get().printOutput = !(result["minimalOutput"].as<bool>());

    Configuration::Get().progressBarWidth = 50;

    Key lowerBound = Key(1);
    Key upperBound = std::numeric_limits<Key>::max()-1;

    printSettings();

    Dataset<Key, Value>* dataset;
    Workload<Key, Value>* workload;
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
        workload->add_point_uniform(100000);
    }

    Generator<Key, Value>* gen = new FullHorizontalPartitionGenerator<Key, Value>(dataset, workload);

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

    auto startPop = std::chrono::system_clock::now();
    std::time_t tPop = std::chrono::system_clock::to_time_t(startPop);
    std::cout << "Starting population creation at: " << std::ctime(&tPop) << std::endl;
    Population<Key, Value>* pop = new Population<Key, Value>(gen);
    std::pair<const Individual<Key, Value>*, double> p = pop->tournamentSelection(pop->size(), 0.5);
    const HorizontalPartitionIndividual<Key, Value>* bestIndividual = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(p.first);
    double averageFitness = p.second;
    if (Configuration::Get().saveDotFiles)
        bestIndividual->toDotFile("bestIndividualGeneration" + std::to_string(0) + ".dot", Configuration::Get().includePartitions, Configuration::Get().includeEntries);
    if (Configuration::Get().printOutput) {
        std::cout << "Initial population created" << std::endl;
        std::cout << "Best individual of the initial population: " << std::endl << std::endl;
        std::cout << bestIndividual->describe(csvPerLevel) << std::endl << std::endl;
    } else
        bestIndividual->describe(csvPerLevel);
    csv << "0" << "," << bestIndividual->getID() << "," << bestIndividual->getFitness() << "," << averageFitness << "," << 0 << std::endl;

    std::set<Mutation<Key, Value>*> mutations;
    
    /* Basic mutations */
    mutations.insert(new PartitionTypeChangeMutation<Key, Value>(1.0, lowerBound, upperBound));
    mutations.insert(new PartitionSearchChangeMutation<Key, Value>(1.0, lowerBound, upperBound));
    mutations.insert(new EntrySearchChangeMutation<Key, Value>(1.0, lowerBound, upperBound));
    mutations.insert(new HorizontalMergeMutation<Key, Value>(1.0, upperBound));
    mutations.insert(new VerticalMergeMutation<Key, Value>(1.0, upperBound));
    mutations.insert(new HorizontalSplitMutation<Key, Value>(1.0, lowerBound, upperBound));
    mutations.insert(new VerticalSplitMutation<Key, Value>(1.0, lowerBound, upperBound));
    
    std::vector<double> probabilities;
    for(auto mut : mutations)
        probabilities.push_back(mut -> getProbability());
    std::mt19937 randomGen(std::random_device{}());
    std::discrete_distribution<int> dist(probabilities.begin(), probabilities.end());

    unsigned int bestIndividualId = bestIndividual->getID();
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    auto start = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(start);
    std::cout << "Starting mutations at: " << std::ctime(&t) << std::endl;
    printProgressBar(0, Configuration::Get().numGenerations, Configuration::Get().progressBarWidth, Configuration::Get().quantile, averageFitness, r_usage.ru_maxrss);
    for(unsigned int g=0; g<Configuration::Get().numGenerations; g++) {
        p = pop->tournamentSelection(Configuration::Get().tournamentSelectionSize, Configuration::Get().quantile);
        const Individual<Key, Value>* sample = p.first;
        averageFitness = p.second;
        
        /* Multi-threaded case */
        if (Configuration::Get().parallel) {
            std::vector<std::future<const Individual<Key, Value>*>> futures(Configuration::Get().numMutantsPerGeneration);
            for(unsigned int m=0; m<Configuration::Get().numMutantsPerGeneration; m++) {
                int randomChoice = dist(randomGen);
                auto it = mutations.begin();
                std::advance(it, randomChoice);
                futures[m] = tp.enqueue_task(&Mutation<Key, Value>::mutate, *it, sample);
            }
            std::vector<const Individual<Key, Value>*> indToAdd, indToDelete;
            for(unsigned int m=0; m<Configuration::Get().numMutantsPerGeneration; m++) {
                auto mutant = futures[m].get();
                if (DEBUG)
                    std::cout << mutant->describe() << std::endl;
                if (mutant == sample) {
                    continue;
                }
                if ((Configuration::Get().maxFitness ? mutant->getFitness() >= averageFitness : mutant->getFitness() <= averageFitness) || Configuration::Get().quantile == 0.0) {
                    indToAdd.push_back(mutant);
                } else {
                    indToDelete.push_back(mutant);
                }
            }

            for (auto it = indToAdd.begin(); it != indToAdd.end(); ) {
                pop->add(*it);
                it = indToAdd.erase(it);
            }
            assert(indToAdd.empty());
            for (auto it = indToDelete.begin(); it != indToDelete.end(); ) {
                delete(*it);
                it = indToDelete.erase(it);
            }
            assert(indToDelete.empty());
        } 
        /* Single-threaded case */
        else {
            std::vector<const Individual<Key, Value>*> mutants(Configuration::Get().numMutantsPerGeneration);
            for(unsigned int m=0; m<Configuration::Get().numMutantsPerGeneration; m++) {
                int randomChoice = dist(randomGen);
                auto it = mutations.begin();
                std::advance(it, randomChoice);
                mutants[m] = (*it)->mutate(sample);
            }
            for (auto it = mutants.begin(); it != mutants.end(); ) {
                if (DEBUG)
                    std::cout << (*it)->describe() << std::endl;
                if ((*it) == sample) {
                    if (DEBUG)
                    std::cout << "Ignoring duplicate: " << (*it)->getID() << std::endl;
                    it = mutants.erase(it);
                } else if ((Configuration::Get().maxFitness ? (*it)->getFitness() >= averageFitness : (*it)->getFitness() <= averageFitness) || Configuration::Get().quantile == 0.0) {
                    if (DEBUG)
                    std::cout << "Adding to population: " << (*it)->getID() << std::endl;
                    pop -> add(*it);
                    it = mutants.erase(it);
                    if (DEBUG)
                        std::cout << "New population size: " << pop->size() << std::endl;
                } else {
                    if (DEBUG)
                        std::cout << "Directly deleting individual: " << (*it)->getID() << std::endl;
                    delete(*it);
                    it = mutants.erase(it);
                }
            }
        }

        getrusage(RUSAGE_SELF, &r_usage);
        printProgressBar(g+1, Configuration::Get().numGenerations, Configuration::Get().progressBarWidth, Configuration::Get().quantile, averageFitness, r_usage.ru_maxrss);

        p = pop -> tournamentSelection(pop -> size(), 0.5);
        const HorizontalPartitionIndividual<Key, Value>* bestIndividual = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(p.first);

        if (bestIndividualId != bestIndividual->getID()) {
            if (Configuration::Get().saveDotFiles)
                bestIndividual->toDotFile("bestIndividualGeneration" + std::to_string(g+1) + ".dot", Configuration::Get().includePartitions, Configuration::Get().includeEntries);
            bestIndividualId = bestIndividual->getID();
            if (Configuration::Get().printOutput) {
                std::cout << std::endl << "Best individual of generation: " << g+1 << std::endl << std::endl;
                std::cout << bestIndividual->describe(csvPerLevel) << std::endl << std::endl;
            } else
                bestIndividual->describe(csvPerLevel);
        }
        std::chrono::duration<double> duration = std::chrono::system_clock::now() - start;
        auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
        csv << g+1 << "," << bestIndividual->getID() << "," << bestIndividual->getFitness() << "," << p.second << "," << diff << std::endl;
        
        /* Early abort */
        auto bestHp = bestIndividual->getHorizontalPartition();
        if (bestHp->get_entries()->size() == dataset->size() && bestHp->getPartitionType() == PartitionType::SoAPartition && bestHp->get_entry_search_method()->getSearchType() == SearchType::InterpolationSearchMethod)
            break;
    }

    std::cout << std::endl << std::endl;
    std::chrono::duration<double> duration = std::chrono::system_clock::now() - start;
    std::cout << "Finished in: ";
    printDuration(std::cout, duration);
    std::cout << std::endl << std::endl;

    p = pop -> tournamentSelection(pop -> size(), Configuration::Get().quantile);
    const HorizontalPartitionIndividual<Key, Value>* bestIndividualFinalPopulation = dynamic_cast<const HorizontalPartitionIndividual<Key, Value>*>(p.first);
    averageFitness = p.second;
    if (Configuration::Get().printOutput) {
        std::cout << "Best individual of the final population:" << std::endl << std::endl;
        std::cout << bestIndividualFinalPopulation->describe() << std::endl << std::endl;
        //std::cout << "Size of the final population: " << pop->size() << ", average fitness: " << averageFitness << std::endl << std::endl;
        if (Configuration::Get().printFinalPopulation)
            std::cout << "Final population: " << *pop << std::endl;

        std::cout << "Mutation execution counts:" << std::endl;
        for(Mutation<Key, Value>* mut : mutations) {
            std::cout << mut->getName() << " : " << mut->getExecutionCounter() << std::endl;
        }
        std::cout << std::endl;
    }

    csv.close();
    delete(gen);
    delete(pop);
    while(!mutations.empty()) {
        auto it = mutations.begin();
        mutations.erase(it);
        delete(*it);
    }
    delete(dataset);
    delete(workload);
    return 0;
}
