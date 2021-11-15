#pragma once

#include "generator.hpp"
#include "individual.hpp"
#include "util/globals.hpp"
#include <math.h>
#include <ostream>
#include <random>
#include <limits>
#include <cassert>
#include <set>

const bool DEBUG_POPULATION = false;

/**
 * implementation of the population containing a set of individuals ordered by their fitness
 */
template<typename Key, typename Value>
class Population {
    public:
        Population(Generator<Key, Value>* gen) : maxPopulationSize(Configuration::Get().populationSize), maxFitness(Configuration::Get().maxFitness)
        {
            unsigned int initialSize = Configuration::Get().initialPopulationSize;
            if (initialSize > maxPopulationSize)
                throw std::invalid_argument("Initial population size cannot be larger than maximum size");
            population = gen->generateIndividuals(initialSize);
        };

        ~Population() {
            if (DEBUG_POPULATION)
                std::cout << "Population expected to delete " << population.size() << " entries" << std::endl;
            for(auto it = population.begin(); it != population.end(); ) {
                if (DEBUG_POPULATION)
                    std::cout << "Population deleting individual in destructor: " << (*it)->getID() << std::endl;
                delete(*it);
                it = population.erase(it);
            }
        }

        /** perform tournament selection, i.e.
        * - select a random subset of size num
        * - determine the best individual within this subset
        * - determine the fitness quantile within this subset
        * returns a pair consisting of the best individual and the fitness quantile
        */
        std::pair<const Individual<Key, Value>*, double> tournamentSelection(int num, double quantile) {
            std::set<const Individual<Key, Value>*> selection;
            std::sample(population.begin(), population.end(), std::inserter(selection, selection.begin()), num, std::mt19937{std::random_device{}()});
            double bestFitness = maxFitness ? std::numeric_limits<int>::min() : std::numeric_limits<int>::max();
            std::vector<double> fitnessValues;
            const Individual<Key, Value>* bestIndividual = NULL;
            for(const Individual<Key, Value>* elem : selection) {
                double fitness = elem -> getFitness();
                fitnessValues.push_back(fitness);
                if(maxFitness ? fitness > bestFitness : fitness < bestFitness) {
                    bestFitness = fitness;
                    bestIndividual = elem;
                }
            }
            std::vector<double>::iterator b = fitnessValues.begin();
            std::vector<double>::iterator e = fitnessValues.end();
            std::vector<double>::iterator nth = b;
            const std::size_t pos = maxFitness ? quantile * std::distance(b, e) : (1-quantile) * std::distance(b, e);
            std::advance(nth, pos);
            std::nth_element(b, nth, e);
            return std::make_pair(bestIndividual, *nth);
        }

        /**
         * performs a tournament selection as specified before,
         * additionally writes information to the given stringstream
         */
        std::pair<const Individual<Key, Value>*, double> tournamentSelection(int num, double quantile, std::stringstream &ss) {
            std::set<const Individual<Key, Value>*> selection;
            std::sample(population.begin(), population.end(), std::inserter(selection, selection.begin()), num, std::mt19937{std::random_device{}()});
            double bestFitness = maxFitness ? std::numeric_limits<int>::min() : std::numeric_limits<int>::max();
            std::vector<double> fitnessValues;
            const Individual<Key, Value>* bestIndividual = NULL;
            ss << "IDs=(";
            double sum_fitness = 0.0;
            for(const Individual<Key, Value>* elem : selection) {
                ss << elem->getID() << ",";
                double fitness = elem -> getFitness();
                sum_fitness += fitness;
                fitnessValues.push_back(fitness);
                if(maxFitness ? fitness > bestFitness : fitness < bestFitness) {
                    bestFitness = fitness;
                    bestIndividual = elem;
                }
            }
            ss << ")";
            double mean = sum_fitness / fitnessValues.size();
            std::vector<double>::iterator b = fitnessValues.begin();
            std::vector<double>::iterator e = fitnessValues.end();
            std::vector<double>::iterator nth = b;
            const std::size_t pos = maxFitness ? quantile * std::distance(b, e) : (1-quantile) * std::distance(b, e);
            std::advance(nth, pos);
            std::nth_element(b, nth, e);
            std::sort(fitnessValues.begin(), fitnessValues.end());
            double sum_nominator = 0.0;
            for (auto f : fitnessValues) sum_nominator += (f - mean) * (f - mean);
            double stddev_fitness = sqrt((sum_nominator) / fitnessValues.size());
            ss << " | median_fitness=" << *nth
               << " | stddev_fitness=" << stddev_fitness
               << " | min_fitness=" << *fitnessValues.begin()
               << " | max_fitness=" << *std::prev(fitnessValues.end())
               << " | best_ind_ID=" << bestIndividual->getID();
            return std::make_pair(bestIndividual, *nth);
        }

        /** returns the best individual of the population */
        const Individual<Key, Value> * best_individual() {
            if (maxFitness) return *std::prev(population.end());
            else return *population.begin();
        }
        /** returns the worst individual of the population */
        const Individual<Key, Value> * worst_individual() {
            if (not maxFitness) return *std::prev(population.end());
            else return *population.begin();
        }

        /** compute the mean fitness of the population */
        double mean_fitness() {
            double sum_fitness = 0.0;
            for (auto p : population) {
                auto fitness = p->getFitness();
                sum_fitness += fitness;
            }
            return sum_fitness / population.size();
        }


        /** Compute the standard deviation defined as:
         * sqrt(sum(x_i - mean)^2/N) */
        double stddev_fitness() {
            double mean = mean_fitness();
            double sum_nominator = 0.0;
            for (auto p : population) {
                auto fitness = p->getFitness();
                sum_nominator += (fitness - mean) * (fitness - mean);
            }
            return sqrt((sum_nominator) / population.size());
        }

        /**
         * Add the given individual to the population
         * In case the population exceeds its maximum size, remove the worst individual afterwards
         */
        void add(const Individual<Key, Value>* elem) {
            if(not population.insert(elem).second) {
                if (DEBUG_POPULATION)
                    std::cout << "Population deleting individual in add duplicate: " << elem->getID() << std::endl;
                delete(elem);
            }
            while(population.size() >= maxPopulationSize) {
                if(maxFitness) {
                    auto it = population.begin();
                    assert((*it)->getFitness() <= elem->getFitness());
                    if (DEBUG_POPULATION)
                        std::cout << "Population deleting individual in add: " << (*it)->getID() << std::endl;
		            delete((*it));
                    population.erase(it);
                } else {
                    auto it = std::prev(population.end());
                    assert((*it)->getFitness() >= elem->getFitness());
                    if (DEBUG_POPULATION)
                        std::cout << "Population deleting individual in add: " << (*it)->getID() << std::endl;
                    delete((*it));
                    population.erase(it);
                }
            }
        }

        /** return the current population size */
        unsigned int size() { return population.size(); }

        /** Print out the population */
        friend std::ostream& operator << (std::ostream& stream, Population<Key, Value>& pop) {
            stream << "Population of size " << pop.population.size() << " with members:" << std::endl;
            for(const Individual<Key, Value>* elem : pop.population) {
                stream << *elem << std::endl;
            }
            return stream;
        }

    /* Iterator. */
    using iterator = typename std::set<const Individual<Key, Value>*, IndividualFitnessComparator<Key, Value>>::iterator;
    using const_iterator = typename std::set<const Individual<Key, Value>*, IndividualFitnessComparator<Key, Value>>::const_iterator;

    iterator begin() { return population.begin(); }
    iterator end()   { return population.end(); }
    const_iterator cbegin() const { return population.cbegin(); }
    const_iterator cend()   const { return population.cend(); }

    private:
        std::set<const Individual<Key, Value>*, IndividualFitnessComparator<Key, Value>> population;
        unsigned int maxPopulationSize;
        bool maxFitness;
};
