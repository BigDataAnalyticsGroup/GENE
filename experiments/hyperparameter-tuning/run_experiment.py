#!/usr/bin/python3

import subprocess as sp
import pandas as pd
import itertools as it
import numpy as np
from datetime import datetime
from pathlib import Path
import os

pathToProgram = "../../build/release/bin/main_grid_search"

# Build in release mode
process = sp.run(["ninja",  "-C", "../../build/release/"])
process.check_returncode()

Path("./results").mkdir(parents=True, exist_ok=True)

generations = [400]
mutants = [10, 50, 250]
population = [100, 500, 2000]
tournament = [0.1, 0.5, 1.0]
initial = [10, 100]
keys = [100000]
maximize = [False]
seed = [1]
penaltyMissingKeys = [0]
penaltyEmptyPartitions = [0]
repetitions = [5]
partition_capacity = keys
entry_capacity = keys
quantile = [0.0, 0.25, 0.5, 0.75, 1.0]
resultFile = ["results/results.csv"]
includePartitions = [False]
includeEntries = [False]
parallel = [True]
cores = [15]
sparse = [False]
printFinalPopulation = [False]
loadFactorEntriesLeaves = [0.1]
loadFactorPartitionsInner = [0.00005]

runs = np.arange(5)

columns = ["generations", "mutants", "population", "tournament", "initial",
           "keys", "maximize", "seed", "penaltyMissingKeys",
           "penaltyEmptyPartitions", "repetitions", "partition_capacity",
           "entry_capacity", "quantile",
           "resultFile", "includePartitions", "includeEntries", "parallel",
           "cores", "sparse", "printFinalPopulation",
           "loadFactorEntriesLeaves", "loadFactorPartitionsInner",
           "Run", "Generation", "Time"]
dtypes = [np.uint, np.uint, np.uint, np.uint, np.uint, np.uint, np.bool,
          np.int, np.double, np.double, np.uint, np.uint, np.uint,
          np.double, pd.StringDtype(), np.bool, np.bool,
          np.bool, np.uint, np.bool, np.bool, np.double, np.double, np.uint, np.uint, np.uint]

runtimes = pd.DataFrame()
for c, d in zip(columns, dtypes):
    runtimes[c] = pd.Series(dtype=d)

num_combinations = len(generations) * len(mutants) * len(population) * len(tournament) * len(initial) * \
    len(keys) * len(maximize) * len(seed) * len(penaltyMissingKeys) * \
    len(penaltyEmptyPartitions) * len(repetitions) * len(partition_capacity) * len(entry_capacity) * \
    len(quantile) * len(resultFile) * \
    len(includePartitions) * len(includeEntries) * len(parallel) * len(cores) * len(sparse) * \
    len(printFinalPopulation) * \
    len(loadFactorEntriesLeaves) * len(loadFactorPartitionsInner) * len(runs)
combinations = it.product(generations, mutants, population, tournament, initial, keys, maximize, seed,
                          penaltyMissingKeys, penaltyEmptyPartitions, repetitions,
                          partition_capacity, entry_capacity, quantile,
                          resultFile, includePartitions, includeEntries, parallel, cores, sparse,
                          printFinalPopulation, loadFactorEntriesLeaves,
                          loadFactorPartitionsInner, runs)

counter = 0
start = datetime.now()
print("Starting exploration of", num_combinations, "possible settings at", start)
for setting in combinations:
    counter += 1
    print("Combination:", counter, "/",
          num_combinations, "Current setting:", setting)

    arguments = []
    arguments.append("-g")
    arguments.append(str(setting[0]))
    arguments.append("-m")
    arguments.append(str(setting[1]))
    arguments.append("-p")
    arguments.append(str(setting[2]))
    arguments.append("-t")
    arguments.append(str(int(setting[3] * setting[2])))
    arguments.append("-i")
    arguments.append(str(setting[4]))
    arguments.append("-k")
    arguments.append(str(setting[5]))
    if setting[6]:
        arguments.append("--maximize")
    arguments.append("-s")
    arguments.append(str(setting[7]))
    arguments.append("--penaltyMissingKeys")
    arguments.append(str(setting[8]))
    arguments.append("--penaltyEmptyPartitions")
    arguments.append(str(setting[9]))
    arguments.append("-r")
    arguments.append(str(setting[10]))
    arguments.append("--partition_capacity")
    arguments.append(str(setting[11]))
    arguments.append("--entry_capacity")
    arguments.append(str(setting[12]))
    arguments.append("-q")
    arguments.append(str(setting[13]))
    arguments.append("--resultFile")
    arguments.append(str(setting[14]))
    if setting[15]:
        arguments.append("--includePartitions")
    if setting[16]:
        arguments.append("--includeEntries")
    if setting[17]:
        arguments.append("--parallel")
        arguments.append("-c")
        arguments.append(str(setting[18]))
    if setting[19]:
        arguments.append("--sparse")
    if setting[20]:
        arguments.append("--printFinalPopulation")
    arguments.append("--loadFactorEntriesLeaves")
    arguments.append(str(setting[21]))
    arguments.append("--loadFactorPartitionsInner")
    arguments.append(str(setting[22]))
    arguments.append("--noDot")

    process = sp.run([pathToProgram] + arguments,
                     stdout=sp.PIPE, universal_newlines=True)
    process.check_returncode()

    resultCsv = pd.read_csv(resultFile[0], header=0)
    lastGeneration = resultCsv.tail(1)

    row = {}
    for i in range(24):
        row[columns[i]] = setting[i]
    row[columns[-2]] = lastGeneration["Generation"].values[0]
    row[columns[-1]] = lastGeneration["Runtime"].values[0]
    runtimes = runtimes.append(row, ignore_index=True)
    runtimes.to_csv("results/grid_search.csv")
    os.remove(resultFile[0])
    os.remove(resultFile[0].replace(".csv", "perLevel.csv"))
    print("Finished setting", "runtime:",
          row[columns[-1]], "ms", "generations:", row[columns[-2]])

end = datetime.now()
print("Finished exploration at", end, "\nTotal Runtime:", end-start)

runtimes.to_csv("results/grid_search.csv")
median_runtimes = runtimes.groupby(
    columns[:-3]).median().sort_values(columns[-1])
print("\nTop Median Runtimes:\n", median_runtimes[columns[-1]].head(5))
mean_runtimes = runtimes.groupby(columns[:-3]).mean().sort_values(columns[-1])
print("\nTop Mean Runtimes:\n", mean_runtimes[columns[-1]].head(5))
