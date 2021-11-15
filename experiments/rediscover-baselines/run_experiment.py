#!/usr/bin/python3

import subprocess as sp
from pathlib import Path
import sys
import os

RELEASE_BUILD = "../../build/release"
RELEASE_BUILD_ROOT = "./build/release"
EXPERIMENT_FOLDER = "experiments/rediscover-baselines/"
pathToDataGeneration = "/bin/generate_data"
pathToWorkloadGeneration = "/bin/generate_workload"
pathToProgram = RELEASE_BUILD + "/bin/main"
pathToBtreeBaseline = RELEASE_BUILD + "/bin/btree_baseline"
pathToHashtableBaseline = RELEASE_BUILD + "/bin/hashtable_baseline"
pathToMixBaseline = RELEASE_BUILD + "/bin/btree_baseline"
baselinePrograms = [pathToBtreeBaseline,
                    pathToHashtableBaseline, pathToMixBaseline]

# Build in release mode
process = sp.run(["ninja",  "-C", "../../build/release/"])
process.check_returncode()

# dataset, workload & result files
datasetFileUniDense = "data/uni-dense_100000.data"
datasetFileBooks = "data/books_100000.data"
datasetFilesScalingBooks = ["data/books_1000000.data",
                            "data/books_10000000.data"]
scalingStepsize = 1

workloadFileBtreeUniDense = "workloads/range_uni-dense_100000_10000_0.000000_0.010000.wkl"
workloadFileHashtableUniDense = "workloads/point_uni-dense_100000_10000_1.000000_0.010000.wkl"
workloadFileMixUniDense = "workloads/mix_uni-dense_100000_10000_0.800000_0.010000.wkl"
workloadFilesUniDense = [workloadFileBtreeUniDense,
                         workloadFileHashtableUniDense, workloadFileMixUniDense]
workloadFileBtreeBooks = "workloads/range_books_100000_10000_0.000000_0.010000.wkl"
workloadFileHashtableBooks = "workloads/point_books_100000_10000_1.000000_0.010000.wkl"
workloadFileMixBooks = "workloads/mix_books_100000_10000_0.800000_0.010000.wkl"
workloadFilesBooks = [workloadFileBtreeBooks,
                      workloadFileHashtableBooks, workloadFileMixBooks]
workloadPrefixes = ["btree", "hashtable", "mix"]

resultFileBtreeUniDense = "results/uni-dense_btree_100000.csv"
resultFileHashtableUniDense = "results/uni-dense_hashtable_100000.csv"
resultFileMixUniDense = "results/uni-dense_mix_100000.csv"
resultFilesUniDense = [resultFileBtreeUniDense,
                       resultFileHashtableUniDense, resultFileMixUniDense]
resultFileBtreeBooks = "results/books_btree_100000.csv"
resultFileHashtableBooks = "results/books_hashtable_100000.csv"
resultFileMixBooks = "results/books_mix_100000.csv"
resultFilesBooks = [resultFileBtreeBooks,
                    resultFileHashtableBooks, resultFileMixBooks]

# generate data & workloads
Path("./data").mkdir(parents=True, exist_ok=True)
Path("./workloads").mkdir(parents=True, exist_ok=True)
Path("./results").mkdir(parents=True, exist_ok=True)

os.chdir("../..")
process = sp.run([RELEASE_BUILD_ROOT + pathToDataGeneration,
                  "uni-dense", "100000", EXPERIMENT_FOLDER + datasetFileUniDense])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToWorkloadGeneration, EXPERIMENT_FOLDER + datasetFileUniDense,
                  "10000", "0", "0.01", EXPERIMENT_FOLDER + workloadFileBtreeUniDense])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToWorkloadGeneration, EXPERIMENT_FOLDER + datasetFileUniDense,
                  "10000", "1", "0.01", EXPERIMENT_FOLDER + workloadFileHashtableUniDense])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToWorkloadGeneration, EXPERIMENT_FOLDER + datasetFileUniDense,
                  "10000", "0.8", "0.01", EXPERIMENT_FOLDER + workloadFileMixUniDense])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToDataGeneration,
                  "books", "100000", EXPERIMENT_FOLDER + datasetFileBooks])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToDataGeneration, "books",
                  "1000000", EXPERIMENT_FOLDER + datasetFilesScalingBooks[0]])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToDataGeneration, "books",
                  "10000000", EXPERIMENT_FOLDER + datasetFilesScalingBooks[1]])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToWorkloadGeneration, EXPERIMENT_FOLDER + datasetFileBooks,
                  "10000", "0", "0.01", EXPERIMENT_FOLDER + workloadFileBtreeBooks])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToWorkloadGeneration, EXPERIMENT_FOLDER + datasetFileBooks,
                  "10000", "1", "0.01", EXPERIMENT_FOLDER + workloadFileHashtableBooks])
process.check_returncode()
process = sp.run([RELEASE_BUILD_ROOT + pathToWorkloadGeneration, EXPERIMENT_FOLDER + datasetFileBooks,
                  "10000", "0.8", "0.01", EXPERIMENT_FOLDER + workloadFileMixBooks])
process.check_returncode()
os.chdir(EXPERIMENT_FOLDER)

# Set genetic program parameters
generations = [8000, 8000, 8000]    # number of generations g_max
mutants = 10                        # number of mutants per generation s_max
population = 50                     # max size of population
tournament = 25                     # sample size (in tournament selection)
initial = 20                        # initial size of population
keys = 100000                       # dataset size
# maximize fitness value (we minimize the runtime)
maximize = False
seed = 42                           # seed for random generator
# num of repetitions (median fitness value is taken)
repetitions = 5
partition_capacity_inner = 100000             # max partition capacity
partition_capacity_leaf = 0             # max partition capacity
entry_capacity_inner = 0               # max entry_capacity
entry_capacity_leaf = 100000              # max entry_capacity

quantile = 0.5                      # quantile to be reached for population insertion
includePartitions = True           # in dot files
includeEntries = False              # in dot files
parallel = True                     # multi-threading
cores = 2                           # num cores (capped at num mutants)
printFinalPopulation = False        # print all individuals at the end
loadFactorEntriesLeaves = 0.01      # load factor for leaves
loadFactorPartitionsInner = 1e-4  # load factor for inner nodes
prioritySoA = 4                     # specify priorities
priorityTree = 1
priorityHash = 4
priorityLinear = 0
priorityBinary = 2
priorityDefault = 2
priorityInterpolation = 2
priorityExponential = 1
priorityLinearRegression = 0
noDot = False
minimalOutput = True
tracing = False
trainTestSplit = False
chainMutationLength = 0

# Set baseline parameters
btree_fanout = 10
btree_leafsize = 1000
hashtable_size = 100000
split_num_leaves = 11

# Assert that a workload and result file is defined for each baseline
assert(len(baselinePrograms) == len(workloadFilesUniDense))
assert(len(baselinePrograms) == len(resultFilesUniDense))
assert(len(baselinePrograms) == len(workloadFilesBooks))
assert(len(baselinePrograms) == len(resultFilesBooks))

# Loop over baselines, execute genetic as well as baseline program for uni-dense dataset
for i in [0, 1, 2]:

    # Define baseline arguments
    if i == 0:
        # btree baseline
        arguments = [
            datasetFileUniDense,
            workloadFilesUniDense[i],
            str(btree_fanout),
            str(int(btree_leafsize)),
            str(repetitions),
            "results/uni-dense_" + workloadPrefixes[i] + "_",
        ]
    elif i == 1:
        # hashtable baseline
        arguments = [
            datasetFileUniDense,
            workloadFilesUniDense[i],
            str(int(hashtable_size)),
            str(repetitions),
            "results/uni-dense_" + workloadPrefixes[i] + "_",
        ]
    elif i == 2:
        # mix baseline
        arguments = [
            datasetFileUniDense,
            workloadFilesUniDense[i],
            str(btree_fanout),
            str(int(btree_leafsize)),
            str(repetitions),
            "results/uni-dense_" + workloadPrefixes[i] + "_",
        ]

    print("\nExecuting baseline algorithm for workload",
          i, "\n")

    # Run baseline experiment
    process = sp.run([baselinePrograms[i]] + arguments)
    process.check_returncode()

    # Define genetic arguments
    arguments = [
        "-g", str(generations[i]),
        "-m", str(mutants),
        "-p", str(population),
        "-t", str(tournament),
        "-i", str(initial),
        "-k", str(keys),
        "-s", str(seed),
        "-r", str(repetitions),
        "--partition_capacity_inner", str(
            int(partition_capacity_inner)),
        "--partition_capacity_leaf", str(partition_capacity_leaf),
        "--entry_capacity_inner", str(entry_capacity_inner),
        "--entry_capacity_leaf", str(
            int(entry_capacity_leaf)),
        "-q", str(quantile),
        "--resultFile", resultFilesUniDense[i],
        "-c", str(cores),
        "--datasetFile", datasetFileUniDense,
        "--workloadFile", workloadFilesUniDense[i],
        "--loadFactorEntriesLeaves", str(loadFactorEntriesLeaves),
        "--loadFactorPartitionsInner", str(
            loadFactorPartitionsInner),
        "--prioritySoA", str(prioritySoA),
        "--priorityTree", str(priorityTree),
        "--priorityHash", str(priorityHash),
        "--priorityLinear", str(priorityLinear),
        "--priorityBinary", str(priorityBinary),
        "--priorityDefault", str(priorityDefault),
        "--priorityInterpolation", str(priorityInterpolation),
        "--priorityExponential", str(priorityExponential),
        "--priorityLinearRegression", str(priorityLinearRegression),
        "--outputPrefix", "results/uni-dense_" + workloadPrefixes[i] + "_",
        "--chainMutationLength", str(chainMutationLength),
    ]
    if (maximize):
        arguments.append("--maximize")
    if (includePartitions):
        arguments.append("--includePartitions")
    if (includeEntries):
        arguments.append("--includeEntries")
    if (parallel):
        arguments.append("--parallel")
    if (printFinalPopulation):
        arguments.append("--printFinalPopulation")
    if (noDot):
        arguments.append("--noDot")
    if (minimalOutput):
        arguments.append("--minimalOutput")
    if (tracing):
        arguments.append("--tracing")
    if (not trainTestSplit):
        arguments.append("--disableTrainTestSplit")

    print("\nExecuting genetic algorithm for workload",
          i, "\n")

    # Run genetic experiment
    process = sp.run([pathToProgram] + arguments)
    process.check_returncode()

# Loop over baselines, execute genetic as well as baseline program for books dataset
for i in [0, 1, 2]:

    # Define baseline arguments
    if i == 0:
        # btree baseline
        arguments = [
            datasetFileBooks,
            workloadFilesBooks[i],
            str(btree_fanout),
            str(int(btree_leafsize)),
            str(repetitions),
            "results/books_" + workloadPrefixes[i] + "_",
        ]
    elif i == 1:
        # hashtable baseline
        arguments = [
            datasetFileBooks,
            workloadFilesBooks[i],
            str(int(hashtable_size)),
            str(repetitions),
            "results/books_" + workloadPrefixes[i] + "_",
        ]
    elif i == 2:
        # mix baseline
        arguments = [
            datasetFileBooks,
            workloadFilesBooks[i],
            str(btree_fanout),
            str(int(btree_leafsize)),
            str(repetitions),
            "results/books_" + workloadPrefixes[i] + "_",
        ]

    print("\nExecuting baseline algorithm for workload",
          i, "\n")

    # Run baseline experiment
    process = sp.run([baselinePrograms[i]] + arguments)
    process.check_returncode()

    # Define genetic arguments
    arguments = [
        "-g", str(generations[i]),
        "-m", str(mutants),
        "-p", str(population),
        "-t", str(tournament),
        "-i", str(initial),
        "-k", str(keys),
        "-s", str(seed),
        "-r", str(repetitions),
        "--partition_capacity_inner", str(
            int(partition_capacity_inner)),
        "--partition_capacity_leaf", str(partition_capacity_leaf),
        "--entry_capacity_inner", str(entry_capacity_inner),
        "--entry_capacity_leaf", str(
            int(entry_capacity_leaf)),
        "-q", str(quantile),
        "--resultFile", resultFilesBooks[i],
        "-c", str(cores),
        "--datasetFile", datasetFileBooks,
        "--workloadFile", workloadFilesBooks[i],
        "--loadFactorEntriesLeaves", str(loadFactorEntriesLeaves),
        "--loadFactorPartitionsInner", str(
            loadFactorPartitionsInner),
        "--prioritySoA", str(prioritySoA),
        "--priorityTree", str(priorityTree),
        "--priorityHash", str(priorityHash),
        "--priorityLinear", str(priorityLinear),
        "--priorityBinary", str(priorityBinary),
        "--priorityDefault", str(priorityDefault),
        "--priorityInterpolation", str(priorityInterpolation),
        "--priorityExponential", str(priorityExponential),
        "--priorityLinearRegression", str(priorityLinearRegression),
        "--outputPrefix", "results/books_" + workloadPrefixes[i] + "_",
        "--chainMutationLength", str(chainMutationLength),
        "--scalingStepsize", str(scalingStepsize),
    ]
    if len(datasetFilesScalingBooks) > 0:
        arguments.append("--scalingDatasets")
        arguments.append(",".join(datasetFilesScalingBooks))
    if (maximize):
        arguments.append("--maximize")
    if (includePartitions):
        arguments.append("--includePartitions")
    if (includeEntries):
        arguments.append("--includeEntries")
    if (parallel):
        arguments.append("--parallel")
    if (printFinalPopulation):
        arguments.append("--printFinalPopulation")
    if (noDot):
        arguments.append("--noDot")
    if (minimalOutput):
        arguments.append("--minimalOutput")
    if (tracing):
        arguments.append("--tracing")
    if (not trainTestSplit):
        arguments.append("--disableTrainTestSplit")

    print("\nExecuting genetic algorithm for workload",
          i, "\n")

    # Run genetic experiment
    process = sp.run([pathToProgram] + arguments)
    process.check_returncode()
