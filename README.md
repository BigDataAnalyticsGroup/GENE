# GENE

This repository contains code used for our VLDB paper ["The next 50 Years in Database Indexing or: The Case for Automatically Generated Index Structures"](https://arxiv.org/abs/2009.10669) (Pre-print).

## Download & Prerequisites
To use this code, we recommend to clone the repository using `git`.  You should use the `--recurse-submodules` flag in
your `git clone` command to automatically download the submodules while cloning the repository.
```sh
git clone --recurse-submodules git@github.com:BigDataAnalyticsGroup/GENE.git
```

The execution and visualization of experiments is based on the following tools:
* C++ compiler supporting C++17
* `bash>=4`
* `cmake>=3.1`
* `Python>=3.5`
* `md5sum`
* `wget`
* `zstd`
* `graphviz` (optional)

To install the required Python modules, you can use Python's package installer `pip` in combination with the
`requirements.txt`.
```sh
pip install -r requirements.txt
```

To download the data required by the experiments, run the following shell script **from the root folder** of the project.
```sh
./data/download.sh
```

## Build
Our build system is based on `cmake`.  See the following instructions for an example build process in `release` mode.

**NOTE:** The following example requires the execution from the root folder. In addition, reproducing our experimental
results with the provided scripts requires exactly this folder structure.
```sh
mkdir -p build/release
cd build/release
cmake -DCMAKE_BUILD_TYPE=Release ../..
make
cd ../..
```

## Executing Genetic Searches
After executing the build instructions, the generated binaries are located in `./build/release/bin/`.  The executable
for the genetic search is called `main`. It offers a multitude of different parameters to specify nearly all important
hyperparameters and inputs. Calling the executable with the `-h` or `--help` flag will display all available options
With the default settings, the executable will run a minimal example with 100 keys and 10 generations, which allows to
test if everything works as expected.
```sh
./build/release/bin/main
```

## Experiments
In the `./experiments/` folder, we provide scripts to reproduce the experimental results form our paper. In particular,
we provide the following experiments.

* `hyperparameter-tuning`: Section 6.1 "Hyperparameter Tuning", computes the best hyperparameters for our genetic
  search.
* `rediscover-baselines`: Section 6.2 "Rediscover Suitable Baseline Indexes", demonstrates that our genetic algorithm is
  capable of reproducing the performance of various basline indexes.
* `optimized-vs-heuristic`: Section 6.3 "Optimized vs Heuristic Indexes", compares the performance of GENE with
  representatives of different prevalent heuristic index types.

To run the different experiments, navigate to the corresponding folder and execute the Python
`run_experiment.py` or Shell `run_experiment.sh` script inside the respective folder. For example, to run the
hyperparameter search:
```sh
cd experiments/hyperparameter-tuning/
python3 run_experiment.py
```
Depending on the specific experiment, the script generates the following folders:
* `data`: contains the underlying dataset files
* `workloads`: contains the underlying workload files
* `results`: contains the corresponding result files in `csv` and/or `dot` format

To visualize the results, execute the accompanying `visualization.py` script (except for the `hyperparameter-tuning`
experiment). This will produce `pdf` files containing the plots as shown in the paper.
To export the dot files using `graphviz` as `pdf`, execute the following command:
```sh
dot -Tpdf path/to/dot/file -o /path/to/pdf/file
```
