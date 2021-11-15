# GENE

This repository contains code used for our VLDB paper ["The next 50 Years in Database Indexing or: The Case for Automatically Generated Index Structures"](https://arxiv.org/abs/2009.10669) (Pre-print).

## Download & Prerequisites

To use this code, simply clone the repository using git.
As we use submodules to include code for the [TLX-tree](https://github.com/tlx/tlx), a [Robin-Hood-map](https://github.com/Tessil/robin-map) and others.
You should use the `--recurse-submodules` flag in your `git clone` command to automatically download the submodules while cloning our repository.

To automate the execution of experiments and to plot the results, we used scripts based on Python 3 and common tools such as jupyter notebooks and libraries such as matplotlib.
To use the available scripts, you will therefore need to install Python 3 and the necessary packages which are all available via Python's package manager pip.

## Building

We us the [Ninja build system](https://ninja-build.org) to automate code generation in our project.
After installing ninja and cmake on your local machine, you should create folders to store the binary targets in.
We recommend a folder structure such as "build/debug" and "build/release" to receive a single folder containing all binary targets, but separated by their build type.
After entering the folder for your binary targets, simply execute the following command:
`cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug path/to/repository/root/folder` (for the debug build, substiture `Debug` by `Release` for a release build)
If you want to specify a compiler of your choice, you can do so by adding the flags `DCMAKE_C_COMPILER=path/to/c/compiler` or `DCMAKE_CXX_COMPILER=path/to/c++/compiler`.
We recommend using [LLVM's clang compiler family](https://clang.llvm.org).
To finally build the code, simply execute `ninja` in the corresponding `build/debug` or `build/release` folder.
Alternatively, you can also specify which folder to build by specifying the corresponding flag: `ninja -C path/to/folder` will build the corresponding folder given by the path.

## Executing Genetic Searches

After building the 
The main executable of our repository is the so called `main_genetic_general` and should be available under `build/release/bin` (assuming the default folder structure).
It offers a multitude of different parameters to specify nearly all important hyperparameters and inputs.
Calling the executable with the `-h` or `--help` flag will display all available options with a short description.
With the default settings, the executable will run a minimal example with a small amount of keys and generations which allows to test if everything works as expected.

## Reproducing Experiments

In the experiments folder, we offer different Python scripts to reproduce the experiments shown in our paper.
Note however that the paths in these scripts assume our recommended folder structure, especially the existence of a `build/release` folder containing the executable compiled in release mode.
If you used a different build system or a different folder structure, you might have to adapt the corresponding paths in the scripts before you can execute them.

The `experiments/hyperparameter-tuning` folder contains a script to conduct the hyperparameter search as described in Section 6.1 of our paper.
The `experiments/rediscover-baselines` folder contains the script to reproduce the experiment "Rediscover Suitable Baseline Indexes" as described in Section 6.2 of our paper.
The `experiments/optimized-vs-heuristic` folder contains a shell script to reproduce the experiment "Optimized vs.
Heuristic Indexes" as described in Section 6.3 of our paper. [NOTE: TOOD]

All these scripts will produce several outputs by default:
1) A CSV-file containing the results achieved in each generation of the genetic algorithm. The exact names of these files depend on the experiment and are specified in the corresponding script as `resultFile`.
2) A CSV-file containing a textual description of the best individuals found during the search process. The exact names of these files again depend on the concrete experiment and consist of the name stored in the `resultFile` filed extended by the suffix `perLevel`.
Note: We only store this textual representation once per individual found and do not repeat it several times if the corresponding individual remains the best one over multiple generations of the search.
You can however join the two CSV-files on the ID of the individual if you need information about the best individual in a specific generation.
3) Graphviz Dot-files containing a visual representation of the best individuals found during search. 
Each individual is stored exactly once, namely in the first generation it appeared as best performing individual.
The file names have always the form `bestIndividualGenerationX.dot` with `X` being the generation this individual first appeared.
Depending on the concrete experiment, the file names might be prefixed to distinguish different runs performed by the same Python script.
You can again match that information with the CSV files by scanning the CSV-file for the first generation a specific individual appeared as the best performing one.
To visualize the dot files, you first need to install graphviz on your local machine.
You can then create a PDF-file containing the visualization using the following command: `dot -Tpdf path/to/dot/file -o /path/to/pdf/file`.

To create the visualizations of the different experiments (e.g. Figure 6), we then used simple jupyter notebooks and the Python matplotlib library.
To read the CSV-files, we relied on the pandas and numpy libraries.
All of these libraries were freely available via Python's package manager pip at the time of writing.

