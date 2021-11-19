# import necessary libraries

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# define constants and helping functions

results_folder = "./results/"
figsize = (8, 5)
fontsize_large = 14
fontsize_small = 13
algo_name = 'GENE'
original_dataset_size = 100000
marker = ["-", "--", "-.", ":"]
marker = [(0, ()), (0, (5, 7)), (0, (4, 4, 1, 4)), (0, (1, 1))]

def size_to_string(size):
    if size < 1e3:
        return str(int(size))
    elif size < 1e6:
        return str(int(size // 1e3)) + "K"
    else:
        return str(int(size // 1e6)) + "M"
    
# read results for uni-dense dataset

filename_genetic_point = "uni-dense_point_100000.csv"
filename_baseline_point = "uni-dense_point_baseline.csv"
filename_genetic_mix = "uni-dense_mix_100000.csv"
filename_baseline_mix = "uni-dense_mix_baseline.csv"
filename_genetic_range = "uni-dense_range_100000.csv"
filename_baseline_range = "uni-dense_range_baseline.csv"

results_genetic_point = pd.read_csv(results_folder + filename_genetic_point, header=0)
results_baseline_point = pd.read_csv(results_folder + filename_baseline_point, header=0)
results_genetic_mix = pd.read_csv(results_folder + filename_genetic_mix, header=0)
results_baseline_mix = pd.read_csv(results_folder + filename_baseline_mix, header=0)
results_genetic_range = pd.read_csv(results_folder + filename_genetic_range, header=0)
results_baseline_range = pd.read_csv(results_folder + filename_baseline_range, header=0)

# plot uni-dense with point workload

bestFitness = results_genetic_point["Fitness"].min()
generations = results_genetic_point[results_genetic_point["Fitness"] == bestFitness]["Generation"].min() + 1
median_point = results_baseline_point.median()["Duration"]

fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.axhline(median_point, color="blue", linestyle=":", label="baseline")
ax.plot(results_genetic_point["Generation"].values[:generations], results_genetic_point["Fitness"].values[:generations], color="blue", label=algo_name)
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("uni-dense_point.pdf", bbox_inches="tight")
plt.close(fig)

# plot uni-dense with mix workload

bestFitness = results_genetic_mix["Fitness"].min()
generations = results_genetic_mix[results_genetic_mix["Fitness"] == bestFitness]["Generation"].min() + 1
median_mix = results_baseline_mix.median()["Duration"]

fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.axhline(median_mix, color="red", linestyle=":", label="baseline")
ax.plot(results_genetic_mix["Generation"].values[:generations], results_genetic_mix["Fitness"].values[:generations], color="red", label=algo_name)
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("uni-dense_mix.pdf", bbox_inches="tight")
plt.close(fig)

# plot uni-dense with range workload

bestFitness = results_genetic_range["Fitness"].min()
generations = results_genetic_range[results_genetic_range["Fitness"] == bestFitness]["Generation"].min() + 1
median_range = results_baseline_range.median()["Duration"]

fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.axhline(median_range, color="green", linestyle=":", label="baseline")
ax.plot(results_genetic_range["Generation"].values[:generations], results_genetic_range["Fitness"].values[:generations], color="green", label=algo_name)
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("uni-dense_range.pdf", bbox_inches="tight")
plt.close(fig)

# read results for books dataset

results_folder = "./results/"
filename_genetic_point = "books_point_100000.csv"
filename_genetic_point_scaled = "books_point_100000_Scaled.csv"
filename_baseline_point = "books_point_baseline.csv"
filename_genetic_mix = "books_mix_100000.csv"
filename_genetic_mix_scaled = "books_mix_100000_Scaled.csv"
filename_baseline_mix = "books_mix_baseline.csv"
filename_genetic_range = "books_range_100000.csv"
filename_genetic_range_scaled = "books_range_100000_Scaled.csv"
filename_baseline_range = "books_range_baseline.csv"

results_genetic_point = pd.read_csv(results_folder + filename_genetic_point, header=0)
results_genetic_point_scaled = pd.read_csv(results_folder + filename_genetic_point_scaled, header=0)
results_baseline_point = pd.read_csv(results_folder + filename_baseline_point, header=0)
results_genetic_mix = pd.read_csv(results_folder + filename_genetic_mix, header=0)
results_genetic_mix_scaled = pd.read_csv(results_folder + filename_genetic_mix_scaled, header=0)
results_baseline_mix = pd.read_csv(results_folder + filename_baseline_mix, header=0)
results_genetic_range = pd.read_csv(results_folder + filename_genetic_range, header=0)
results_genetic_range_scaled = pd.read_csv(results_folder + filename_genetic_range_scaled, header=0)
results_baseline_range = pd.read_csv(results_folder + filename_baseline_range, header=0)

# plot books with point workload

bestFitness = results_genetic_point["Fitness"].min()
generations = results_genetic_point[results_genetic_point["Fitness"] == bestFitness]["Generation"].min() + 1
median_point = results_baseline_point.median()["Duration"]

fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.axhline(median_point, color="blue", linestyle=":", label="baseline")
ax.plot(results_genetic_point["Generation"].values[:generations], results_genetic_point["Fitness"].values[:generations], color="blue", label=algo_name)
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("books_point.pdf", bbox_inches="tight")
plt.close(fig)

# plot books with point workload and scaling

bestFitness = results_genetic_point["Fitness"].min()
generations = results_genetic_point[results_genetic_point["Fitness"] == bestFitness]["Generation"].min() + 1
grouped_dataframe = results_genetic_point_scaled.groupby("Dataset")

marker_style = 0
fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.plot(results_genetic_point["Generation"].values[:generations], results_genetic_point["Fitness"].values[:generations] / results_genetic_point["Fitness"].values[0], color="blue", linestyle=marker[marker_style % len(marker)], label=size_to_string(original_dataset_size))
marker_style += 1
for name, group in grouped_dataframe:
    size = group["Size"].unique()[0]
    ax.plot(group["Generation"].values[:generations], group["Fitness"].values[:generations] / group["Fitness"].values[0], color="blue", linestyle=marker[marker_style % len(marker)], label=size_to_string(size))
    marker_style += 1
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("books_point_scaled.pdf", bbox_inches="tight")
plt.close(fig)

# plot books with mix workload

bestFitness = results_genetic_mix["Fitness"].min()
generations = results_genetic_mix[results_genetic_mix["Fitness"] == bestFitness]["Generation"].min() + 1
median_mix = results_baseline_mix.median()["Duration"]

fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.axhline(median_mix, color="red", linestyle=":", label="baseline")
ax.plot(results_genetic_mix["Generation"].values[:generations], results_genetic_mix["Fitness"].values[:generations], color="red", label=algo_name)
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("books_mix.pdf", bbox_inches="tight")
plt.close(fig)

# plot books with mix workload and scaling

bestFitness = results_genetic_mix["Fitness"].min()
generations = results_genetic_mix[results_genetic_mix["Fitness"] == bestFitness]["Generation"].min() + 1
grouped_dataframe = results_genetic_mix_scaled.groupby("Dataset")

marker_style = 0
fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.plot(results_genetic_mix["Generation"].values[:generations], results_genetic_mix["Fitness"].values[:generations] / results_genetic_mix["Fitness"].values[0], color="red", linestyle=marker[marker_style % len(marker)], label=size_to_string(original_dataset_size))
marker_style += 1
for name, group in grouped_dataframe:
    size = group["Size"].unique()[0]
    ax.plot(group["Generation"].values[:generations], group["Fitness"].values[:generations] / group["Fitness"].values[0], color="red", linestyle=marker[marker_style % len(marker)], label=size_to_string(size))
    marker_style += 1
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("books_mix_scaled.pdf", bbox_inches="tight")
plt.close(fig)

# plot books with range workload

bestFitness = results_genetic_range["Fitness"].min()
generations = results_genetic_range[results_genetic_range["Fitness"] == bestFitness]["Generation"].min() + 1
median_range = results_baseline_range.median()["Duration"]

fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.axhline(median_range, color="green", linestyle=":", label="baseline")
ax.plot(results_genetic_range["Generation"].values[:generations], results_genetic_range["Fitness"].values[:generations], color="green", label=algo_name)
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("books_range.pdf", bbox_inches="tight")
plt.close(fig)

# plot books with range workload and scaling

bestFitness = results_genetic_range["Fitness"].min()
generations = results_genetic_range[results_genetic_range["Fitness"] == bestFitness]["Generation"].min() + 1
grouped_dataframe = results_genetic_range_scaled.groupby("Dataset")

marker_style = 0
fig, ax = plt.subplots(figsize=figsize)
ax.set_xlabel("generation", fontsize=fontsize_large)
ax.set_ylabel("median runtime [ms]", fontsize=fontsize_large)
ax.tick_params(axis='x', labelsize=fontsize_small)
ax.tick_params(axis='y', labelsize=fontsize_small)
ax.plot(results_genetic_range["Generation"].values[:generations], results_genetic_range["Fitness"].values[:generations] / results_genetic_range["Fitness"].values[0], color="green", linestyle=marker[marker_style % len(marker)], label=size_to_string(original_dataset_size))
marker_style += 1
for name, group in grouped_dataframe:
    size = group["Size"].unique()[0]
    ax.plot(group["Generation"].values[:generations], group["Fitness"].values[:generations] / group["Fitness"].values[0], color="green", linestyle=marker[marker_style % len(marker)], label=size_to_string(size))
    marker_style += 1
ax.legend(prop={'size': fontsize_small})
ax.set_ylim(bottom=0)
plt.savefig("books_range_scaled.pdf", bbox_inches="tight")
plt.close(fig)