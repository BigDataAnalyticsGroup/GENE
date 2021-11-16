import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration parameters
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)   # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)   # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE) # fontsize of the figure title

# Read result file
df = pd.read_csv('results/optimized_vs_heuristic.out')
df['time'] = df['time'] / df['workload_size']

# Create a barplot based on data in `df`
def barplot(df, dataset, axis, xlabel="", ylabel=""):
    ax = sns.barplot(x="workload_type", y="time", hue="index_type", estimator=np.median, ci="sd", capsize=.02, data=df[df['distribution'] == dataset], ax=axis)
    ax.set_title(dataset)
    ax.set_ylim(ymin=0)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set(xticklabels=[])
    return ax

# Create figure and corresponding axes
fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, sharey=True, figsize=(8.5,3))

# UNIFORM DENSE
ax = barplot(df, "uniform_dense", ax1, "", "avg. lookup time [ns]")
ax.tick_params(axis='y')
ax.legend()

# BOOKS
ax = barplot(df, "books", ax2)
ax.tick_params(axis='y')
ax.get_legend().remove()

# OSM
ax = barplot(df, "osm", ax3)
ax.tick_params(axis='y')
ax.get_legend().remove()

plt.savefig('optimized_vs_heuristic.pdf', bbox_inches='tight', format='pdf')
plt.close()
