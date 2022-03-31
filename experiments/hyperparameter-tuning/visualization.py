import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# enable use of latex labels in matplotlib
plt.rcParams.update({
  "text.usetex": True
})

# define constants
results_folder = "./results/"
width = 10
height_factor = 0.23
filename_grid_search = "grid_search.csv"
hyperparameters = ["mutants", "population", "tournament", "initial", "quantile"]
rows_to_print = 25

# read results
df = pd.read_csv(results_folder + filename_grid_search, header=0, index_col=0)

median_df = df[hyperparameters + ["Time"]].groupby(hyperparameters).median().sort_values("Time").reset_index()
median_df["Rank"] = median_df.index.values + 1
median_df = median_df[["Rank"] + hyperparameters + ["Time"]]
median_df["Time"] = median_df["Time"] / 1000
median_df["tournament"] = median_df["tournament"] * 100
median_df.update(median_df[['Rank', 'mutants', 'population', 'initial']].astype(int).applymap('{:d}'.format))
median_df.update(median_df[['tournament']].astype(int).applymap('{:d}\%'.format))
median_df.update(median_df[['quantile', 'Time']].astype(float).applymap('{:.2f}'.format))
median_df.rename(columns={"Rank": r"\textbf{Rank}", "mutants": r"\textbf{mutants (}$\mathbf{s_{max}}$\textbf{)}", "population": r"\textbf{population (}$\mathbf{s_\Pi}$\textbf{)}", "tournament": r"\textbf{tournament (}$\mathbf{s_T}$\textbf{)}", "initial": r"\textbf{initial (}$\mathbf{s_{init}}$\textbf{)}", "quantile": r"\textbf{quantile (}$\mathbf{q}$\textbf{)}", "Time": r"\textbf{Median Time [s]}"}, inplace=True)

mean_df = df[hyperparameters + ["Time"]].groupby(hyperparameters).mean().sort_values("Time").reset_index()
mean_df["Rank"] = mean_df.index.values + 1
mean_df = mean_df[["Rank"] + hyperparameters + ["Time"]]
mean_df["Time"] = mean_df["Time"] / 1000
mean_df["tournament"] = mean_df["tournament"] * 100
mean_df.update(mean_df[['Rank', 'mutants', 'population', 'initial']].astype(int).applymap('{:d}'.format))
mean_df.update(mean_df[['tournament']].astype(int).applymap('{:d}\%'.format))
mean_df.update(mean_df[['quantile', 'Time']].astype(float).applymap('{:.2f}'.format))
mean_df.rename(columns={"Rank": r"\textbf{Rank}", "mutants": r"\textbf{mutants (}$\mathbf{s_{max}}$\textbf{)}", "population": r"\textbf{population (}$\mathbf{s_\Pi}$\textbf{)}", "tournament": r"\textbf{tournament (}$\mathbf{s_T}$\textbf{)}", "initial": r"\textbf{initial (}$\mathbf{s_{init}}$\textbf{)}", "quantile": r"\textbf{quantile (}$\mathbf{q}$\textbf{)}", "Time": r"\textbf{Average Time [s]}"}, inplace=True)

# plot top configurations based on median runtime
figsize = (width, math.ceil(height_factor * rows_to_print))
fig, ax =plt.subplots(figsize=figsize)
ax.axis('tight')
ax.axis('off')
ax.set_title("Top " + str(min(rows_to_print, median_df.shape[0])) + " Configurations by Median Time")
table = ax.table(cellText=median_df.values[:rows_to_print, :], colLabels=median_df.columns, loc='center', cellLoc='right')
plt.savefig("best_configs_median_time.pdf", bbox_inches='tight')
plt.close()

# plot all configurations based on median runtime
figsize = (width, math.ceil(height_factor * median_df.shape[0]))
fig, ax =plt.subplots(figsize=figsize)
ax.axis('tight')
ax.axis('off')
ax.set_title("All Configurations by Median Time")
table = ax.table(cellText=median_df.values, colLabels=median_df.columns, loc='center', cellLoc='right')
plt.savefig("all_configs_median_time.pdf", bbox_inches='tight')
plt.close()

# plot top configurations based on average runtime
figsize = (width, math.ceil(height_factor * rows_to_print))
fig, ax =plt.subplots(figsize=figsize)
ax.axis('tight')
ax.axis('off')
ax.set_title("Top " + str(min(rows_to_print, mean_df.shape[0])) + " Configurations by Average Time")
table = ax.table(cellText=mean_df.values[:rows_to_print, :], colLabels=mean_df.columns, loc='center', cellLoc='right')
plt.savefig("best_configs_average_time.pdf", bbox_inches='tight')
plt.close()

# plot all configurations based on average runtime
figsize = (width, math.ceil(height_factor * mean_df.shape[0]))
fig, ax =plt.subplots(figsize=figsize)
ax.axis('tight')
ax.axis('off')
ax.set_title("All Configurations by Average Time")
table = ax.table(cellText=mean_df.values, colLabels=mean_df.columns, loc='center', cellLoc='right')
plt.savefig("all_configs_average_time.pdf", bbox_inches='tight')
plt.close()