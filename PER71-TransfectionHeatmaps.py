### Molecular recording of sequential cellular events into DNA ###
##  PER71 Transfection and Chemical Epochs, Bigram Frequency Pt. 2B
##  Bigram Heatmap and orientation fraction barplots

#   Written by Matt Demelo 
#   04 Apr 2024

######

# Imports expanded bigrams data.frame generated from the PER71_ExpBigram.py script
# for only transfection epochs (PER58 and PER60), and plots a heatmap of Position 1
# vs Position 2 of every possible bigram, with 'heat' determined by bigram count.
# # Used for generating data in Ext. Data Fig 7A.


# import dependencies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math as math
import os
from matplotlib import rcParams

pd.options.mode.chained_assignment = None  # default='warn', gets rid of chain copy warnings

#set wd
os.chdir('/Users/mattdemelo/Documents/LovelessLab/peChyronPaper_python/PER71_F_sig_dataframes/')

# get lists of directories and sort by epoch num


dateofint = '20240404' # for choosing the specific .csv that was generated from the ExpBigram_AllEpoch.py script

# read bigram data.frames
exp_bigrams = pd.read_csv(dateofint+'_expanded_TransEpoch-PER58-60_bigrams.csv',header = 0, index_col=0)

# split bigram into position 1 and 2
exp_bigrams['Bigram'] = exp_bigrams['Bigram'].str.upper()
exp_bigrams[["P1", "P2"]] = exp_bigrams.Bigram.str.extract(r"(.)(.)")



## Draw heatmap figure for epochs.
#  Primarily for transfection epochs, can also do chemical

# Set font family
plt.rcParams['font.family'] = 'Arial'
rcParams['xtick.labelsize'] = 28  # Size of x-axis tick labels
rcParams['ytick.labelsize'] = 28  # Size of y-axis tick labels
rcParams['axes.labelsize'] = 28   # Size of axis labels
rcParams['axes.titlesize'] = 48   # Size of plot titles

# Loop over the facets
for title, sample_data in zip(['RPUI', 'URSP'], exp_bigrams.groupby('SampleID')):
    sample_id, data = sample_data
    plt.figure(figsize=(8, 8))  # Set the size of the figure to be square
    d = data.pivot_table(index='P1', columns='P2', values='Count', aggfunc=np.sum).fillna(0)  # Pivot the data
    sns.heatmap(d, annot=False, cmap="viridis", vmin=0, vmax=20000).set_facecolor('black')
    plt.title(title)  # Set plot title
    plt.xlabel('P2')  # Set x-axis label
    plt.ylabel('P1')  # Set y-axis label
    plt.tight_layout()  # Adjust layout to prevent overlap
    plt.savefig(f"{dateofint}_{title}_transepoch_bigrams.pdf", dpi=1200)  # Save the plot
    plt.savefig(f"{dateofint}_{title}_transepoch_bigrams.png", dpi=1200)  # Save the plot


