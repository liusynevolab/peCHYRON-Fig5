### Molecular recording of sequential cellular events into DNA ###
##  PER71 Bigram Frequency Analysis Pt. 2A -- Expanded Bigram Treatment Order Analysis

#   Written by Matt Demelo 
#   23 Mar 2023 - 04 Apr 2024

######

# Imports expanded bigrams data.frame generated from the PER71_ExpBigram.py script
# and computes ratios of bigrams and their reversed counterpart bigram for DI/ID,
# DC/CD, and IC/CI bigram pairs. Performs univariate kmeans clustering on bigram ratio 
# data, and plots bigram ratios colored by cluster for one biological replicate (k value
# determined as inflection point of elbow plots of k = 2 - 10. Saves a .csv file containing
# clustering information for each bigram comparison.
# Used for generating data in Fig 5 and Ext. Data Fig 9.


# import dependencies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math as math
import os
from matplotlib import rcParams
from sklearn.cluster import KMeans
from scipy.signal import find_peaks

pd.options.mode.chained_assignment = None  # default='warn', gets rid of chain copy warnings

#set wd, same as prior script.
#os.chdir('/Users/mattdemelo/Documents/LovelessLab/peChyronPaper_python/PER71_F_sig_dataframes/')

dateofint = '20240406' # for choosing the specific .csv that was generated from the PER71-ExpBigram.py script

# read bigram data.frames
exp_bigrams = pd.read_csv(dateofint+'_expanded_per71_bigrams.csv',header = 0, index_col=0)

# split bigram into position 1 and 2
exp_bigrams['Bigram'] = exp_bigrams['Bigram'].str.upper()
exp_bigrams[["P1", "P2"]] = exp_bigrams.Bigram.str.extract(r"(.)(.)")

# set font to arial globally
rcParams['font.family'] = 'Arial'


## XY/YX Bigram orientation comparison barplots

# pivot table with bigram combos for chemical epochs
chem = exp_bigrams
chemTable = (chem.pivot_table(index = 'SampleID',columns ='Bigram',values ='Count'))

chemTable = pd.DataFrame(chemTable)


# df of 'comparisons' of interest with XY/YX combos
d = {'SampleID':chemTable.index,'DI/ID':chemTable['DI']/chemTable['ID'],
     'DC/CD':chemTable['DC']/chemTable['CD'],'IC/CI':chemTable['IC']/chemTable['CI']}

chemCompare = pd.DataFrame(data=d).reset_index(drop=True)


# melted down comparison table
chemMelt = chemCompare.melt(id_vars='SampleID', var_name='Comparison', value_name='Value')

# log2 transformation of fraction
chemMelt['Log2Countcomp'] = np.log2(chemMelt['Value'])

#convert SampleID to int
chemMelt['SampleID'] = chemMelt['SampleID'].astype(int)

#remove unwanted samples
chemMelt = chemMelt[chemMelt['SampleID'] >= 27].reset_index(drop=True)

# read SampleID-Treatment key
idkey = pd.read_csv('./per71_IDs.csv',header = 0).dropna()
idkey['Trt'] = idkey['Trt'].str.replace('C', '-')

# merge Treatment into sampleIDs of expanded bigram ratio data.frame
idTable = chemMelt.merge(idkey, how = 'left',on = ['SampleID'])

# read data from total signature clustering (Part 1A) for Rep1
clustable2 = pd.read_csv('./'+dateofint+'_PER71-Clustered.csv', header =0, index_col=0).dropna().reset_index(drop=True)

#merge cluster data into 'identified' data.frame for Rep1 only

clustID = idTable.merge(clustable2, how = 'left',on = ['SampleID']).reset_index(drop=True)

# create empty final dataframe
klusterdf = pd.DataFrame()

# define sample info
reps = ['R1','R2']

sigs = ['I','D']

comparisons =  ['DI/ID','DC/CD','IC/CI']

inducers = ['IPTG','DOX']

# Define range of kmeans clusters to use
cluster_counts = list(range(2,11))

#loop over reps and signatures to produce plots
offoff = clustID[clustID['DOXclust'].str.contains('-', regex=False) & clustID['IPTGclust'].str.contains('-', regex=False)]

for r in reps:

    #### Single Replicate Clustering
    currrep = clustID[clustID['RepID'] == r]
    
    for comp in comparisons:
        
        # filter for unigram clusters of interest
        if comp == comparisons[0]:
            filtered_df = currrep[currrep['DOXclust'].str.contains('+', regex=False) & currrep['IPTGclust'].str.contains('+', regex=False)]
        elif comp == comparisons[1]:
            filtered_df = currrep[currrep['DOXclust'].str.contains('+', regex=False)]
        elif comp == comparisons[2]:
            filtered_df = currrep[currrep['IPTGclust'].str.contains('+', regex=False)]
        
        
        # filter for comparison of interest
        currcomp = filtered_df[filtered_df['Comparison'] ==comp]

            
        ## Elbow plot, DI/ID

        # Initialize an empty list to store inertia values
        inertia_values = []

        # Iterate over cluster counts
        for kl in cluster_counts:
            # Initialize KMeans model
            kmeans = KMeans(n_clusters=kl, random_state=42,n_init=10)

            # Fit KMeans model
            kmeans.fit(currcomp[['Log2Countcomp']])

            # Append the inertia to the list
            inertia_values.append(kmeans.inertia_)

        # Plot the elbow plot
        plt.figure(figsize=(8, 6))
        plt.plot(cluster_counts, inertia_values, marker='o', linestyle='-')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Inertia')
        plt.title('Elbow Plot for Univariate KMeans Clustering of '+comp+', ' + r)
        plt.xticks(cluster_counts)

        # Calculate the second derivative of the inertia values
        second_derivative = np.gradient(np.gradient(inertia_values))

        # Find peaks in the second derivative
        peaks, _ = find_peaks(second_derivative)

        # Plot the elbow plot with identified inflection point(s)

        # Highlight inflection point(s)
        for peak in peaks:
            plt.plot(cluster_counts[peak], inertia_values[peak], 'ro')  # Mark the inflection point(s) in red

        name = comp.replace('/', '-')
        plt.grid(True)
        plt.savefig(dateofint+"-"+name+"-kmeans-elbow_" + r + ".pdf",dpi=100)
        
        #### comparison barplots
        plt.clf()
        
        # Set the number of clusters to appropriate number
        kn = cluster_counts[peak]

         #Initialize KMeans model
        kmeans = KMeans(n_clusters=kn, random_state=42,n_init=10)

        # Fit KMeans model
        kmeans.fit(currcomp[['Log2Countcomp']])

        # Predict cluster labels
        currcomp[name+' Clustering'] = kmeans.labels_

        # Get mean values of Log2Countcomp for each cluster
        cluster_means = currcomp.groupby(name+' Clustering')['Log2Countcomp'].mean()

        # Sort clusters by mean values
        sorted_clusters = cluster_means.sort_values().index.tolist()

       # set appropriate cluster labels
        if comp == 'DI/ID':
            if kn == 3:
                clustlabs = ['I First', 'Middle', 'D First']
            elif kn == 4:
                clustlabs = ['I First','I/Mid', 'D/Mid', 'D First']
        else:
            if kn == 3:
                clustlabs = ['Late', 'Middle', 'Early']
            elif kn == 4:
                clustlabs = ['Latest','mid-Late', 'mid-Early', 'Earliest']


        # Create a mapping dictionary to assign labels based on sorted clusters
        cluster_labels_mapping = {
            cluster: label 
            for idx, (cluster, label) in enumerate(zip(sorted_clusters, clustlabs))
        }

        # Map cluster labels to magnitude labels
        currcomp[name+' Cluster'] = currcomp[name+' Clustering'].map(cluster_labels_mapping)
        
        
        # create plot
        
        ax = sns.barplot(data=currcomp, x='Trt', y='Log2Countcomp', color='#47889D', zorder=1, 
                         order = currcomp.sort_values('Log2Countcomp')['Trt'], palette='mako',
                         hue=name+' Cluster',hue_order=clustlabs,dodge = False, edgecolor = 'black',
                        )
        sns.despine(top=True, right=True)



        # Set custom axis labels
        ax.set_xlabel('Condition', fontsize=20)
        ax.set_ylabel((r'Log$_2$(' +comp+')') , fontsize=20)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        plt.xticks(fontsize=17)  # Adjust the font size of x-axis tick labels
        plt.yticks(fontsize=20)
        plt.legend(fontsize=20)
        #plt.ylim(-2, 2)
        plt.title(comp+' Comparison, ' + r,fontsize=14)

        plt.axhline(y=0, color='black')

        # Adjust layout
        plt.tight_layout()

        plt.savefig(dateofint+"_PER71-log2bigrambarplot_Cluster" + name+"_"+r+".pdf", dpi=1200, bbox_inches='tight')
        plt.savefig(dateofint+"_PER71-log2bigrambarplot_Cluster" + name+"_"+r+".png", dpi=1200, bbox_inches='tight')

        

        # Show the plot
        
        klusterdf = pd.concat([klusterdf,currcomp])


#### Save clustered data to dataframe/csv
ids = ['SampleID']
vals = ['RepID','IPTGclust','DOXclust','DI-ID Cluster','DC-CD Cluster', 'IC-CI Cluster']

# add 'OFF-OFF' clustered samples back to dataframe
klusterdf = pd.concat([klusterdf,offoff])

clustdat4 = klusterdf[ids + vals]

# Aggregate dataframe rows
df_merged = clustdat4.groupby('SampleID').first().reset_index()
df_merged = df_merged.merge(idkey, how = 'right',on = ['SampleID'])

# Save aggregated cluster data to csv
df_merged.to_csv(dateofint+'_PER71-ClusteredDataTable.csv',encoding='utf-8')

