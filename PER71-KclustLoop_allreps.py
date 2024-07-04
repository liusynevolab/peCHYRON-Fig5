### Molecular recording of sequential cellular events into DNA ###
##  PER71 Bigram Frequency Analysis Pt. 1A -- Overall Unigram Frequency Kmeans Clustering

#   Written by Matt Demelo 
#   16 Feb 2024 - 06 Apr 2024

######

# Takes processed reads of peCHYRON insertions and computes overall
# abundance of each signature in each sample, and performs univariate 
# kmeans clustering on the total signature abundance of D or I.
# Clustering information is layered onto scatterplots of D abundance vs
# I abundance, with points colored by associated cluster. (k value 
# determined as inflection point for a k = 2 - 10 elbow plot). Clustering
# information is then saved to a .csv for use in part 2 of this script pipeline.
# Used for generating data in Ext. Data Fig 7A, 8, and 9.


#import dependencies
import pandas as pd
from datetime import date
import os
import numpy as np
import math as math
import re
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from scipy.signal import find_peaks
from collections import Counter

pd.options.mode.chained_assignment = None  # default='warn', gets rid of chain copy warnings


# set file paths
# these should exist as subdirectories or files within the working directory as follows:



#set wd
#os.chdir('/Users/mattdemelo/Documents/LovelessLab/peChyronPaper_python/PER71_F_sig_dataframes/')


# get lists of directories and sort by SampleID, removing anything that is not from PER71
pedat = [filename for filename in os.listdir('./F_sig_dataframes/') if 'PER71' in filename]


## Get all 'expanded' bigrams for a given epoch, obtaining:
# -The 'read' it came from
# -The bigram identity (as an expanded bigram)
# -The read count
# -The fraction of the read populated by the bigram


# preallocate dataframes for all bigrams from all epochs of a treatment category


totalcountsdf = pd.DataFrame() #initialize data.frame



for j in range(0,len(pedat)):

    # Read current epoch as data.frame
    sigdat = pd.read_csv(('./F_sig_dataframes/'+pedat[j]), header=0).dropna()
    sigdat = sigdat.rename(columns={'signatures (1 letter)': 'signature'})

    # Extract 'absolute' signature and convert T and L signatures to D and I signatures respectively
    sigdat['signature'] = sigdat['signature'].str.upper().str.replace('T', 'D').str.replace('L', 'I')

    # Initialize a Counter to store the total counts of each character
    total_counts = Counter()
    
    CondID = pedat[j].split('_')[0]

    

    # Iterate over each row in the DataFrame
    for index, row in sigdat.iterrows():
        # Count the occurrences of each character in the signature
        signature_counts = Counter(row['signature'])

        # Multiply the counts of each character by the value in the 'counts' column for the row
        weighted_counts = {char: count * row['counts'] for char, count in signature_counts.items()}

        # Update the total counts with the weighted counts for this row
        total_counts.update(weighted_counts)

        # Remove sequencing mistakes
        total_counts['.'] = 0

    # Calculate the total sum of counts across all characters
    total_sum = sum(total_counts.values())

    # Calculate the fraction of each character in the total sum
    fractions = {char: count / total_sum for char, count in total_counts.items()}

    # Convert the total counts and fractions to a DataFrame
    sum_counts = pd.DataFrame(list(total_counts.items()), columns=['signature', 'total_counts'])
    fractions_df = pd.DataFrame(list(fractions.items()), columns=['signature', 'PercentSig'])

    # Merge the two DataFrames on the 'character' column
    result_df = pd.merge(sum_counts, fractions_df, on='signature')
    result_df['PercentSig'] = result_df['PercentSig']*100
    result_df['SampleID'] = int(CondID.split('-')[1])
    
    totalcountsdf = pd.concat([totalcountsdf,result_df])
    

dateofint = date.today().strftime("%Y%m%d")


#totalcountsdf.loc[totalcountsdf['PercentSig'] == '.', 'PercentSig'] = 'None'

# Remove irrelevant samples
totalcountsdf = totalcountsdf[totalcountsdf['SampleID']>=27]

# Add CondID identifiers between replicates
totalcountsdf['CondID'] = ((totalcountsdf['SampleID'] - 27) % 27) + 1

# Create a new column 'RepID' to identify which replicate a given row corresponds to
totalcountsdf['RepID'] = ['R1' if epoch <= 53 else 'R2' for epoch in totalcountsdf['SampleID']]

# Sort the DataFrame by 'SampleID' after adding the new columns
totalcountsdf = totalcountsdf.sort_values(by=['CondID', 'RepID']).reset_index(drop=True)

# Pivot the DataFrame to create separate columns for each unique value of 'signature'
melted_df = totalcountsdf.pivot_table(index=['SampleID', 'CondID', 'RepID'], 
                                      columns='signature', values='PercentSig', 
                                      aggfunc='first').reset_index()
print(melted_df)
# Reset column names after pivot
melted_df.columns.name = None

# If there are missing values, replace them with 0
melted_df = melted_df.fillna(0)  
melted_df = melted_df.reset_index(drop=True)

# read SampleID-Treatment key
idkey = pd.read_csv('per71_IDs.csv',header = 0).dropna()
idkey['Trt'] = idkey['Trt'].str.replace('C', '-')

# merge Treatment into sampleIDs of expanded bigram ratio data.frame
melted_df = melted_df.merge(idkey, how = 'left',on = ['SampleID'])
#melted_df = melted_df.set_index(['SampleID', 'CondID', 'RepID','Trt'])

print(melted_df)

# create empty final dataframe
klusterdf = pd.DataFrame()

# define sample info
reps = ['R1','R2']

sigs = ['I','D']

inducers = ['IPTG','DOX']

# Define range of kmeans clusters to use
cluster_counts = list(range(2,11))

#loop over reps and signatures to produce plots
for r in reps:

    #### Single Replicate Clustering
    currrep = melted_df[melted_df['RepID'] == r]
    


    
    for sig in sigs:
        # Elbow plot for k identification
        # Initialize an empty list to store inertia values

        curind = [inducer for inducer in inducers if sig in inducer][0]

        inertia_values = []

        # Iterate over cluster counts
        for kl in cluster_counts:
            # Initialize KMeans model
            kmeans = KMeans(n_clusters=kl, random_state=42, n_init=10)

            # Fit KMeans model
            kmeans.fit(currrep[[sig]])

            # Append the inertia to the list
            inertia_values.append(kmeans.inertia_)

        # Plot the elbow plot
        plt.figure(figsize=(8, 6))
        plt.plot(cluster_counts, inertia_values, marker='o', linestyle='-')
        plt.xlabel('Number of Clusters')
        plt.ylabel('Inertia')
        plt.title('Elbow Plot for Univariate KMeans Clustering of '+curind+', ' + r)
        plt.xticks(cluster_counts)

        # Calculate the second derivative of the inertia values
        second_derivative = np.gradient(np.gradient(inertia_values))

        # Find the index of the peak with the largest magnitude in the second derivative
        largest_peak_index = np.argmax(second_derivative)

        # Set the number of clusters to appropriate number
        kn = cluster_counts[largest_peak_index]

        # Plot the elbow plot with the identified inflection point
        plt.plot(cluster_counts, inertia_values, marker='o', linestyle='-')
        plt.plot(cluster_counts[largest_peak_index], inertia_values[largest_peak_index], 'ro')  # Mark the inflection point in red

        plt.grid(True)
        plt.savefig(dateofint+"-"+curind+"-kmeans-elbow_" + r + ".pdf",dpi=100)


        ## kmeans clustering, D vs I unigram plot

        # Initialize KMeans model
        kmeans = KMeans(n_clusters=kn, random_state=42,n_init=10)

        # Fit KMeans model
        kmeans.fit(currrep[[sig]])

        # Predict cluster labels
        currrep[curind+'clust'] = kmeans.labels_

        # Get mean values of Log2Countcomp for each cluster
        cluster_means = currrep.groupby(curind+'clust')[sig].mean()

        # Sort clusters by mean values
        sorted_clusters = cluster_means.sort_values().index.tolist()
        
        
        # set appropriate cluster labels
        if kn == 3:
            clustlabs = ['-', '+', '++']
        elif kn == 4:
            clustlabs = ['- -','-', '+', '++']


        # Create a mapping dictionary to assign labels based on sorted clusters
        cluster_labels_mapping = {
            cluster: label 
            for idx, (cluster, label) in enumerate(zip(sorted_clusters, clustlabs))
        }


        # Map cluster labels to magnitude labels
        currrep[curind+'clust'] = currrep[curind+'clust'].map(cluster_labels_mapping)

        # Scatterplot

        plt.figure(figsize=(8, 7))
        sns.scatterplot(data=currrep, x= 'I', y= 'D', hue=(curind+'clust'), palette='mako',hue_order=clustlabs, s= 75)
        sns.despine(top=True, right=True)
        plt.xlabel('Proportion of Total Signatures, IPTG (%)',fontsize=20)
        plt.ylabel('Proportion of Total Signatures, DOX (%)',fontsize=20)
        plt.xticks(fontsize=20)  # Adjust the font size of x-axis tick labels
        plt.yticks(fontsize=20)
        plt.legend(fontsize=20)
        plt.ylim(0, 30)
        plt.xlim(0, 50)
        plt.title(f'{kn} Clusters, '+curind+' only, ' + r,fontsize=20)


        # Save the plot as a PDF
        plt.savefig(dateofint + "-"+curind+"-kmeans_" + r + ".pdf", dpi=1200)
        plt.savefig(dateofint + "-"+curind+"-kmeans_" + r + ".png", dpi=1200)

        # create dox plot
        
        ax = sns.barplot(data=currrep, x='Trt', y='D', color='#47889D', zorder=1, 
                         order = currrep.sort_values('D')['Trt'], palette='mako',
                         hue=curind+'clust',hue_order=clustlabs,dodge = False, edgecolor = 'black',
                        )
        sns.despine(top=True, right=True)

        # Set custom axis labels
        ax.set_xlabel('Condition', fontsize=20)
        ax.set_ylabel('Proportion of Total Signatures, dox (%)' , fontsize=20)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        plt.xticks(fontsize=17)  # Adjust the font size of x-axis tick labels
        plt.yticks(fontsize=20)
        plt.legend(fontsize=20)
        #plt.ylim(-2, 2)
        #plt.title(comp+' Comparison, ' + r,fontsize=14)

        plt.axhline(y=0, color='black')

        # Adjust layout
        plt.tight_layout()

        plt.savefig(dateofint+"_PER71-doxplot_Cluster" + r + ".pdf", dpi=1200, bbox_inches='tight')
        plt.savefig(dateofint+"_PER71-doxplot_Cluster" +r+".png", dpi=1200, bbox_inches='tight')

    # append replicate data to larger data.frame
    klusterdf = pd.concat([klusterdf,currrep])
   
# save klusterdf to .csv file
klusterdf.to_csv(dateofint+'_PER71-Clustered.csv',encoding='utf-8')