### Molecular recording of sequential cellular events into DNA ###
##  PER71 Bigram Frequency Analysis Pt. 1B -- Expanded bigram identification

#   Written by Matt Demelo 
#   23 Mar 2023 - 30 Mar 2024

######

# Takes processed reads of peCHYRON insertions and identifies
# 'expanded bigrams' within each insertion and their frequency
# (corresponding to the number of reads for that insertion) for
# each all transfection and chemical epochs. See methods section
# for an in-depth definition of the 'expanded bigram' identification.
# Used for generating data in Ext. Data Fig 8 and 9.


#import dependencies
import pandas as pd
from datetime import date
import os
import numpy as np
import math as math



#set wd
os.chdir('/Users/tl115/Dropbox/peCHYRON paper 2/NatureComms/Figure_Panels/Figure 5/')

# get lists of directories and sort by epoch num
pedat = os.listdir('./F_sig_dataframes/')


# Get all 'expanded' bigrams for a given epoch, obtaining:
# -The 'read' it came from
# -The bigram identity (as an expanded bigram)
# -The read count
# -The fraction of the read populated by the bigram


# Function for 'reducing' homopolymeric signature stretches within reads to one signature
def bigram_reducer(string):
    reduced_string = string[0]  # Initialize the reduced string with the first character
    
    for i in range(1, len(string)):
        if string[i] != string[i-1]:
            reduced_string += string[i]  # Append non-repeated characters to the reduced string
    
    return reduced_string

# preallocate dataframes for all bigrams from all epochs of a treatment category
b = {'Expt':[],'SampleID':[], 'Read': [], 'Bigram': [], 'Count':[], 'BigramFraction': []}

exp_bigrams = pd.DataFrame(data = b) #initialize data.frame




for j in range(0,len(pedat)):

    # read current epoch as data.frame
    seqreaddat = pd.read_csv((pedat[j]),header = 0).dropna()

    #convert A->B and B->A signatures to common signature
    seqreaddat['signatures (1 letter)'] = seqreaddat['signatures (1 letter)'].str.upper()

    # remove errant assignments

    seqreaddat = seqreaddat[-seqreaddat['signatures (1 letter)'].str.contains('\.')] # removes blanks

    
   # seqreaddat['signatures (1 letter)'] = seqreaddat['signatures (1 letter)'].str.replace('C','') #removes constitutive
                                                                                                      #from transfection epochs

    seqreaddat = seqreaddat.loc[seqreaddat['signatures (1 letter)'] != ''] # removes rows that previously only empty read rows

    seqreaddat = seqreaddat.reset_index(drop=True) #reset indices
    
    # Convert T and L signatures to D and I signatures respectively.
    seqreaddat['signatures (1 letter)'] = seqreaddat['signatures (1 letter)'].str.replace('T', 'D').str.replace('L', 'I')


    # Get epoch number
    
    sampleid = pedat[j].split('_')[0]
    epnum = int(sampleid.split('-')[1])
    #preallocate dataframes for bigrams from epoch
    bigrams = pd.DataFrame(data = b)

    # for an epoch, iterate over all reads in epoch
    for k in range(0,len(seqreaddat['counts'])):


        #initialize read
        expt = pedat[j][0:5]
        read =  seqreaddat.iloc[k].str.upper() #get row for a signature set
        readlength = seqreaddat['length'][k] # length of all signatures
        readcount = int(seqreaddat['counts'][k]) # count of this signature in the epoch
        readgram = seqreaddat['signatures (1 letter)'][k].upper() #get signatures
        readgram_og = readgram # save initial read with signatures before reduction
        readgram = bigram_reducer(readgram) #reduce repetitive signatures to 1


        #Get bigrams
        bigramlist = []

        # build a list of bigrams found within the whole 'read'
        # with expanded bigrams, a bigram is P1 and any possible unigram found after
        # within the read, even with separation.
        if len(readgram)>2: # if len >2, find possible bigrams
            for w in range(len(readgram)-1): # iterate through every possible P1 position
                P1 = readgram[w]
                for v in range (w+1,len(readgram)): # iterate through every position after P1
                    P2 = readgram[v] # define P2 unigram
                    bigramlist.append(P1+P2) #build bigram from P1 and P2        
        else: #if len < 2, only 1 possible.
            bigramlist = [readgram] # if read is not greater than 2, only 1 bigram is possible.


        # Store bigrams

        if(len(bigramlist)>0)&(len(''.join(set(bigramlist)))>1):

            # build a dictionary which gives each bigram found in a read the read count and accessory info
            d = {'Expt':expt,'SampleID':epnum, 'Read': [readgram_og]*len(bigramlist), 
                 'Bigram': bigramlist, 'Count':readcount*len(bigramlist), 
                 'BigramFraction': [1/len(bigramlist)]*len(bigramlist)}

            bigram_df = pd.DataFrame(data=d) # put in data.frame
            
            # append to main df
            bigrams = pd.concat([bigrams,bigram_df])
            
        
        else:
            bigram_df = pd.DataFrame(data = b) # if no bigrams found, set as empty and do not append


    #concat to all epoch bigram df
    exp_bigrams = pd.concat([exp_bigrams,bigrams]) #build all epoch bigram data.frame
    exp_bigrams = exp_bigrams.reset_index(drop=True)

    
# write data.frames to csv for safe keeping
# saves with date title for csv
dateofint = date.today().strftime("%Y%m%d")

exp_bigrams = exp_bigrams.drop_duplicates()

per71_bigrams = exp_bigrams[exp_bigrams['Expt'] == 'PER71']

per71_bigrams.to_csv(dateofint+'_expanded_per71_bigrams.csv',encoding='utf-8')

tepoch_bigrams = exp_bigrams[exp_bigrams['Expt'] != 'PER71']

tepoch_bigrams.to_csv(dateofint+'_expanded_TransEpoch-PER58-60_bigrams.csv',encoding='utf-8')
