# Guide for Bigram Frequency Analysis

This is a step by step guide for performing the bigram frequency analysis in Figure 5b, Extended Data Figures 8-9, and Supplementary Figure 6a of 'Open-ended molecular recording of sequential cellular events into DNA'. Note that this analysis was performed on a Apple Silicon (ARM-based) MacBook. 

## Setup

### Python setup
Before running the scripts in this pipeline, make sure you have installed Python 3 and have the following Python 3 Libraries (many of which are automatically installed if Python 3 is installed using [Anaconda](https://www.anaconda.com/download)):

 - pandas
 - numpy
 - matplotlib
 - seaborn
 - math
 - os
 - scikit-learn (sklearn)
 - scipy

### Files and Directories
1) Create a working directory to perform the analysis in a location of your choosing (ex. /Downloads/BigramFun)

2) Download the **F_sig_dataframes** directory, and move this to your working directory so that it is a subdirectory of your working directory (ex. /Downloads/BigramFun/**F_sig_dataframes**).

3) Download the following python files and put them in your working directory:
	- PER71_ExpBigram.py
	- PER71_KclustLoop_allreps.py
	- PER71_ExpBigram-KClustBarPlots-AllReps.py
	- PER71_TransfectionHeatmaps.py

4) Open up a terminal instance, and navigate to your working directory (customize for your filepath):
`username$ cd /Downloads/BigramFun/`

5) Check your working directory with the `ls` command, which should show the following:

`username$ ls`
`F_sig_dataframes`
`PER71_ExpBigram.py`
`PER71_KclustLoop_allreps.py`
`PER71_ExpBigram-KClustBarPlots-AllReps.py`
`PER71_TransfectionHeatmaps.py`

6) Using a text editor of choice, open each .py file and edit the working directory path to match the one above:
Change:
`os.chdir('/placeholder/filepath/change/as/needed/')`
To your working directory:
`os.chdir('/Downloads/BigramFun/')`
At these lines:
`PER71_ExpBigram.py` -- Line 27
`PER71_KclustLoop_allreps.py` -- Line 41
`PER71_ExpBigram-KClustBarPlots-AllReps.py` -- Line 32
`PER71_TransfectionHeatmaps.py` -- Line 28

7) At this point, the scripts should be ready to run.

## Running the script pipeline
1) In a terminal window, run the `PER71_KclustLoop_allreps.py`. This script performs univariate K-means clustering for each signature to define clusters of samples by unigram frequency (see script for more info):
`username$ python PER71_KclustLoop_allreps.py`

2) In your working directory verify that the following file was created:
`YYYYMMDD-inducer-kmeans-elbow_R#.pdf; YYYYMMDD-inducer-kmeans_R#.pdf; YYYYMMDD-inducer-kmeans_R#.png; YYYYMMDD_PER71-Clustered.csv`

	
3) In a terminal window, run the `PER71_ExpBigram.py`. This script creates plots for Figure 5A and identifies bigrams in each sequencing read and enumerates unique bigram combinations in each read (see Methods and script for more details):
`username$ python PER71_ExpBigram.py`

4) In your working directory verify that the following files were created
`YYYYMMDD_expanded_per71_bigrams.csv; YYYYMMDD_expanded_TransEpoch-PER58-60_bigrams.csv `

5) In a terminal window, run the `PER71_ExpBigram-KClustBarPlots-AllReps.py`. This script generates plots for Figure 5B using univariate K-means clustering performed previously to group samples into IPTG or Dox-positive groups, and then perform univariate K-means clustering on bigram-anagram comparisons to identify treatment epochs (see script and methods for more details):
`username$ python PER71_ExpBigram-KClustBarPlots-AllReps.py`

6) In your working directory verify that the plot files were created and dataframes are saved to .csv:
`YYYYMMDD_PER71-log2bigrambarplot_ClusterXY-YX_R#.pdf; YYYYMMDD_PER71-log2bigrambarplot_ClusterXY-YX_R#.png; YYYYMMDD_PER71-ClusteredDataTable.csv; YYYYMMDD_PER71-AllColumns_ClusteredDataTable.csv`

7) In a terminal window, run the `PER71_TransfectionHeatmaps.py`. This script generates plots for Extended data figure 7c, comparing the frequency of unigrams in bigrams for two given expected orders of induction (RPUI, URSP):
`username$ python PER71_TransfectionHeatmaps.py`

8) In your working directory verify that the plot files were created and dataframes are saved to .csv:
`YYYYMMDD_URSP_transepoch_bigrams.pdf; YYYYMMDD_URSP_transepoch_bigrams.png;
YYYYMMDD_RPUI_transepoch_bigrams.pdf;
YYYYMMDD_RPUI_transepoch_bigrams.png`
