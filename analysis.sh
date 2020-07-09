#!/bin/bash
# To call the two python scripts Jaccard_index and Cluster_analysis and implements the algorithms with the python scripts.

#Contains Jaccard Index algorithm that reads, cleanse and scores the chemical compounds and genes and writes the output to a csv file
python3 Jaccard_Index.py

#Clustering analysis algorithm that reads the score files, data merge,matrix conversion and hierarchical clustering.
python3 Cluster_analysis.py

