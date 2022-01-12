

#!/usr/bin/python
import pandas as pd
import numpy as np
import sys 
import argparse
import os
import subprocess


# Print out a message when the program is initiated.
print('----------------------------------------------------------------\n')
print('             Add alignment informations to cluster  .\n')
print('----------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow to make mmseqs jobs for slurm')
parser.add_argument("-c", "--Cluster_file", help="The cluster file")
parser.add_argument("-r", "--Relaxed_file", help="The relaxed output from mmseqs2 analysis")
parser.add_argument("-o", "--out",help="The ouptut files")

args = parser.parse_args()

Cluster_file=args.Cluster_file
Relaxed_file=args.Relaxed_file
Output_file=args.out

print("Open Relaxed file...")
Relaxed=pd.read_csv(Relaxed_file,sep="\t",header=None)
Relaxed.columns=["query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","tlen"]
print("Open Cluster file...")
Clusters=pd.read_csv(Cluster_file,sep="\t")
print(Clusters)
try:
        Clusters.columns = ["0","Clustername","Names"]
except:
        Clusters=Clusters[["Clusternames","Names"]]
        Clusters.columns = ["Clustername","Names"]
#Cluster_and_Relaxed=pd.merge(Relaxed,Clusters,how='inner', left_on='query', right_on='Names').drop_duplicates(subset='query')

x = Relaxed.merge(Clusters, left_on="query", right_on="Names")
y = Relaxed.merge(Clusters, left_on="target", right_on="Names")

Cluster_and_Relaxed = x.merge(y, on=["query", "target", "Clustername"])


Cluster_and_Relaxed=Cluster_and_Relaxed[["Clustername","query","target","pident_x","alnlen_x","mismatch_x","gapopen_x","qstart_x","qend_x","tstart_x","tend_x","evalue_x","bits_x","tlen_x"]]
Cluster_and_Relaxed.columns=["Clustername","query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","tlen"]

Cluster_and_Relaxed=Cluster_and_Relaxed[["Clustername","query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","tlen"]]
Cluster_and_Relaxed.to_csv(Output_file,sep=";")

print("Cluster informations added with alignment informations\n")
print("Output file written to : ",Output_file)

print("Nb remaining clusters :", len(Cluster_and_Relaxed['Clustername'].unique()),"/", len(Clusters['Clustername'].unique()))
print("Nb remaining queries :", len(Cluster_and_Relaxed['query'].unique()),"/", len(Clusters.loc[Clusters['Names'].str.contains(":")]['Names'].unique()))
print("Nb remaining viral seq :", len(Cluster_and_Relaxed['target'].unique()),"/", len(Clusters.loc[~Clusters['Names'].str.contains(":")]['Names'].unique()))
