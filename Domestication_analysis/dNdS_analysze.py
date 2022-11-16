#!/usr/bin/python
import sys 
#First take the ete3 package modified to get the Muse and Gant model:
sys.path.insert(0, "/beegfs/data/bguinet/anaconda3/pkgs/ete3-3.1.1-pyhf5214e1_0/site-packages/")
from ete3 import EvolTree, TreeStyle , TextFace
from ete3.treeview.layouts import evol_clean_layout
from ete3 import faces
import pandas as pd
import numpy as np
import sys
import argparse
import re
import os
import subprocess
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
import os
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import statistics
from io import StringIO

# Print out a message when the program is initiated.
print('------------------------------------------------------------------------------------------------------------------------------\n')
print('#dNdS_analyzer.\n')
print('\n')
print('This script will create bash scripts in order to run dNdS analysis for cluster of genes.\n') 
print('Each analysis will fill a dataframe with \n')
print('\n')
print(' Cluster_name  Mean_dNdS   Pvalue_dNdS\n')
print(' cluster1      0.34        0.20           <------- foreground branches are not significantly different from 1.    \n') 
print(' cluster2      0.002       0.000003       <------- foreground branches are significantly different from 1.\n')
print('\n')
print('foreground branches being the hymenoptera candidate loci for a viral domestication. \n')
print('If foreground branches are significantly different from 1, then the dN/dS values are under purifying or adaptative selection. \n')
print('------------------------------------------------------------------------------------------------------------------------------\n')
#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow to make busco jobs for slurm')
parser.add_argument("-tree", "--tree",help="Gene phylogeny of the cluster in nwk format")
parser.add_argument("-aln", "--alignment",help="Codon alignment of the cluster (in fasta format)")
parser.add_argument("-o", "--output_file", help="The directory where the output will be written")
parser.add_argument("-c", "--cluster_name", help="The cluster name")
parser.add_argument("-l", "--list_species",nargs='+',help="List of species to test")
parser.add_argument("-enb", "--event_number",help ="The number of the event in the phylogenetic gene tree")
parser.add_argument("-m", "--model", help="Model mode")
args = parser.parse_args()

# Variable that stores fasta sequences
Aln=args.alignment
Tree=args.tree
Output_path=args.output_file
list_of_species_to_test=args.list_species
event_number=args.event_number
model_mode=args.model
print("list_of_species_to_test;",list_of_species_to_test)
#Aln="/beegfs/data/bguinet/these/Cluster_alignment/Codon_alignment_clusters/Cluster62448.codon.aln"
#Tree="/beegfs/data/bguinet/these/Phylogeny_cluster/Cluster62448.codon.aln.treefile"
#Output_path="/beegfs/data/bguinet/these/dNdS_analysis/dNdS_results/"

Cluster_name=args.cluster_name

print("----------------------------------------------------------------")
print("Beginning of the dNdS analysis for the cluster :"+Cluster_name)
print("----------------------------------------------------------------")
print("\n")
"""
#print("Charging dataframe..")
#tab = pd.read_csv(DataFrame,sep="\t",header=None)
#list1=["ok","Unnamed: 0","Unnamed:_0","Clustername","query",	"target",	"pident",	"alnlen",	"mismatch",	"gapopen",	"qstart",	"qend",	"tstart",	"tend",	"evalue",	"bits",	"tlen",	"Taxid","species",	"genus","family",	"no_rank",	"superkingdom",	"subfamily",	"subgenus",	"suborder",	"order",	"class",	"subphylum","phylum",	"superfamily",	"infraorder", "cohort",	"infraclass",	"subclass",	"kingdom"]
#tab.columns = list1
#tab=tab.drop(["ok"  ,"Unnamed: 0","Unnamed:_0"],axis=1)
#list_of_species_to_test=[]
#subtab = tab.loc[(tab.Clustername == Cluster_name)]
#for index, row in subtab.iterrows():
#  i= re.sub("([:()])","_",row["query"])
#  list_of_species_to_test.append(i)
#list_of_species_to_test=list(dict.fromkeys(list_of_species_to_test))
#print(list_of_species_to_test)
"""


#Open the tree and the alignment codon file of the cluster of interest
print("- The tree and the alignment codon file of the cluster of interest have been charged -")
tree = EvolTree(Tree,binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")
print(tree)
#Remove leafs from the tree if they where removed from the trimming process in the Codon alignment

List_leaf_node=[]
for leaf in tree:
	List_leaf_node.append(leaf.name)
List_node_to_keep=[]
for record in SeqIO.parse(Aln, "fasta"):
  List_node_to_keep.append(record.id)

if len(List_leaf_node)!=len(List_node_to_keep):
	tree.prune(List_node_to_keep)
	tree.write(format=1, outfile=Aln+".Pruned")
	tree=EvolTree(Aln+".Pruned",binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")
"""
#for leaf in tree:
#	if leaf.name.count('_')>3: #if it is a hymenoptera locus
#		list_of_species_to_test2.append(leaf.name)
""" 
print("\n")
print('Tree of the cluster : ', Cluster_name)
print(tree)
print("\n")
tree.workdir = Output_path

tree.link_to_alignment(Aln)

if model_mode == "fb_neut":
	print("\n")
	print("--------------------------------------")
	print("Running of the the dNdS model analysis")
	print("--------------------------------------")
	print("\n")
	tree.run_model('fb'+'.'+str(Cluster_name)+'_'+str(event_number),CodonFreq=4,estFreq=1)
	print("Model run done")
	#Load the model 
tree.link_to_evol_model('/beegfs/data/bguinet/these/dNdS_analysis2/dNdS_results/fb'+'.'+str(Cluster_name)+'_'+str(event_number)+'/out', 'fb')
for model in  tree._models:
	model=tree.get_evol_model(model)
print(model)
#Here we will get the model ouput in a dataframe :
df = pd.read_csv(StringIO(model.__str__()), skiprows=6,names=['marks','omega','node_ids','name'])
df = df.applymap(lambda x: x.split(":")[1])
#Remove white space
df=df.applymap(str.strip).rename(columns=str.strip)
print(df)
marks=[]
omega=[]
for i in list_of_species_to_test:
    print(i)
    s = df.loc[df['name'] == i,'node_ids']
    print("s:",s)
    #dNdS = df.loc[df['name'] == i,'omega']
    #print(dNdS)
    #dNdS = float(dNdS)
    #s=int(s)
    #print("dNdS of",i," equal: ", dNdS)
    #print(dNdS)
    #omega.append(dNdS)
    marks.append(s)
print("Node marked :")
for i in list_of_species_to_test:
   print("- ",i)
#dNdS_average=statistics.mean(omega)
#print("dN/dS average is :",statistics.mean(omega))

print("\n")
print("-------------------------------------------------------------------------")
print(".  Adding of node marks...                                               ")
print(".  the marks will allow to force the loci to be set to a dN/dS value =1  ")
print("-------------------------------------------------------------------------")
# mark a group of branches of interest
tree.mark_tree (marks, ['#1'])
print("\n")

print("--------------------------------------------------")
print("      Running of the the neutral model         ")
print(" all marked sequence are forced to have a dNdS =1 ")
print("--------------------------------------------------")

if model_mode == "fb_neut":
	tree.run_model('b_neut'+'.'+str(Cluster_name)+'_'+str(event_number),CodonFreq=4,estFreq=1)
	print("neutral model analysis done.")
	print("\n")
if model_mode == "free":
	print("--------------------------------------------------")
	print("      Running of the the free-ratio model            ")
	print(" all marked sequence are allowed to be free       ")
	print("--------------------------------------------------") 
	#tree.run_model('b_free'+'.'+str(Cluster_name)+'_'+str(event_number),CodonFreq=4,estFreq=1)
	print("free-ratio model analysis done.")
	print("\n")
	tree.link_to_evol_model('/beegfs/data/bguinet/these/dNdS_analysis2/dNdS_results/b_neut'+'.'+str(Cluster_name)+'_'+str(event_number)+'/out', 'b_neut')
	for model in  tree._models:
		fb_model_neut=tree.get_evol_model(model)
	#fb_model_neut = tree.get_evol_model('b_neut'+"."+Cluster_name+"_"+event_number)
	#fb_model_free = tree.get_evol_model("b_free"+"."+Cluster_name+"_"+str(event_number))
	tree.link_to_evol_model('/beegfs/data/bguinet/these/dNdS_analysis2/dNdS_results/b_free'+'.'+str(Cluster_name)+'_'+str(event_number)+'/out', 'b_free')
	for model in  tree._models:
		fb_model_free=tree.get_evol_model(model)

	#Charging the free model output
	df_free = pd.read_csv(StringIO(fb_model_free.__str__()), skiprows=6,names=['marks','omega','node_ids','name'])
	df_free = df_free.applymap(lambda x: x.split(":")[1])
	#Remove white space
	df_free=df_free.applymap(str.strip).rename(columns=str.strip)

	print("\n")
	omega=[]
	for i in list_of_species_to_test:
		s = df_free.loc[df['name'] == i,'node_ids']
		#print(s)
		dNdS = df_free.loc[df['name'] == i,'omega']
		#print(dNdS)
		dNdS = float(dNdS)
		s=int(s)
		#print("dNdS of",i," equal: ", dNdS)
		#print(dNdS)
		omega.append(dNdS)

	print("Free model output: ")
	print(df_free)
	print("\n")
	print("dN/dS average is :",omega[0]) #Take the first omega value
	print("\n")

	#Run of the likelihhod ratio test...
	print("------------------------------------------------------------------")
	print("           Running of the likelihood ratio test                   ")
	print("                                                                  ")
	print("         if pvalue 'b_free' vs 'b_neut' < 0.05:                   ")
	print(" the foreground branches are significantly different from 1 '     ")
	print("------------------------------------------------------------------")

	#pvalue=tree.get_most_likely ('b_free'+"."+Cluster_name+"_"+str(event_number), 'b_neut'+"."+Cluster_name+"_"+str(event_number))
	
	try:
    		from scipy.stats import chi2
    		chi_high = lambda x, y: 1 - chi2.cdf(x, y)
	except ImportError:
    		from .utils import chi_high
	def get_most_likely2(altn, null):
 		if hasattr(altn, 'lnL') and hasattr(null, 'lnL'):
  			if  null.lnL - altn.lnL < 0:
    				return chi_high(2 * abs(altn.lnL - null.lnL),float(altn.np - null.np))
  			else:
  				warn("\nWARNING: Likelihood of the alternative model is ""smaller than null's (%f - %f = %f)" % (null.lnL, altn.lnL, null.lnL - altn.lnL) + "\nLarge differences (> 0.1) may indicate mistaken ""assigantion of null and alternative models")
  				return 1
	pvalue=get_most_likely2(fb_model_free,fb_model_neut)
	print("pvalue estimated :", pvalue)
	print("\n")

	#Insert the value into a new dataframe 
	df = pd.DataFrame(columns=("Clustername","Event","Mean_dNdS","Pvalue_dNdS"))
	df=df.append({"Clustername":Cluster_name,"Event":str(event_number),"Mean_dNdS":omega[0],"Pvalue_dNdS":pvalue},ignore_index=True)
	print("Global result :")
	print(df)
	df.to_csv(Output_path+"dNds_"+Cluster_name+"_"+str(event_number)+".out",sep="\t",index=False)
