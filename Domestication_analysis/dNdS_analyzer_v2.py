#!/usr/bin/python
from warnings import warn
import re
import shutil
from tempfile import mkstemp
import sys 
#First take the ete3 package modified to get the Muse and Gant model:
sys.path.insert(0, "/beegfs/data/bguinet/anaconda3/pkgs/ete3-3.1.1-pyhf5214e1_0/site-packages/")
from ete3 import EvolTree
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
import statistics
from io import StringIO
import random 
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
#Output_path="/beegfs/data/bguinet/these/dNdS_analysis/dNdS_results2/"

Cluster_name=args.cluster_name

print("----------------------------------------------------------------")
print("Beginning of the dNdS analysis for the cluster :"+Cluster_name)
print("----------------------------------------------------------------")
print("\n")

def sed(pattern, replace, source, dest=None, count=0):
    """Reads a source file and writes the destination file.
    In each line, replaces pattern with replace.
    Args:
        pattern (str): pattern to match (can be re.pattern)
en the tree and the alignment codon file of the cluster of interest
print("- The tree and the alignment codon file of the cluster of interest have been charged -")

try:
        tree = EvolTree(Tree,binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")
except:
        tree = EvolTree("/beegfs/data/bguinet/these/Cluster_phylogeny_filtred/"+Cluster_name+"_AA.dna.treefile",binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")



        replace (str): replacement str
        source  (str): input filename
        count (int): number of occurrences to replace
        dest (str):   destination filename, if not given, source will be over written.        
    """
    fin = open(source, 'r')
    num_replaced = count
    if dest:
        fout = open(dest, 'w')
    else:
        fd, name = mkstemp()
        fout = open(name, 'w')
    for line in fin:
        out = re.sub(pattern, replace, line)
        fout.write(out)
        if out != line:
            num_replaced += 1
        if count and num_replaced > count:
            break
    try:
        fout.writelines(fin.readlines())
    except Exception as E:
        raise E
    fin.close()
    fout.close()
    if not dest:
        shutil.move(name, source)

#Open the tree and the alignment codon file of the cluster of interest
print("- The tree and the alignment codon file of the cluster of interest have been charged -")

try:
	tree = EvolTree(Tree,binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")
except:
	tree = EvolTree("/beegfs/data/bguinet/these/Cluster_phylogeny_filtred/"+Cluster_name+"_AA.dna.treefile",binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")

#Reduce complexity of the tree to speed up the analysis 

if os.path.exists("/beegfs/data/bguinet/these/Cluster_alignment_filtred/Codon_alignment_clusters/"+Cluster_name+"_NT.dna_without_shift_cleaned2"):
	record_dict = SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Cluster_alignment_filtred/Codon_alignment_clusters/"+Cluster_name+"_NT.dna_without_shift_cleaned2", "fasta"))
elif os.path.exists("/beegfs/data/bguinet/these/Cluster_alignment_filtred/Codon_alignment_clusters/"+Cluster_name+"_NT.dna_without_shift_cleaned"):
	record_dict = SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Cluster_alignment_filtred/Codon_alignment_clusters/"+Cluster_name+"_NT.dna_without_shift_cleaned", "fasta"))

#record_dict=SeqIO.to_dict(SeqIO.parse(Aln,"fasta"))
list_records=[]
list_virus=[]
for i in record_dict:
  list_records.append(record_dict[i].id)
  if '__' not in record_dict[i].id:
   list_virus.append(record_dict[i].id)
print(list_records)

list_records=list({s.rsplit('__', 1)[-1]: s for s in list_records[::-1]}.values())[::-1]

focal_loci=list_of_species_to_test.copy()
list_records=list(set(list_records + focal_loci))

#Only temporary  code 
list_records=[]
for leaf in tree:
	list_records.append(leaf.name)


#If nb viral seq >30: (to speed up the analyse) 
if len(list_virus)>30:
 print("Need to remove some viral sequence to speed up the analyse...")
 species_to_keep= list_of_species_to_test.copy()
 tab=pd.read_csv("/beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono.m8",sep=";")
 subtab=tab.loc[tab['Clustername'].str.contains(Cluster_name) & tab['Event'].eq(event_number)]
 list_all_viruses=tab.loc[tab['Clustername'].str.contains(Cluster_name)]['target'].unique()
 subset_target = random.sample(list(list_all_viruses) , int(len(list_virus)/4))
 for i in subtab['target'].unique():
  if i in list_records:
   species_to_keep.append(i)
 for i in subset_target:
  if i in list_records:
   species_to_keep.append(i)
 tree.prune(species_to_keep)
 list_records=[]
 new_list_virus=[]
 for leaf in tree:
   list_records.append(leaf.name)
   if '__' not in lef:
    new_list_virus.append(leaf)
 print(" Nb viral sequence remove : ",len(list_virus)- len(new_list_virus))

#Remove fasta sequences according to the tree:

with open("/beegfs/data/bguinet/these/Cluster_alignment_filtred/Codon_alignment_clusters/"+Cluster_name+"_NT.dna_without_shift_cleaned.Pruned","w") as outputfile:
  for i in record_dict:
    if record_dict[i].id in list_records:
      print('>',record_dict[i].id)
      print('>',record_dict[i].id,sep="",file=outputfile)
      print(record_dict[i].seq,file=outputfile)


Aln="/beegfs/data/bguinet/these/Cluster_alignment_filtred/Codon_alignment_clusters/"+Cluster_name+"_NT.dna_without_shift_cleaned.Pruned"

new_tree=tree.write (format=0,outfile=Tree+".pruned")
tree = EvolTree(new_tree,binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")

print(tree)
print("\n")

#Remove leafs from the tree if they where removed from the trimming process in the Codon alignment

List_leaf_node=[]
for leaf in tree:
	List_leaf_node.append(leaf.name)
List_node_to_keep=[]
for record in SeqIO.parse(Aln, "fasta"):
  List_node_to_keep.append(record.id)
	
if len(List_leaf_node)!=len(List_node_to_keep):
	print("étree pruned at: ",Aln,".Pruned")
	tree.prune(List_node_to_keep)
	tree.write(format=1, outfile=Aln+".Pruned")
	tree=EvolTree(Aln+".Pruned",binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")
	print(tree)
print("\n")
print('Tree of the cluster : ', Cluster_name)


#Get nodeid 
df_marks = pd.DataFrame(columns=['Node_name', 'Node id'])
ancestor=tree.get_common_ancestor(list_of_species_to_test)
print(list_of_species_to_test)
print(ancestor)
print("Type4")
print(type(ancestor))
for node in ancestor.traverse("levelorder"):	
	try:
		df_marks=df_marks.append({"Node_name":node.name,"Node_id":node.node_id}, ignore_index=True)
	except:
		continue
print(df_marks)
#In order to deal with small phylogeny
print(df_marks.shape[0]) 
if df_marks.shape[0] < 5:
	df_marks=df_marks.loc[df_marks["Node_name"].str.contains("__")]
else:
	df_marks = df_marks.iloc[1:]

print(df_marks)
marks=list(df_marks['Node_id'].astype(int))

tree.workdir = Output_path
tree.link_to_alignment(Aln)


print("\n")
print("-------------------------------------------------------------------------")
print(".  Adding of node marks...                                               ")
print(".  the marks will allow to force the loci to be set to a dN/dS value =1  ")
print("-------------------------------------------------------------------------")

# mark a group of branches of interest
tree.mark_tree(marks, ['#1'])
print('node marked :') 
print(tree.write ())
print("\n")

print("--------------------------------------------------")
print("      Running of the the neutral model         ")
print(" all marked sequence are forced to have a dNdS =1 ")
print("--------------------------------------------------")
import shutil
#remove previous run
if os.path.isdir('/beegfs/data/bguinet/these/dNdS_analysis_filtred/dNdS_results/b_neut'+'.'+str(Cluster_name)+'_'+str(event_number)):
 shutil.rmtree('/beegfs/data/bguinet/these/dNdS_analysis_filtred/dNdS_results/b_neut'+'.'+str(Cluster_name)+'_'+str(event_number)+'/')

if os.path.isdir('/beegfs/data/bguinet/these/dNdS_analysis_filtred/dNdS_results/b_free'+'.'+str(Cluster_name)+'_'+str(event_number)):
 shutil.rmtree('/beegfs/data/bguinet/these/dNdS_analysis_filtred/dNdS_results/b_free'+'.'+str(Cluster_name)+'_'+str(event_number)+'/')

if model_mode == "fb_neut" or model_mode=="ALL":
	tree.run_model('b_neut'+'.'+str(Cluster_name)+'_'+str(event_number),CodonFreq=4,estFreq=1,getSE=1,noisy=3,verbose=1)
	print("neutral model analysis done.")
	print("\n")
if model_mode == "free" or model_mode == "ALL":
	print("--------------------------------------------------")
	print("      Running of the the free-ratio model            ")
	print(" all marked sequence are allowed to be free       ")
	print("--------------------------------------------------") 
	tree.run_model('b_free'+'.'+str(Cluster_name)+'_'+str(event_number),CodonFreq=4,estFreq=1,getSE=1,noisy=3,verbose=1)
	print("free-ratio model analysis done.")
	print("\n")
	tree.link_to_evol_model('/beegfs/data/bguinet/these/dNdS_analysis_filtred/dNdS_results/b_neut'+'.'+str(Cluster_name)+'_'+str(event_number)+'/out', 'b_neut')
	for model in  tree._models:
		fb_model_neut=tree.get_evol_model(model)
	tree.link_to_evol_model('/beegfs/data/bguinet/these/dNdS_analysis_filtred/dNdS_results/b_free'+'.'+str(Cluster_name)+'_'+str(event_number)+'/out', 'b_free')
	for model in  tree._models:
		fb_model_free=tree.get_evol_model(model)
	#Charging the free model output
	df_free = pd.read_csv(StringIO(fb_model_free.__str__()), skiprows=6,names=['marks','omega','node_ids','name'])
	df_free = df_free.applymap(lambda x: x.split(":")[1])
	#Remove white space
	df_free=df_free.applymap(str.strip).rename(columns=str.strip)
	print("\n")
	print(df_free)
	omega=[]
	for i in list_of_species_to_test:
		s = df_free.loc[df_free['name'] == i,'node_ids']
		#print(s)
		dNdS = df_free.loc[df_free['name'] == i,'omega']
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

	#get the Standar Error of the dN/dS estimation : 
	with open('/beegfs/data/bguinet/these/dNdS_analysis_filtred/dNdS_results/b_free'+'.'+str(Cluster_name)+'_'+str(event_number)+'/out', "r") as ifile:
		for line in ifile:
			if line.startswith("SEs for parameters:"):
				SE=next(ifile, ' ').strip()
				SE=re.split('\s+', SE)
				SE=SE[-1]
	#Insert the value into a new dataframe 
	df = pd.DataFrame(columns=("Clustername","Event","Mean_dNdS","Pvalue_dNdS"))
	df=df.append({"Clustername":Cluster_name,"Event":str(event_number),"Mean_dNdS":omega[0],"Pvalue_dNdS":pvalue,"SE_dNdS":SE},ignore_index=True)
	print("Global result :")
	print(df)
	df.to_csv(Output_path+"dNds_"+Cluster_name+"_"+str(event_number)+".out",sep="\t",index=False)
