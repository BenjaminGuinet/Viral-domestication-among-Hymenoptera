

import os
import argparse
import pandas as pd 
import numpy as np
import re 
from taxadb.taxid import TaxID
from taxadb.accessionid import AccessionID
from taxadb.taxid import TaxID
import pandas as pd
import re,os 
import subprocess
from Bio import SeqIO 
import numpy  as np 
import networkx as nx



#This script will try to find additional hits between filamentous ORFS from LbFV, DFV and the newly described PoFV within the genome of Platygaster orseoliae
#Load the files 
PoFV_ORFs="/beegfs/data/bguinet/LbFV_family_project/Genomes/PoFV/Predicted_orfs/Final_ORF_prediction_PoFV.faa"
LbFV_ORFs="/beegfs/data/bguinet/LbFV_family_project/Genomes/LbFV/Predicted_orfs/Previously_predicted_108_ORFs.faa"
DFV_ORFS="/beegfs/data/bguinet/LbFV_family_project/Genomes/DFV/Predicted_orfs/Previously_predicted_68_ORFs.faa"
ALL_ORFs="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/ALL_filamentous_ORFs.aa"
ALL_ORFs_db="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/ALL_filamentous_ORFs_db"


Porseoliae_assembly="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/Platygaster_orseoliae_corrected2.fa"
Porseoliae_assembly_db="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/Platygaster_orseoliae_corrected2_db"

ALL_filamentous_ORFs="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/ALL_filamentous_ORFs.aa" 
ALL_filamentous_ORFs_db="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/ALL_filamentous_ORFs_db" 


#MMseqs results files
ALL_ORFs_vs_Porseoliae_result="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result"
ALL_ORFs_vs_Porseoliae_temp="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_temp"
ALL_ORFs_vs_Porseoliae_table="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result.m8"


# Run first mmseqs analysis 

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Porseoliae_assembly + " "+ Porseoliae_assembly_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ ALL_filamentous_ORFs + " "+ ALL_filamentous_ORFs_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs search "+ ALL_filamentous_ORFs_db + " "+ Porseoliae_assembly_db + " "+ ALL_ORFs_vs_Porseoliae_result  + " "+ ALL_ORFs_vs_Porseoliae_temp  + " "+  -e 0.000001 --threads 20 "  , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs convertalis --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qlen,qcov,tcov' "+ ALL_filamentous_ORFs_db + " "+ Porseoliae_assembly_db + " "+ ALL_ORFs_vs_Porseoliae_result  + " "+ ALL_ORFs_vs_Porseoliae_table  , shell=True)


# Open the result file 
PoFV_vs_filamentous_ORFs=pd.read_csv(ALL_ORFs_vs_Porseoliae_table, sep="\t")
PoFV_vs_filamentous_ORFs.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']


#Filter 
PoFV_vs_filamentous_ORFs=PoFV_vs_filamentous_ORFs[PoFV_vs_filamentous_ORFs['bits'].gt(50)]
#Remove PoFV scaffolds within Porseoliae assembly 
PoFV_scaffolds=SeqIO.to_dict(SeqIO.parse(PoFV_ORFs, "fasta"))
list_PoFV_scaffolds=[]
for scaff in PoFV_scaffolds.keys():
  scaffold=re.sub("PoFV_","",scaff)
  scaffold=re.sub("_orf.*","",scaffold)
  list_PoFV_scaffolds.append(scaffold)
  
#Keep only scaffolds not present within PoFV 
PoFV_vs_filamentous_ORFs=PoFV_vs_filamentous_ORFs[~PoFV_vs_filamentous_ORFs['target'].isin(list_PoFV_scaffolds)] 

###############################################
## ADD coverage to each remaining scaffolds  ##
###############################################
Add coverage 
cov_tab=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/Mapping/cov_GC_new_mean.tab",sep="\t")
#Get mean FL coverage =
#cov_tab.loc[cov_tab['scaffolds'].isin(List_FL_scaffolds)].median()
#Median cov depth Porseoliae BUSCO =33X  (look for 40X max)
#Merge both informations 
PoFV_vs_filamentous_ORFs= PoFV_vs_filamentous_ORFs.merge(cov_tab,left_on="target",right_on="scaffolds",how="left")
BUSCO_tab=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/run_busco/run_BUSCO_v3/full_table_Platygaster_orseoliae_BUSCO_v3_newcoord.tsv",sep=";")
#Keep only complete and fragmented BUSCOs 
BUSCO_tab=BUSCO_tab.loc[BUSCO_tab['Status'].isin(["Complete","fragmented"])]
#Merge both informations 
BUSCO_tab= BUSCO_tab.merge(cov_tab,left_on="Contig",right_on="scaffolds",how="left")
#To remove unknown species scaffolds
BUSCO_tab=BUSCO_tab.loc[~BUSCO_tab['Median_cov_depth'].lt(20)]
BUSCO_tab=BUSCO_tab.loc[~BUSCO_tab['Median_cov_depth'].isna()]
from statistics import mean
list_scaff_count_cov_busco=[]
#For each cov_dept, put it into the list
for i in BUSCO_tab['Median_cov_depth']:
                        list_scaff_count_cov_busco.append(i)
len_scaff_count_cov_busco=len(list_scaff_count_cov_busco)
mean_list_scaff_count_cov_busco=mean(list_scaff_count_cov_busco)
list_cov_scaff_viral=[]
for index, row in PoFV_vs_filamentous_ORFs.iterrows():
        if row['Median_cov_depth'] > mean_list_scaff_count_cov_busco:
                pvalue=(sum(i > row['Median_cov_depth'] for i in list_scaff_count_cov_busco)/len_scaff_count_cov_busco)
        else:
                pvalue=(sum(i < row['Median_cov_depth'] for i in list_scaff_count_cov_busco)/len_scaff_count_cov_busco)
        list_cov_scaff_viral.append({'scaf_name':row['target'],'scaf_length':row['qlen'],'cov_depth_candidat':row['Median_cov_depth'],'cov_depth_BUSCO':mean_list_scaff_count_cov_busco,'pvalue_cov':pvalue})


PoFV_vs_filamentous_ORFs_cov= PoFV_vs_filamentous_ORFs.merge(pd.DataFrame(list_cov_scaff_viral),left_on="target",right_on="scaf_name",how="left")
PoFV_vs_filamentous_ORFs_cov = PoFV_vs_filamentous_ORFs_cov.drop_duplicates(subset=['query','target','tstart','tend'], keep='first')
#Only keep scaffold with typical euk coverage 
PoFV_vs_filamentous_ORFs_cov= PoFV_vs_filamentous_ORFs_cov.loc[PoFV_vs_filamentous_ORFs_cov['cov_depth_candidat'].lt(40) & PoFV_vs_filamentous_ORFs_cov['cov_depth_candidat'].gt(20)]


#########
# save1 #
#########

PoFV_vs_filamentous_ORFs_cov.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov.m8",sep=";",index=False)
PoFV_vs_filamentous_ORFs_cov= pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov.m8",sep=";")


##########################################################
#Look for HSPs to merge manually thanks to the assemblies 
##########################################################

PoFV_vs_filamentous_ORFs_cov['count_duplicate'] = PoFV_vs_filamentous_ORFs_cov.groupby(['query','target'])['query'].transform('size')
Putative_HSPs=PoFV_vs_filamentous_ORFs_cov.loc[PoFV_vs_filamentous_ORFs_cov['count_duplicate'].ge(2)]
Putative_HSPs=Putative_HSPs.loc[Putative_HSPs['query'].str.contains("PoFV")]
Putative_HSPs=Putative_HSPs.sort_values(['query', 'target'], ascending=[True, True])
Putative_HSPs['diff_length']=np.nan
Putative_HSPs['HSP_group']=np.nan
grouped = Putative_HSPs.groupby(['query', 'target'])
number=1
for group_name, group in grouped:
  diff=abs(min(group['qend'])-max(group['qstart']))
  Group=number
  for index, row in group.iterrows():
    #print(group_name)
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']),"diff_length"]= diff
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']),"HSP_group"]= Group
  number+=1

#Count number HSP 
Putative_HSPs['count_HSP']=Putative_HSPs.groupby("HSP_group")["HSP_group"].transform("count")

Putative_HSPs.loc[Putative_HSPs['count_HSP'].gt(2),"diff_length"]= 70

# If the distance between the different HSPs is > 80, do not merge them 
Putative_HSPs=Putative_HSPs.loc[Putative_HSPs['diff_length'].lt(80)]
PoFV_vs_filamentous_ORFs_cov['diff_length']=np.nan
PoFV_vs_filamentous_ORFs_cov['HSP_group']=np.nan
Putative_HSPs['new_start']=np.nan
Putative_HSPs['new_end']=np.nan
Putative_HSPs['strand']="+"
Putative_HSPs.loc[Putative_HSPs['tstart']> Putative_HSPs['tend'],"strand"]="-"
m = Putative_HSPs['strand'].eq('-')
Putative_HSPs[['tstart','tend']] = np.where(m.to_numpy()[:, None], 
                               Putative_HSPs[['tend','tstart']], 
                               Putative_HSPs[['tstart','tend']])
grouped = Putative_HSPs.groupby(['query', 'target'])
#Merge sequences together 
for group_name, group in grouped:
  new_tstart=min(group['tstart'])
  new_tend=max(group['tend'])
  new_qcov=sum(group['qcov'])
  diff=row['diff_length']
  Group=row['HSP_group']
  for index, row in group.iterrows():
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']),"new_tstart"]= new_tstart
    Putative_HSPs.loc[(Putative_HSPs['query']==row['query']) & (Putative_HSPs['target']==row['target']),"new_tend"]= new_tend
    PoFV_vs_filamentous_ORFs_cov.loc[(PoFV_vs_filamentous_ORFs_cov['query']==row['query']) & (PoFV_vs_filamentous_ORFs_cov['target']==row['target']),"diff_length"]= diff
    PoFV_vs_filamentous_ORFs_cov.loc[(PoFV_vs_filamentous_ORFs_cov['query']==row['query']) & (PoFV_vs_filamentous_ORFs_cov['target']==row['target']),"HSP_group"]= Group
    PoFV_vs_filamentous_ORFs_cov.loc[(PoFV_vs_filamentous_ORFs_cov['query']==row['query']) & (PoFV_vs_filamentous_ORFs_cov['target']==row['target']),"qcov"]= new_qcov
    if row['strand']=="+":
      PoFV_vs_filamentous_ORFs_cov.loc[(PoFV_vs_filamentous_ORFs_cov['query']==row['query']) & (PoFV_vs_filamentous_ORFs_cov['target']==row['target']),"tstart"]= new_tstart
      PoFV_vs_filamentous_ORFs_cov.loc[(PoFV_vs_filamentous_ORFs_cov['query']==row['query']) & (PoFV_vs_filamentous_ORFs_cov['target']==row['target']),"tend"]= new_tend
    if row['strand']=="-":
      PoFV_vs_filamentous_ORFs_cov.loc[(PoFV_vs_filamentous_ORFs_cov['query']==row['query']) & (PoFV_vs_filamentous_ORFs_cov['target']==row['target']),"tend"]= new_tstart
      PoFV_vs_filamentous_ORFs_cov.loc[(PoFV_vs_filamentous_ORFs_cov['query']==row['query']) & (PoFV_vs_filamentous_ORFs_cov['target']==row['target']),"tstart"]= new_tend
#Filter on query coverage
PoFV_vs_filamentous_ORFs_cov=PoFV_vs_filamentous_ORFs_cov[PoFV_vs_filamentous_ORFs_cov['qcov'].ge(0.25)]


PoFV_vs_filamentous_ORFs_cov['strand']=np.where(PoFV_vs_filamentous_ORFs_cov["tstart"]>PoFV_vs_filamentous_ORFs_cov["tend"],'-','+')

m = PoFV_vs_filamentous_ORFs_cov['strand'].eq('-')
PoFV_vs_filamentous_ORFs_cov.loc[m, ['tstart','tend']] = PoFV_vs_filamentous_ORFs_cov.loc[m, ['tend','tstart']].to_numpy()

#Remove overlapping loci 
is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
PoFV_vs_filamentous_ORFs_cov['overlapping_group'] = PoFV_vs_filamentous_ORFs_cov.sort_values(['target', 'tstart', 'tend']) \
                .groupby('target').apply(is_overlapped).droplevel(0).cumsum()
#
PoFV_vs_filamentous_ORFs_cov = PoFV_vs_filamentous_ORFs_cov.sort_values(['overlapping_group', 'evalue', 'alnlen'], ascending=[True, True, False]) \
        .groupby(PoFV_vs_filamentous_ORFs_cov['overlapping_group']).head(1)
#
# Putative_HSPs.loc[(Putative_HSPs['target'].str.contains("scaffold_983")]
# Extract all predicted sequences 
Assembly_record = SeqIO.to_dict(SeqIO.parse(Porseoliae_assembly, "fasta"))
All_predicted_sequences="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_sequences.faa"
with open(All_predicted_sequences,"w") as output:
  for index, row in PoFV_vs_filamentous_ORFs_cov.drop_duplicates(subset=['query', 'target','tstart','tend'], keep='first').iterrows():
      if row['tstart'] > row['tend']:
        #minus strand
        print(">PoFV_",row['target']+ ":" + str(row['tstart'])+"-"+str(row['tend'])+"(-)",sep="",file=output)
        print(str(Assembly_record[row['target']][row['tend']-1:row['tstart']].seq.reverse_complement().translate()),file=output)
      if row['tstart'] < row['tend']:
        #plus strand 
        print(">PoFV_",row['target']+ ":" + str(row['tstart'])+"-"+str(row['tend'])+"(+)",sep="",file=output)
        print(str(Assembly_record[row['target']][row['tstart']-1:row['tend']].seq.translate()),file=output)



#########
# save2 #
#########

PoFV_vs_filamentous_ORFs_cov_HSPs= PoFV_vs_filamentous_ORFs_cov.copy()
PoFV_vs_filamentous_ORFs_cov_HSPs.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs.m8",sep=";",index=False)
PoFV_vs_filamentous_ORFs_cov_HSPs= pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs.m8",sep=";")


### run NR analysis ### 

All_predicted_sequences_db="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_sequences_db"
All_predicted_sequences_vs_NR_resut="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_sequences_vs_NR_result"
All_predicted_sequences_vs_NR_temp="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_sequences_vs_NR_tpm"
All_predicted_sequences_vs_NR_table="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_sequences_vs_NR_result.m8"
NR_db="/beegfs/data/bguinet/these/NR_db"

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ All_predicted_sequences + " "+ All_predicted_sequences_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search + " "+  All_predicted_sequences_db NR_db + " "+  All_predicted_sequences_vs_NR_resut + " "+ All_predicted_sequences_vs_NR_temp + " -e 0.000001 --threads 10", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,taxid,taxname,taxlineage' + " "+  All_predicted_sequences_db NR_db + " "+  All_predicted_sequences_vs_NR_resut + " "+ All_predicted_sequences_vs_NR_table + " --threads 10", shell=True)


#############################
### Filter non-viral loci ###
#############################

NR_blast=pd.read_csv(All_predicted_sequences_vs_NR_table,sep="\t")
NR_blast.columns=["query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","qlen","tlen","tcov","taxid","taxname","taxlineage"]
NR_blast=NR_blast.loc[NR_blast['bits'].ge(50)]

#Take into account only non-hymenoptera matches 

NR_blast=NR_blast.loc[~NR_blast['taxlineage'].str.contains("Hymeno")]

NR_blast['count_Bacteria'] = NR_blast['taxlineage'].str.contains('Bacter').groupby(NR_blast['query']).transform('sum')
NR_blast['count_Archaea'] = NR_blast['taxlineage'].str.contains('Archaea').groupby(NR_blast['query']).transform('sum')
NR_blast['count_Eukaryota'] = NR_blast['taxlineage'].str.contains('Eukaryota').groupby(NR_blast['query']).transform('sum')
NR_blast['count_Virus'] = NR_blast['taxlineage'].str.contains('Virus').groupby(NR_blast['query']).transform('sum')
NR_blast['Count_Bacteria_Eukaryota_Archeae'] = NR_blast['count_Eukaryota'] + NR_blast['count_Bacteria']+ NR_blast['count_Archaea']
NR_blast['Perc_viral_vs_other'] = NR_blast['count_Virus'] / NR_blast['Count_Bacteria_Eukaryota_Archeae']
NR_blast['Viral_confidence'] = "NA"
NR_blast.loc[NR_blast['Count_Bacteria_Eukaryota_Archeae'].eq(0),"Viral_confidence"] = "Viral"
NR_blast.loc[NR_blast['Count_Bacteria_Eukaryota_Archeae'].gt(10) & NR_blast['Perc_viral_vs_other'].gt(1),"Viral_confidence"] = "Uncertain"
NR_blast.loc[NR_blast['Perc_viral_vs_other'].gt(1) & NR_blast['Viral_confidence'].str.contains("NA") ,"Viral_confidence"] = "Viral"
NR_blast.loc[NR_blast['Perc_viral_vs_other'].lt(1) & NR_blast['count_Virus'].eq(0) ,"Viral_confidence"] = "Not_viral"
NR_blast.loc[ NR_blast['count_Virus'].eq(0)  & NR_blast['Count_Bacteria_Eukaryota_Archeae'].eq(0) ,"Viral_confidence"] = "Viral"
NR_blast.loc[ NR_blast['count_Virus'].eq(0)  & NR_blast['Count_Bacteria_Eukaryota_Archeae'].lt(5) ,"Viral_confidence"] = "Viral"
NR_blast.loc[NR_blast['Viral_confidence'].str.contains("NA") ,"Viral_confidence"] = "Uncertain"
NR_blast.loc[NR_blast['count_Eukaryota'].gt(100) ,"Viral_confidence"] = "Not_viral"
NR_blast.loc[NR_blast['Viral_confidence'].isna() ,"Viral_confidence"] = "Viral"

Loci_to_remove=NR_blast.loc[NR_blast['Viral_confidence'].eq("Not_viral")]['query'].unique()

#Add full name 

for index, row in PoFV_vs_filamentous_ORFs_cov_HSPs.iterrows():
  if row['strand']=="-":
    full_loci_name="PoFV_"+row['target']+":"+str(row['tend'])+"-"+str(row['tstart'])+"(-)"
  else: 
    full_loci_name="PoFV_"+row['target']+":"+str(row['tstart'])+"-"+str(row['tend'])+"(+)"
  PoFV_vs_filamentous_ORFs_cov_HSPs.loc[(PoFV_vs_filamentous_ORFs_cov_HSPs['target']==row['target']) & (PoFV_vs_filamentous_ORFs_cov_HSPs['tstart']==row['tstart']),'full_loci_name']= full_loci_name

PoFV_vs_filamentous_ORFs_cov_HSPs=PoFV_vs_filamentous_ORFs_cov_HSPs.loc[~PoFV_vs_filamentous_ORFs_cov_HSPs['full_loci_name'].isin(Loci_to_remove)]


#########
# save3 #
#########
PoFV_vs_filamentous_ORFs_cov_HSPs_NR=PoFV_vs_filamentous_ORFs_cov_HSPs.copy()
PoFV_vs_filamentous_ORFs_cov_HSPs_NR.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR.m8",sep=";",index=False)
PoFV_vs_filamentous_ORFs_cov_HSPs_NR=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR.m8",sep=";")

# Write new filtred sequences 

PoFV_ORFs="/beegfs/data/bguinet/LbFV_family_project/Genomes/PoFV/Predicted_orfs/Final_ORF_prediction_PoFV.faa"
LbFV_ORFs="/beegfs/data/bguinet/LbFV_family_project/Genomes/LbFV/Predicted_orfs/Previously_predicted_108_ORFs.faa"
DFV_ORFs="/beegfs/data/bguinet/LbFV_family_project/Genomes/DFV/Predicted_orfs/Previously_predicted_68_ORFs.faa"
ALL_virus_ORFs="/beegfs/data/bguinet/LbFV_family_project/Clustering/All_known_viral_ORFs.faa"

Porseoliae_assembly=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/Platygaster_orseoliae_corrected2.fa", "fasta"))

Filtred_loci="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_filtred_sequences.faa"

with open(Filtred_loci,"w") as output:
  for index, row in PoFV_vs_filamentous_ORFs_cov_HSPs_NR.iterrows():
    if row['strand']=="+":
      print(">",row['full_loci_name'],sep="",file=output)
      print(Porseoliae_assembly[row['target']].seq[row['tstart']-1:row['tend']].translate(),file=output)
    if row['strand']=="-":
      print(">",row['full_loci_name'],sep="",file=output)
      print(Porseoliae_assembly[row['target']].seq[row['tstart']-1:row['tend']].reverse_complement().translate(),file=output)
  for fasta in SeqIO.parse(PoFV_ORFs, "fasta"):
    print(">",fasta.id,sep="",file=output)
    print(fasta.seq,file=output)
  for fasta in SeqIO.parse(LbFV_ORFs, "fasta"):
    print(">",fasta.id,sep="",file=output)
    print(fasta.seq,file=output)
  for fasta in SeqIO.parse(DFV_ORFs, "fasta"):
    print(">",fasta.id,sep="",file=output)
    print(fasta.seq,file=output)
  for fasta in SeqIO.parse(ALL_virus_ORFs, "fasta"):
    if "DFV" in fasta.id:
      continue 
    elif "LbFV" in fasta.id:
      continue 
    else:  
      print(">",fasta.id,sep="",file=output)
      print(fasta.seq,file=output)


##################
# Run clustering #
##################

Filtred_loci_db="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_filtred_sequences_db"
Cluster_result="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clustering_sequences_result"
Cluster_temp="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clustering_sequences_tpm"
Cluster_table="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clustering_sequences.m8"


subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Filtred_loci + " "+ Filtred_loci_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search  "+  Filtred_loci_db  +" "+ Filtred_loci_db + " "+  Cluster_result + " "+ Cluster_temp + " -e 0.000001 --threads 10", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,qcov'  "+  Filtred_loci_db +" "+ Filtred_loci_db + " "+  Cluster_result + " "+ Cluster_table + " --threads 10", shell=True)


# Open cluster table 

df=pd.read_csv(Cluster_table,sep="\t",header=None)
df.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen','tcov','qcov']

df=df.loc[df['qcov'].ge(0.15) | df['tcov'].ge(0.15) ]
df=df.loc[df['bits'].gt(50)]


# Create the graph from the dataframe
g = nx.Graph()

g.add_edges_from(df[['query','target']].itertuples(index=False))
new = list(nx.connected_components(g))

mapped =  {node: f'Cluster{cid + 1}' for cid, component in enumerate(new) for node in component}

cluster = pd.DataFrame({'Cluster': mapped.values(), 'Names':mapped.keys()})

#Add alignment informations 
# keep on query side all non PoEFV loci and on target side all PoEFV 
cluster_with_ali_info=cluster.copy()

df2=df.copy()
df2.columns=['target','query','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']

df=df.loc[~df['query'].str.contains("\\(\\+|\\(-")]
df=df.loc[df['target'].str.contains("\\(\\+|\\(-")]
df2=df2.loc[~df2['query'].str.contains("\\(\\+|\\(-")]
df2=df2.loc[df2['target'].str.contains("\\(\\+|\\(-")]

df3=df.append(df2)
df3=df3.loc[~(df3['query'] == df3['target'])]

df3 = df3.drop_duplicates(subset=['query', 'target'], keep='first')


#Merge informations 
df4=df3.merge(cluster,left_on="target",right_on="Names",how="left")

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters=df4.merge(PoFV_vs_filamentous_ORFs_cov_HSPs_NR[['scaffolds', 'GC', 'Arithmetic_cov_depth', 'Scaf_length', 'Median_cov_depth', 'scaf_name', 'scaf_length', 'cov_depth_candidat', 'cov_depth_BUSCO', 'pvalue_cov', 'count_duplicate', 'diff_length', 'HSP_group', 'strand', 'overlapping_group', 'full_loci_name']],right_on="full_loci_name",left_on="Names")
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters = PoFV_vs_filamentous_ORFs_cov_HSPs_clusters.drop_duplicates()

######
#save4
######
cluster.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters_for_msa_files.tab",sep=";")
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters.m8",sep=";",index=False)
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters.m8",sep=";")


#####################################################
# Find ORFs along the scaffolds with candidate loci #
#####################################################


# Write PoEFV EVEs scaffolds into a file 

PoEFV_scaffolds= "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds.aa"
Porseoliae_assembly=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/Platygaster_orseoliae_corrected2.fa", "fasta"))

with open(PoEFV_scaffolds,"w") as output:
  for scaff in PoFV_vs_filamentous_ORFs_cov_HSPs_clusters['scaf_name'].unique():
    print(">",Porseoliae_assembly[scaff].id,sep="",file=output)
    print(Porseoliae_assembly[scaff].seq,file=output)

#

Filtred_PoEFV_loci="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Porseoliae_filtred_sequences.faa"


for loci in SeqIO.parse(Filtred_PoEFV_loci,"fasta"):
with open(Filtred_PoEFV_loci,"w") as output:
  for index, row in PoFV_vs_filamentous_ORFs_cov_HSPs_NR.iterrows():
    if row['strand']=="+":
      print(">",row['full_loci_name'],sep="",file=output)
      print(Porseoliae_assembly[row['target']].seq[row['tstart']-1:row['tend']].translate(),file=output)
    if row['strand']=="-":
      print(">",row['full_loci_name'],sep="",file=output)
      print(Porseoliae_assembly[row['target']].seq[row['tstart']-1:row['tend']].reverse_complement().translate(),file=output)
##

# Run python ORF finder 
import pathlib  

subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/orfipy "+ PoEFV_scaffolds +" --ignore-case  --procs 5 --outdir /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/ --bed /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/scaffold_orfipy.bed --min 150 --start ATG --dna /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/scaffold_orfipy.fna", shell=True)
# Translate 

ORF_bed=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/scaffold_orfipy.bed",sep="\t",header=None)
ORF_bed.columns=['Scaffold_name','ORF_start','ORF_end','ORF_name','zero','ORF_strand']

ORF_bed['ORF_name2']=ORF_bed['ORF_name'].str.replace(";.*","")
ORF_bed['ORF_name2']=ORF_bed['ORF_name2'].str.replace("ID=","")

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/scaffold_orfipy.faa","w") as output:
        record_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/scaffold_orfipy.fna","fasta"))
        for species in ORF_bed['Scaffold_name'].unique():
                for index, row in ORF_bed.loc[ORF_bed['Scaffold_name'].str.contains(species)].iterrows():
                        print(">",row['ORF_name2'],';',row['ORF_start'],"-",row['ORF_end'],"(",row['ORF_strand'],")",sep="",file=output)
                        print(record_dict[row['ORF_name2']].seq.translate(),file=output)

# Run mmseqs between ORFs and previously selected candidates

Predicted_ORFs="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/scaffold_orfipy.faa"
Predicted_ORFs_db="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/scaffold_orfipy_db"
Filtred_PoEFV_loci="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Porseoliae_filtred_sequences.faa"
Filtred_PoEFV_loci_db="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Porseoliae_filtred_sequences_db"

Filtred_PoEFV_loci_vs_Predicted_ORFs_result="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filtred_PoEFV_loci_vs_Predicted_ORFs_result"
Filtred_PoEFV_loci_vs_Predicted_ORFs_temp="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filtred_PoEFV_loci_vs_Predicted_ORFs_tpm"
Filtred_PoEFV_loci_vs_Predicted_ORFs_table="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filtred_PoEFV_loci_vs_Predicted_ORFs_result.m8"

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Predicted_ORFs + " "+ Predicted_ORFs_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Filtred_PoEFV_loci + " "+ Filtred_PoEFV_loci_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search  "+  Filtred_PoEFV_loci_db  +" "+ Predicted_ORFs_db + " "+  Filtred_PoEFV_loci_vs_Predicted_ORFs_result + " "+ Filtred_PoEFV_loci_vs_Predicted_ORFs_temp + " -e 0.000001 --threads 10", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,qcov'  "+  Filtred_PoEFV_loci_db +" "+ Predicted_ORFs_db + " "+  Filtred_PoEFV_loci_vs_Predicted_ORFs_result + " "+ Filtred_PoEFV_loci_vs_Predicted_ORFs_table + " --threads 10", shell=True)

### open the mmseqs file 

ORF_vs_EVE_table=pd.read_csv(Filtred_PoEFV_loci_vs_Predicted_ORFs_table,sep="\t",header=None)
ORF_vs_EVE_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen','tcov','qcov']
# Keep only self matching hits 
ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['bits'].ge(50)]

ORF_vs_EVE_table['query_EVE_species']=ORF_vs_EVE_table['query'].str.replace(".*:","")
ORF_vs_EVE_table['query_EVE_species'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV"),ORF_vs_EVE_table['query'].str.replace(".*_",""),inplace=True)
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['target'].str.replace(".*;","")
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("\\(\\+\\)","")
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("\\(-\\)","")
ORF_vs_EVE_table['ORF_strand']="NA"
ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query'].str.contains("\\+"),'ORF_strand'] ="+"
ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query'].str.contains("-"),'ORF_strand'] ="-"
#ORF_vs_EVE_table['target'].str.replace(".*;","")
ORF_vs_EVE_table['ORF_end']=ORF_vs_EVE_table['ORF_start'].str.replace(".*-","").astype(int)
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("-.*","").astype(int)
ORF_vs_EVE_table['target_ORF_species']=ORF_vs_EVE_table['target'].str.replace(":.*","")
ORF_vs_EVE_table['query_EVE_scaffold']=ORF_vs_EVE_table['query'].str.replace(":.*","")
ORF_vs_EVE_table['query_EVE_scaffold'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV") , ORF_vs_EVE_table['query'].str.replace("-.*",""),inplace=True)
ORF_vs_EVE_table['query_EVE_scaffold'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV") , ORF_vs_EVE_table['query_EVE_scaffold'].str.split('_').str[:-1].str.join('_'),inplace=True)
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target'].str.replace(".*:","")
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target_ORF_scaffold'].str.replace("_ORF.*","")
ORF_vs_EVE_table['query_EVE_scaffold']=ORF_vs_EVE_table['query_EVE_scaffold'].str.replace("PoFV_","")

ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query_EVE_scaffold'] == ORF_vs_EVE_table['target_ORF_scaffold']]

# Remove duplicate 
ORF_vs_EVE_table = ORF_vs_EVE_table.drop_duplicates()

# find overlapping ORFs within the same scaffold
ORF_vs_EVE_table.rename(columns={'query': 'ORF_query',
                   'target': 'ORF_target'},
          inplace=True, errors='raise')

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters.merge(ORF_vs_EVE_table[['ORF_query','ORF_target','ORF_start','ORF_end','ORF_strand']],left_on="Names",right_on="ORF_query",how="left")


import numpy as np

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['Overlapp_ORF_EVEs']= np.nan
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_perc']= np.nan

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['start-end']= PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['target'].str.replace("\\(.*","")
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['start-end']= PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['start-end'].str.replace(".*:","")
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['start']= PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['start-end'].str.replace("-.*","")
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['end']= PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['start-end'].str.replace(".*-","")

m = PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['strand'].eq('-')
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs[['start','end']] = np.where(m.to_numpy()[:, None], 
                               PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs[['end','start']], 
                               PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs[['start','end']])
                               
for index, row in PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_start'].ge(0)].iterrows():
        overlapp="no"
        if (int(row['start']) <= int(row['ORF_start'])) & (int(row['end']) <= int(row['ORF_end'])) & (int(row['end']) > int(row['ORF_start'])):
                #print("left")
                #print(row['start']," : ", row['end'])
                #print(row['ORF_start']," : ", row['ORF_end'])
                overlapp = "yes-left"
        elif (int(row['start']) >= int(row['ORF_start'])) &  (int(row['end']) >= int(row['ORF_end'])) & (int(row['start']) < int(row['ORF_end'])):
                #print("right")
                #print(row['start']," : ", row['end'])
                #print(row['ORF_start']," : ", row['ORF_end'])
                overlapp = "yes-right"
        elif (int(row['start']) >= int(row['ORF_start'])) &  (int(row['end']) <= int(row['ORF_end'])):
                overlapp = "yes-inside"
        elif (int(row['start']) <= int(row['ORF_start'])) & (int(row['end']) >= int(row['ORF_end'])) :
                #print("outside")
                #print(row['start']," : ", row['end'])
                #print(row['ORF_start']," : ", row['ORF_end'])
                overlapp = "yes-outside"
        else:
                overlapp = "no"
        PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_query'].eq(row['ORF_query']) & PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_target'].eq(row['target']) & PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_end'].eq(row['ORF_end']) & PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['end'].eq(row['end']),"Overlapp_ORF_EVEs"]=overlapp
        if overlapp =="no":
                continue
        else:
                PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_query'].eq(row['ORF_query']) & PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_target'].eq(row['ORF_target']) & PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_end'].eq(row['ORF_end']) & PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['end'].eq(row['end']),"ORF_perc"]= (int(row['ORF_end'])- int(row['ORF_start']))/ (int(row['end'])- int(row['start']))
        #print("\n")
#

# remove non-overlapping ORFS with EVEs 
import numpy as np
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['Overlapp_ORF_EVEs'].eq("no"),"ORF_start"] = np.nan
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['Overlapp_ORF_EVEs'].eq("no"),"ORF_end"] = np.nan

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['Overlapp_ORF_EVEs'].eq("yes-outside") & PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_perc'].lt(0.5),"ORF_start"] = np.nan
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['Overlapp_ORF_EVEs'].eq("yes-outside") & PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_perc'].lt(0.5),"ORF_end"] = np.nan

# find overlapping ORFs within the same scaffold

is_overlapped = lambda x: x['ORF_start'] >= x['ORF_end'].shift(fill_value=-1)
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['Overlapp_group'] = PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.sort_values(['scaf_name', 'ORF_start', 'ORF_end']) \
                .groupby(['scaf_name'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['ORF_start'].isna(),"Overlapp_group"] = np.nan

# Add if loci is pseudogenized of not 
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['Pseudogenized']="no"
Filtred_PoEFV_loci_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Porseoliae_filtred_sequences.faa","fasta"))

for index, row in PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.iterrows():
  if "*" in str(Filtred_PoEFV_loci_dict[row['full_loci_name']].seq):
    PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs['full_loci_name']==row['full_loci_name'],"Pseudogenized"]="yes"

#########
# save5 #
#########

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs.m8",sep=";",index=False)
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs.m8",sep=";")

#############
####Run Augustus and search against uniprot  #####
####Run search in scaffolds against Repeatpeps ####

################
###Add repeat ##
################
import pandas as pd 
import numpy as np 

PoEFV_scaffolds= "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds.aa"
PoEFV_scaffolds_db= "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_db"
Repeat_db="/beegfs/data/bguinet/these/Repeat_env_analysis2/RepeatPeps_db"

PoEFV_loci_vs_repeat_result="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_repeat_result"
PoEFV_loci_vs_repeat_temp="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_repeat_tpm"
PoEFV_loci_vs_repeat_table="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_repeat_result.m8"

# Run mmseqs against repeat db 

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ PoEFV_scaffolds + " "+ PoEFV_scaffolds_db , shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search  "+  PoEFV_scaffolds_db  +" "+ Repeat_db + " "+  PoEFV_loci_vs_repeat_result + " "+ PoEFV_loci_vs_repeat_temp + " -e 0.000001 --threads 10", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,qcov'  "+  PoEFV_scaffolds_db +" "+ Repeat_db + " "+  PoEFV_loci_vs_repeat_result + " "+ PoEFV_loci_vs_repeat_table + " --threads 10", shell=True)

Repeat_tab=pd.read_table(PoEFV_loci_vs_repeat_table,sep="\t")
Repeat_tab.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen','tcov','qcov']
sLength=Repeat_tab.shape[0]
Repeat_tab['strand']=np.where(Repeat_tab["qstart"]>Repeat_tab["qend"],'-','+')
Repeat_tab[['qstart', 'qend']] = np.sort(Repeat_tab[['qstart', 'qend']].values, axis=1)
m = Repeat_tab['strand'].eq('-')
Repeat_tab['Newqstart'] = np.where(m, Repeat_tab['qlen'].sub(Repeat_tab['qend']), Repeat_tab['qstart'])
Repeat_tab['Newqend'] = np.where(m, Repeat_tab['qlen'].sub(Repeat_tab['qstart']), Repeat_tab['qend'])

Repeat_tab = Repeat_tab.sort_values(['query', 'qend', 'qstart'])

c1 = Repeat_tab['query'].shift() != Repeat_tab['query']
c2 = Repeat_tab['qend'].shift() - Repeat_tab['qstart'] < 0
Repeat_tab['overlap'] = (c1 | c2).cumsum()
#Finally, we get the row with the maximum sum in each group using groupby.

Repeat_tab['qmatchlen'] = Repeat_tab['Newqend'].astype(int) - Repeat_tab['Newqstart'].astype(int)
Repeat_tab=Repeat_tab.sort_values(['evalue'], ascending=True).groupby('overlap').first()

#Remove every match < 100 pb 
Repeat_tab=Repeat_tab.loc[Repeat_tab['qmatchlen'].ge(100)]
Repeat_tab=Repeat_tab.loc[Repeat_tab['evalue'].lt(0.0000000001)]
#Remove unknown annotations 
Repeat_tab=Repeat_tab.loc[~Repeat_tab['target'].str.contains("Unknown")]

Repeat_tab['count_repeat'] = Repeat_tab.groupby('query')['query'].transform('count')
Repeat_tab=Repeat_tab[['query','count_repeat']]
Repeat_tab.columns=['scaf_name','count_repeat']
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs.merge(Repeat_tab,on=['scaf_name'],how="outer")

#########
# save6 #
#########
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat.m8",sep=";",index=False)
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat.m8",sep=";")

######################
###Add metaeuk  #####
######################

PoEFV_scaffolds= "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds.aa"
PoEFV_scaffolds_db= "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_db"
UniProtKB_db="/beegfs/data/bguinet/TOOLS/UniProtKB"

PoEFV_scaffolds_db_vs_metaeuk_result="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_metaeuk_result"
PoEFV_scaffolds_db_vs_metaeuk_temp="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_metaeuk_tpm"
PoEFV_scaffolds_db_vs_metaeuk_prediction="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_metaeuk_result.fas"
PoEFV_scaffolds_db_vs_metaeuk_map="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_metaeuk_result.headersMap.tsv"


subprocess.run("/beegfs/data/bguinet/TOOLS/metaeuk/build/bin/metaeuk  createdb "+ PoEFV_scaffolds  + " "+ PoEFV_scaffolds_db, shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/metaeuk/build/bin/metaeuk easy-predict "+ PoEFV_scaffolds_db  + " "+ UniProtKB_db  + " "+ PoEFV_scaffolds_db_vs_metaeuk_result  + " "+ PoEFV_scaffolds_db_vs_metaeuk_temp  +" --threads 10", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/metaeuk/build/bin/metaeuk taxtocontig "+ PoEFV_scaffolds_db  + " "+ PoEFV_scaffolds_db_vs_metaeuk_prediction  + " "+ PoEFV_scaffolds_db_vs_metaeuk_map  + " "+ UniProtKB_db  + " "+ PoEFV_scaffolds_db_vs_metaeuk_result  + " "+ PoEFV_scaffolds_db_vs_metaeuk_temp +" --majority 0.5 --tax-lineage 1 --lca-mode 2 --threads 10", shell=True)

#incorporate metaeuk result 

# Filter metaeuk results 
Metaeuk_per_pred_table=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_metaeuk_result_tax_per_pred.tsv",sep="\t",header=None)
Metaeuk_per_pred_table.columns=['Protein','taxid','rank','prediction','taxlineages']
Metaeuk_per_pred_table=Metaeuk_per_pred_table.loc[Metaeuk_per_pred_table['taxlineages'].str.contains("Euk",na=False)]
Metaeuk_per_pred_table['Scaffold_and_Species_name']=Metaeuk_per_pred_table['Protein'].str.replace("\\|\\+.*","")
Metaeuk_per_pred_table['Scaffold_and_Species_name']=Metaeuk_per_pred_table['Scaffold_and_Species_name'].str.replace("\\|-.*","")
Metaeuk_per_pred_table['Scaffold_and_Species_name']=Metaeuk_per_pred_table['Scaffold_and_Species_name'].str.split('|', n=1).str.get(-1)
Metaeuk_per_pred_table_count=Metaeuk_per_pred_table.groupby(['Scaffold_and_Species_name']).size().reset_index(name='count_euk')

Metaeuk_per_contig_table=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds_vs_metaeuk_result_tax_per_contig.tsv",sep="\t",header=None)
Metaeuk_per_contig_table=Metaeuk_per_contig_table.loc[Metaeuk_per_contig_table[8].str.contains("Eukaryot")]
Metaeuk_per_contig_table['Metaeuk_vote']="yes"
Metaeuk_per_contig_table=Metaeuk_per_contig_table[[0,'Metaeuk_vote']]
Metaeuk_per_contig_table.columns=['Scaffold_and_Species_name','Metaeuk_vote']
Metaeuk_per_pred_table_count=Metaeuk_per_pred_table_count.merge(Metaeuk_per_contig_table,on="Scaffold_and_Species_name",how="outer")

Metaeuk_per_pred_table_count=Metaeuk_per_pred_table_count[['Scaffold_and_Species_name','count_euk']]
Metaeuk_per_pred_table_count.columns=['scaf_name','count_euk']
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat.merge(Metaeuk_per_pred_table_count,on=['scaf_name'],how="outer")

#########
# save7 #
#########
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat_euk.m8",sep=";",index=False)
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat_euk.m8",sep=";")

################
### Add scores #
################

# Merge with the main table result 
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.count_repeat.fillna(0,inplace=True)
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.count_euk.fillna(0,inplace=True)

# Calculate the scaffolds scores 

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['count_euk'].eq(0) | PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['count_euk'].eq()  ,'Scaffold_score']='C'
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['count_repeat'].eq(0) |PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['count_repeat'].eq(0)  ,'Scaffold_score']='C'

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['count_euk'].ge(1) | PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['count_euk'].ge(1)  ,'Scaffold_score']='A'
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['count_repeat'].ge(1) |PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['count_repeat'].ge(1)  ,'Scaffold_score']='A'

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.copy()


#########
# save8 #
#########
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat_euk.m8",sep=";",index=False)

import pandas as pd 
import re,os 
import subprocess
from Bio import SeqIO 

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat_euk.m8",sep=";")

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['loci_length']=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['end']-PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['start']


#Remove overlapping loci 
is_overlapped = lambda x: x['start'] >= x['end'].shift(fill_value=-1)
PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['overlapping_group'] = PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.sort_values(['scaf_name', 'start', 'end']) \
                .groupby('scaf_name').apply(is_overlapped).droplevel(0).cumsum()
#
Loci_to_keep = PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.sort_values(['overlapping_group', 'loci_length'], ascending=[True, False]) \
        .groupby(PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['overlapping_group']).head(1)

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['full_loci_name'].isin(Loci_to_keep['full_loci_name'].unique())]
      
############################
##### Create MSA files #####
############################

Cluster_file=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters_for_msa_files.tab",sep=";")
grouped = PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.groupby('Cluster')

Filtred_loci=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_filtred_sequences.faa","fasta"))

PoEFV_scaffolds= "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds.aa"
Scaffold_assembly_record = SeqIO.to_dict(SeqIO.parse(PoEFV_scaffolds, "fasta"))

from Bio import SeqIO 
def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

PoFV_assembly = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/PoFV/Predicted_orfs/Final_ORF_prediction_PoFV.fna","fasta"))

from Bio import Entrez

#This table was generate using : https://www.uniprot.org/id-mapping/
uniprot_id_to_refseq_id_tab=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Uniprot_id_to_refseq_id.tab",sep="\t",header=0)
# iterate over each group
for group_name, group in grouped:
    print("Adding ...", group_name)
    sub_Cluster_file=Cluster_file.loc[Cluster_file['Cluster'].eq(group_name)]
    #Add DNA sequences 
    with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+group_name+".dna",'w') as output:
      #ADD PoFV loci 
      for row_index, row in sub_Cluster_file.iterrows():
        if "_orf" in row['Names']:
          if "PoFV" in row['Names']:
            print('>',row['Names'],sep="",file=output)
            print(PoFV_assembly[row['Names']].seq,file=output)
        else: 
          if "(+)" in row['Names']:
            continue 
          elif "(-)" in row['Names']:
            continue 
          else:
            #ADD other virus loci 
            for name in ['DFV','LbFV','DmNV_tom','LdMNPV', 'AcMNPV', 'OrNV', 'CpV', 'NeseNPV', 'ToNV', 'DhNV', 'HzNV-1','MdSGHV', 'GbNV', 'GpSGHV', 'PmNV', 'CuniNPV', 'CoBV', 'AmFV', 'WSSV']:
              if name in row['Names']:
                id=re.sub(name,"",row['Names'])
                id=re.sub("_YP_","YP)",id)
                id=re.sub("_","",id)
                id=re.sub("\\..*","",id)
                id=re.sub("YP\\)","YP_",id)
                try:
                  handle=Entrez.efetch(db="sequences", id=id, rettype="fasta_cds_na", retmode="text")
                  record = SeqIO.read(handle, "fasta")
                except:  
                    try:
                      id=re.sub("\\..*","",uniprot_id_to_refseq_id_tab.loc[uniprot_id_to_refseq_id_tab['UNIPROT'].eq(id)]['REFSEQ'].iloc[0])
                      handle=Entrez.efetch(db="sequences", id=id, rettype="fasta_cds_na", retmode="text")
                      record = SeqIO.read(handle, "fasta")
                    except:
                      print("Sequence from ", name ," : ", id, " not found" )
                print(">",row['Names'],sep="",file=output)
                print(record.seq,file=output)
      #ADD PoEFV loci 
      for loci in group['full_loci_name'].unique():
        scaffold_name = group.loc[group['full_loci_name'].eq(loci)]['scaf_name'].iloc[0]
        loci_start = group.loc[group['full_loci_name'].eq(loci)]['start'].iloc[0]
        loci_end = group.loc[group['full_loci_name'].eq(loci)]['end'].iloc[0]
        strand = group.loc[group['full_loci_name'].eq(loci)]['strand'].iloc[0]
        Scaffold_score = PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['full_loci_name']==loci]['Scaffold_score'].iloc[0]
        new_loci_name=Filtred_loci[loci].id
        new_loci_name=re.sub("PoFV","PoEFV",new_loci_name)
        if group.loc[group['full_loci_name'].eq(loci)]['ORF_perc'].iloc[0] > 0.5:
          if group.loc[group['full_loci_name'].eq(loci)]['Pseudogenized'].iloc[0] == "no":
            new_loci_name=new_loci_name+"|[ORF]|"+Scaffold_score
        if group.loc[group['full_loci_name'].eq(loci)]['Pseudogenized'].iloc[0] == "yes":
          new_loci_name=new_loci_name+"|[Pseudogenized]|"+Scaffold_score
        if pd.isna(group.loc[group['full_loci_name'].eq(loci)]['ORF_perc'].iloc[0]):
            new_loci_name=new_loci_name+"|[noORF]|"+Scaffold_score
        new_loci_name=re.sub("\\(\\+\\)","_+_",new_loci_name)
        new_loci_name=re.sub("\\(-\\)","_-_",new_loci_name)
        new_loci_name=re.sub(":","_",new_loci_name)
        print(">",new_loci_name,sep="",file=output)
        if strand == '-':
          print(str(Scaffold_assembly_record[scaffold_name].seq[loci_start-1:loci_end].reverse_complement()),file=output)
        else:
          print(str(Scaffold_assembly_record[scaffold_name].seq[loci_start-1:loci_end]),file=output)
    #Add AA sequences 
    with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+group_name+".aa",'w') as output:
      #ADD PoFV loci 
      for row_index, row in sub_Cluster_file.iterrows():
        if "_orf" in row['Names']:
          if "PoFV" in row['Names']:
            print('>',row['Names'],sep="",file=output)
            print(PoFV_assembly[row['Names']].seq.translate(),file=output)
        else: 
          if "(+)" in row['Names']:
            continue 
          elif "(-)" in row['Names']:
            continue 
          else:
            #ADD other virus loci 
            for name in ['DFV','LbFV','DmNV_tom','LdMNPV', 'AcMNPV', 'OrNV', 'CpV', 'NeseNPV', 'ToNV', 'DhNV', 'HzNV-1','MdSGHV', 'GbNV', 'GpSGHV', 'PmNV', 'CuniNPV', 'CoBV', 'AmFV', 'WSSV']:
              if name in row['Names']:
                id=re.sub(name,"",row['Names'])
                id=re.sub("_YP_","YP)",id)
                id=re.sub("_","",id)
                id=re.sub("\\..*","",id)
                id=re.sub("YP\\)","YP_",id)
                try:
                  handle=Entrez.efetch(db="sequences", id=id, rettype="fasta_cds_na", retmode="text")
                  record = SeqIO.read(handle, "fasta")
                except:  
                    try:
                      id=re.sub("\\..*","",uniprot_id_to_refseq_id_tab.loc[uniprot_id_to_refseq_id_tab['UNIPROT'].eq(id)]['REFSEQ'].iloc[0])
                      handle=Entrez.efetch(db="sequences", id=id, rettype="fasta_cds_na", retmode="text")
                      record = SeqIO.read(handle, "fasta")
                    except:
                      print("Sequence from ", name ," : ", id, " not found" )
                print(">",row['Names'],sep="",file=output)
                print(record.seq.translate(),file=output)
      #ADD PoEFV loci 
      for loci in group['full_loci_name'].unique():
        scaffold_name = group.loc[group['full_loci_name'].eq(loci)]['scaf_name'].iloc[0]
        loci_start = group.loc[group['full_loci_name'].eq(loci)]['start'].iloc[0]
        loci_end = group.loc[group['full_loci_name'].eq(loci)]['end'].iloc[0]
        strand = group.loc[group['full_loci_name'].eq(loci)]['strand'].iloc[0]
        Scaffold_score = PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['full_loci_name']==loci]['Scaffold_score'].iloc[0]
        new_loci_name=Filtred_loci[loci].id
        new_loci_name=re.sub("PoFV","PoEFV",new_loci_name)
        if group.loc[group['full_loci_name'].eq(loci)]['ORF_perc'].iloc[0] > 0.5:
          if group.loc[group['full_loci_name'].eq(loci)]['Pseudogenized'].iloc[0] == "no":
            new_loci_name=new_loci_name+"|[ORF]|"+Scaffold_score
        if group.loc[group['full_loci_name'].eq(loci)]['Pseudogenized'].iloc[0] == "yes":
          new_loci_name=new_loci_name+"|[Pseudogenized]|"+Scaffold_score
        if pd.isna(group.loc[group['full_loci_name'].eq(loci)]['ORF_perc'].iloc[0]):
            new_loci_name=new_loci_name+"|[noORF]|"+Scaffold_score
        print(">",new_loci_name,sep="",file=output)
        if strand == '-':
          print(str(Scaffold_assembly_record[scaffold_name].seq[loci_start-1:loci_end].reverse_complement().translate()),file=output)
        else:
          print(str(Scaffold_assembly_record[scaffold_name].seq[loci_start-1:loci_end].translate()),file=output)
    print(group_name, " - Added")

# Run alignment with clustal 
subprocess.run("for file in /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/*.aa; do /beegfs/data/bguinet/TOOLS/clustalo -i $file -o $file.aln -v --threads=15 ; done" , shell=True)

# Run alignment with phylogeny  
subprocess.run("for file in /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/*.aa.aln; do /beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s $file -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15; done" , shell=True)

# Count number of cluster with more than 3 sequences 
Cluster_file=Cluster_file.merge(Cluster_file.groupby(['Cluster']).size().reset_index(name='Nb_loci'),on="Cluster")

len(Cluster_file.loc[Cluster_file['Nb_loci'].gt(3) & Cluster_file['Names'].str.contains("\\(")]['Cluster'].unique())
len(Cluster_file.loc[Cluster_file['Nb_loci'].eq(3) & Cluster_file['Names'].str.contains("\\(")]['Cluster'].unique())

# Manually assign event to clusters 

Cluster_with_2_paralogs=['Cluster149','Cluster183','Cluster207','Cluster301','Cluster379','Cluster406','Cluster407',
'Cluster417','Cluster419','Cluster427','Cluster445','Cluster670','Cluster697','Cluster922','Cluster1183','Cluster1245',
'Cluster1248','Cluster1263','Cluster1417']

Cluster_with_corrected_topology=['Cluster1','Cluster9','Cluster63','Cluster85','Cluster101',
'Cluster108', 'Cluster131','Cluster297','Cluster444', 'Cluster454', 'Cluster536']

Cluster_with_too_complicated_topology=['Cluster56','Cluster148']

Cluster_with_good_topology =["Cluster69","Cluster296",'Cluster339','Cluster355','Cluster477']

#### Manual stuff to constrain the trees and test them 

#Run trees with constrain 
#Cluster1

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster1_constrain.nwk","w") as output:
  print("((DFV_QKN22497.1,LbFV_YP_009345624.1),(((PoEFV_scaffold_27772_877-3070___|_Pseudogenized_|B,PoEFV_scaffold_32496_818-2053___|_ORF_|B)100/100,(PoEFV_scaffold_18363_4071-5051___|_Pseudogenized_|B,PoEFV_scaffold_48333_0_100-681___|_ORF_|B)95.1/96)86.4/98,PoFV_scaffold_21913_orf012)100/100);",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster1.aa.aln  --prefix Cluster1.constr1 -m VT+F+I  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster1_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster1.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster1.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster1_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster1.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster1_all.treefile -n 0 -zb 1000 -au --prefix Cluster1.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 

########
#Cluster9
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster9_constrain.nwk","w") as output:
  print("(LbFV_YP_009345682.1,(((PoEFV_scaffold_6401_16593-18479___|_Pseudogenized_|A,PoEFV_scaffold_7180_4204-1607_-_|_ORF_|B)100/100,PoFV_scaffold_15158_orf005)69/69,((B0YLI7_GpSGHV,B2YG51_MdSGHV)85.5/92,B0YLI6_GpSGHV:0.4899820241)100/100:0.9824638465)99.6/100:0.5452482146,DFV_QKN22481.1);",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster9.aa.aln  --prefix Cluster9.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster9_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster9.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster9.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster9_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster9.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster9_all.treefile -n 0 -zb 1000 -au --prefix Cluster9.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 


########
#Cluster63
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster63_constrain.nwk","w") as output:
  print("((((((((S5GFR9_AcMNPV,P41427_AcMNPV)100/100,Q9YMK8_LdMNPV)97.8/96,(A0A097P0N0_CpV,Q91F07_CpV)100/100)92.1/77,Q919N1_CuniNPV)33.5/48,Q9QAB4_NeseNPV)100/100,((A0A0B4VFN2_ToNV,(AAN04416_HzNV-1,(A0A7D5UN60_DhNV,A0A076FCX1_PmNV)99.1/100)93.4/76)21.1/46,((B7SV38_OrNV,ATY70204_DmNV_tom)99.7/100,A4L229_GbNV)99.9/100)98.1/99)96.2/99,AKY03169_AmFV),((LbFV_YP_009345656.1,((PoEFV_scaffold_75772_1153-689_-_|_ORF_|B,PoEFV_scaffold_8315_10729-10394_-_|_noORF_|B)2.9/59,PoFV_scaffold_7638_orf005)99.9/100)97.4/100,(B0YLK7_GpSGHV,B2YG66_MdSGHV)100/100))83.4/91;",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster63.aa.aln  --prefix Cluster63.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster63_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster63.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster63.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster63_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster63.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster63_all.treefile -n 0 -zb 1000 -au --prefix Cluster63.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 


########
#Cluster85
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster85_constrain.nwk","w") as output:
  print("(ATY70228_DmNV_tom,((((A0A7D5YSI3_DhNV,(A0A0B4VGM8_ToNV,AAN04447_HzNV-1)98.3/98)5/36,A0A076FCA4_PmNV)83.4/82,(((PoEFV_scaffold_7821_12565-13134___|_ORF_|B,PoEFV_scaffold_983_9781-10272___|_ORF_|B)0/44,PoFV_scaffold_2702_orf008)95.6/90,(LbFV_YP_009345623.1,DFV_QKN22470.1)77.6/57)99.8/100)99.3/93,A4L1W4_GbNV)95.6/94,B7SVA8_OrNV);",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster85.aa.aln  --prefix Cluster85.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster85_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster85.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster85.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster85_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster85.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster85_all.treefile -n 0 -zb 1000 -au --prefix Cluster85.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 


########
#Cluster101
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster101_constrain.nwk","w") as output:
  print("(Q9YMS8_LdMNPV,(((LbFV_YP_009345641.1,(DFV_QKN22460.1,((PoEFV_scaffold_32496_320-736___|_ORF_|B,PoEFV_scaffold_18363_3575-3991___|_ORF_|B)85.6/67,PoFV_scaffold_21913_orf011)99.9/100)43.1/47)96/93,B0YLK2_GpSGHV)92.2/89,(A0A097P1H6_CpV,Q91EY6_CpV)96.5/100)90.7/96,(A0A097PV38_AcMNPV,P21290_AcMNPV)99.9/100);",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster101.aa.aln  --prefix Cluster101.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster101_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster101.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster101.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster101_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster101.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster101_all.treefile -n 0 -zb 1000 -au --prefix Cluster101.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 

########
#Cluster108
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster108_constrain.nwk","w") as output:
  print("(A0A110AQ78_GpSGHV,(((((PoEFV_scaffold_7180_10129-9815_-_|_noORF_|B,PoEFV_scaffold_150917_435-734___|_ORF_|B)59.8/93:0.1415117877,PoFV_scaffold_15158_orf010)100/100,((((Q9YMI5_LdMNPV,(A0A097PUZ0_AcMNPV,P41668_AcMNPV)98.2/100)95.9/94,Q91F18_CpV)74.7/49,Q6JK91_NeseNPV)77.3/49,(Q77GU5_CuniNPV,Q99GR4_CuniNPV)99.9/100)87.5/78)84.7/70,(((ATY70238_DmNV_tom,B7SVC8_OrNV)97.2/99,A4L1W6_GbNV)96.4/100,(((AAN04382_HzNV-1,A0A7D5UN42_DhNV)32.4/45,A0A076FJ13_PmNV)98.4/96,A0A0B4VGI0_ToNV:0.7215693330)74.1/63)99.9/100)26/65,B2YG83_MdSGHV)99.5/100,B0YLN0_GpSGHV);",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster108.aa.aln  --prefix Cluster108.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster108_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster108.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster108.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster108_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster108.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster108_all.treefile -n 0 -zb 1000 -au --prefix Cluster108.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 


########
#Cluster131
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster131_constrain.nwk","w") as output:
  print("((DFV_QKN22511.1,LbFV_YP_009345701.1),(((PoEFV_scaffold_983_21453-21136_-_|_noORF_|B,(PoEFV_scaffold_22338_3694-4368___|_noORF_|B,PoEFV_scaffold_39178_1770-2375___|_ORF_|B)91/71)96.7/66,PoEFV_scaffold_31167_4430-4083_-_|_ORF_|B)34.5/46,PoFV_scaffold_10252_orf006))99.2/97;",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster131.aa.aln  --prefix Cluster131.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster131_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster131.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster131.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster131_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster131.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster131_all.treefile -n 0 -zb 1000 -au --prefix Cluster131.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 


########
#Cluster297
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster297_constrain.nwk","w") as output:
  print("(DFV_QKN22505.1,(PoFV_scaffold_7638_orf006,(PoEFV_scaffold_8315_10320-9370_-_|_ORF_|B,(PoEFV_scaffold_18363_5057-5305___|_noORF_|B,PoEFV_scaffold_75772_601-2_-_|_noORF_|B))));",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster297.aa.aln  --prefix Cluster297.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster297_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster297.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster297.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster297_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster297.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster297_all.treefile -n 0 -zb 1000 -au --prefix Cluster297.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 

########
#Cluster444
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster444_constrain.nwk","w") as output:
  print("(PoFV_scaffold_7638_orf010,(PoEFV_scaffold_8315_6633-7148___|_ORF_|B,(PoEFV_scaffold_11969_699-178_-_|_ORF_|B,PoEFV_scaffold_44098_2959-2447_-_|_ORF_|B)));",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster444.aa.aln  --prefix Cluster444.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster444_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster444.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster444.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster444_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster444.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster444_all.treefile -n 0 -zb 1000 -au --prefix Cluster444.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 


########
#Cluster454
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster454_constrain.nwk","w") as output:
  print("(B2YG10_MdSGHV,(((((PoEFV_scaffold_6401_20213-19599_-_|_ORF_|A,PoEFV_scaffold_64707_1543-704_-_|_ORF_|B)94.6/89,PoFV_scaffold_49244_orf001)94.7/53,((PoFV_scaffold_21671_orf001_putative_ATPase__CDC48__LbFV_ORF81,PoEFV_scaffold_33452_1_2005-2787___|_noORF_|B)100/72,LbFV_YP_009345685.1)14.6/32)94.4/95,AKY03092_AmFV)49.9/63,B0YLR2_GpSGHV)58.4/52,B0YLR1_GpSGHV);",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster454.aa.aln  --prefix Cluster454.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster454_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster454.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster454.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster454_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster454.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster454_all.treefile -n 0 -zb 1000 -au --prefix Cluster454.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 


########
#Cluster477
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster477_constrain.nwk","w") as output:
  print("(LbFV_YP_009345682.1,(((PoEFV_scaffold_6401_16593-18479___|_Pseudogenized_|A,PoEFV_scaffold_7180_4204-1607_-_|_ORF_|B)100/100,PoFV_scaffold_15158_orf005)69/69,((B0YLI7_GpSGHV,B2YG51_MdSGHV)85.5/92,B0YLI6_GpSGHV:0.4899820241)100/100:0.9824638465)99.6/100:0.5452482146,DFV_QKN22481.1);",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster477.aa.aln  --prefix Cluster477.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster477_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster477.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster477.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster477_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster477.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster477_all.treefile -n 0 -zb 1000 -au --prefix Cluster477.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 

########
#Cluster536
########

with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster536_constrain.nwk","w") as output:
  print("(PoFV_scaffold_2702_orf005,(PoEFV_scaffold_7821_14388-15077___|_noORF_|B,(PoEFV_scaffold_983_14655-15641___|_ORF_|B,PoEFV_scaffold_17380_2-487___|_ORF_|B)0/59));",file=output)

subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster536.aa.aln  --prefix Cluster536.constr1 -m MFP -alrt 1000  -bb 1000 -bnni  -nt 15 -g /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster536_constrain.nwk" , shell=True)

#Concatenate the 2 trees 
subprocess.run("cat /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster536.constr1.treefile /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster536.aa.aln.treefile >> /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster536_all.treefile", shell=True)

#Test the set of trees 
subprocess.run("/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2  -s /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster536.aa.aln -m MFP -z /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Cluster536_all.treefile -n 0 -zb 1000 -au --prefix Cluster536.constr1_test ", shell=True)

#Look at the resulting .iqtree file:
# Resut: no difference 



####################################################
####### Calculate some statistical observations ####
####################################################

# Calculate AT content of scaffolds containing the ORFs 

PoEFV_scaffolds= "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_scandidate_scaffolds.aa"
Scaffold_assembly_record = SeqIO.to_dict(SeqIO.parse(PoEFV_scaffolds, "fasta"))

AT_scaffolds_content=[]
for scaff in PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['scaf_name'].unique():
  AT_content=((Scaffold_assembly_record[scaff].seq.upper().count('A')+Scaffold_assembly_record[scaff].seq.upper().count('T'))/ (Scaffold_assembly_record[scaff].seq.upper().count('A')+Scaffold_assembly_record[scaff].seq.upper().count('T')+Scaffold_assembly_record[scaff].seq.upper().count('G')+Scaffold_assembly_record[scaff].seq.upper().count('C')))
  AT_scaffolds_content.append(AT_content)

import numpy as np
np.mean(AT_scaffolds_content)

# Calculate CD% of scaffolds containing the ORFs 

CD_scaffolds_perc=[]
for scaff in PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['scaf_name'].unique():
  scaff_length= len(str(Scaffold_assembly_record[scaff].seq))
  Tot_loci_length=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['scaf_name'].eq(scaff)]['loci_length'].sum()
  CD_perc=Tot_loci_length/scaff_length
  CD_scaffolds_perc.append(CD_perc)

np.mean(CD_scaffolds_perc)



#####################
### dN/dS analysis ##
#####################

Cluster_to_run_dNdS=['Cluster149','Cluster183','Cluster207','Cluster301','Cluster379','Cluster406','Cluster407',
'Cluster417','Cluster419','Cluster427','Cluster445','Cluster670','Cluster697','Cluster922','Cluster1183','Cluster1245',
'Cluster1248','Cluster1263','Cluster1417','Cluster1','Cluster9','Cluster63','Cluster85','Cluster101',
'Cluster108', 'Cluster131','Cluster297','Cluster444', 'Cluster454', 'Cluster536','Cluster56','Cluster148',"Cluster69","Cluster296",'Cluster339','Cluster355','Cluster477']



#Do a codon alignment 

for cluster in Cluster_to_run_dNdS:
  subprocess.run("java -jar /beegfs/data/bguinet/TOOLS/macse_v2.05.jar -prog alignSequences -seq /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+cluster+".dna" , shell=True)

# Trim codon alignmnent
for cluster in Cluster_to_run_dNdS:
  subprocess.run("/beegfs/data/bguinet/TOOLS/trimal/source/readal -in /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+cluster+"_NT.dna -out /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+cluster+"_NT.dna.unaligned -onlyseqs ", shell=True)
  subprocess.run("/beegfs/data/bguinet/TOOLS/trimal/source/trimal -backtrans /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+cluster+"_NT.dna.unaligned -in /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+cluster+"_AA.dna -automated1 -resoverlap 0.30 -seqoverlap 30 -fasta -ignorestopcodon -out /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+cluster+"_NT.dna.trimmed ", shell=True)


#################################
#### Pairwise dN/dS analysis ####
#################################

# For each Cluster with only 3 sequences, create a phylogeny manually and create the codeml_file 

Cluster_with_2_paralogs=['Cluster149','Cluster183','Cluster207','Cluster301','Cluster379','Cluster406','Cluster407',
'Cluster417','Cluster419','Cluster427','Cluster445','Cluster670','Cluster697','Cluster922','Cluster1183','Cluster1245',
'Cluster1248','Cluster1263','Cluster1417']

for Clusters in Cluster_with_2_paralogs :
  with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+".aa.aln.treefile","w") as output:
    subtab=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['Cluster'].eq(Clusters)]
    list_virus_loci=[]
    list_EVE_loci=[]
    for loci in subtab['query']:
      list_virus_loci.append(loci)
    for loci in subtab['target']:
      loci=re.sub("\\(\\+\\)","_+_",loci)
      loci=re.sub("\\(-\\)","_-_",loci)
      loci=re.sub(":","_",loci)
      loci=re.sub("PoFV","PoEFV",loci)
      list_EVE_loci.append(loci)
    list_virus_loci=list(dict.fromkeys(list_virus_loci))
    list_EVE_loci=list(dict.fromkeys(list_EVE_loci))
    print("("+ list_virus_loci[0]+ ",("+ list_EVE_loci[0]+","+list_EVE_loci[1]+"));",sep="",file=output)

for Clusters in Cluster_with_2_paralogs :
  with open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Pairwise_dNdS/Codeml_job_"+Clusters+".cmd","w") as output:
    print("runmode = -3",file=output)
    print("CodonFreq = 4",file=output)
    print("model = 0 ",file=output)
    print("NSsites = 0",file=output)
    print("kappa = 1",file=output)
    print("omega = 0.5",file=output)
    print("cleandata = 0",file=output)
    print("noisy = 9",file=output)
    print("verbose = 1",file=output)
    print("seqtype = 1",file=output)
    print("icode = 0",file=output)
    print("fix_kappa = 1",file=output)
    print("fix_kappa = 1",file=output)
    print("seqfile =  /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed",file=output)
    print("outfile = /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Pairwise_dNdS/Pairwise_"+Clusters+"_results.txt",file=output)
    print("treefile = /beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+".aa.aln.treefile",file=output)

#

#Run bayesian pairwise analysis 
import os
path = "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Pairwise_dNdS/"
os.chdir(path)

for Clusters in Cluster_with_2_paralogs :
  subprocess.run("/beegfs/data/bguinet/TOOLS/paml4.9i/bin/codeml Codeml_job_"+Clusters+".cmd",shell=True)

# Extract the pairwise results 
 
from ete3 import EvolTree
import pandas as pd 
import os,re
from collections import defaultdict

Pairwise_dNdS_table= pd.DataFrame(columns=['Cluster','SP1', 'SP2', 'SE_dNdS','Pvalue_dNdS','dNdS'])
for filename in os.listdir("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Pairwise_dNdS/"):
        if filename.endswith("_results.txt"):
            print(filename)
            Cluster=re.sub("Pairwise_","",filename)
            Cluster=re.sub("_results.txt","",Cluster)
            if os.path.exists("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Pairwise_dNdS/"+filename): #
                #print(filename)
                #Parse the codeml file in order to get the dS mean between pairwise loci in list_loci :
                status = 0
                genome_dnds = defaultdict(list)
                #The codeml file 
                codeml = open("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Pairwise_dNdS/"+filename, "r") # codeml pairwise ML output here
                take_it ="no"
                for i in codeml.readlines():
                    if i.startswith("pairwise comparison, codon frequencies"):
                           status = 1
                    if status == 1:
                        if i[0].isdigit():
                            line = i.rstrip()
                            line2 = re.sub("\(", "", line)
                            line3 = re.sub("\)", "", line2)
                            spaces = line3.split(" ")
                            first = spaces[1]
                            second = spaces[4]
                            first_split = first.split("..")
                            second_split = second.split("..")
                            g1 = first_split[0]
                            g2 = second_split[0]
                        if take_it =="yes": 
                                #print(i)
                                take_it ="no"
                                line = i.rstrip()
                                line = re.sub("\s+", "\t", line)
                                tabs = line.split("\t")
                                SE_dNdS = tabs[4]
                                Pvalue_dNdS = tabs[7]
                                dNdS = tabs[2]
                                #Note that not all gene-pairs will be printed out. This is because the script filters out all pairs for 
                                #which dS was < 0.01 or > 2. Values < 0.01 indicate that we may not get a reliable estimate of dN/dS, 
                                #since the sequences are so similar. dS values > 9 indicate that the sequences are quite divergent 
                                #and multiple substitutions have likely occured at most sites, so dN/dS estimates will again be compromised.
                                #if float(ds) < 9 and float(ds) > 0.01 and float(dnds) < 10: #Only for busco genes 
                                Pairwise_dNdS_table=Pairwise_dNdS_table.append({"Cluster":Cluster,"SP1":first,"SP2":second,"SE_dNdS":float(SE_dNdS) , "Pvalue_dNdS": float(Pvalue_dNdS),"dNdS": float(dNdS)}, ignore_index=True)
                                Pairwise_dNdS_table.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Pairwise_dNdS/dN_dS_table.tab",sep=";",index=False)
                                #print(dNdS_table)
                        if "E[t" in i:
                            take_it="yes"
#
Pairwise_dNdS_table=Pairwise_dNdS_table.loc[Pairwise_dNdS_table['SP1'].str.contains("\\|") & Pairwise_dNdS_table['SP2'].str.contains("\\|")]


####
#################################
#### Classic dN/dS analysis ####
#################################
from ete3 import EvolTree
from io import StringIO

if os.path.exists("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Classic_dNdS/dN_dS_table.tab") :
  Classical_dNdS_table=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Classic_dNdS/dN_dS_table.tab",sep=";")
else:
  Classical_dNdS_table= pd.DataFrame(columns=['Cluster','SP1', 'SP2', 'SE_dNdS','Pvalue_dNdS','dNdS'])

Cluster_with_good_topology =["Cluster69","Cluster296",'Cluster339','Cluster355','Cluster477']
Cluster_with_corrected_topology=['Cluster1','Cluster9','Cluster85', 'Cluster131','Cluster297','Cluster444', 'Cluster454', 'Cluster536','Cluster63','Cluster101','Cluster108']
Cluster_classical_dNdS=["Cluster69","Cluster296",'Cluster339','Cluster355','Cluster477','Cluster1','Cluster477','Cluster1','Cluster9','Cluster85', 'Cluster131','Cluster297','Cluster444', 'Cluster454', 'Cluster536','Cluster63','Cluster101','Cluster108','Cluster69',"Cluster296","Cluster339",'Cluster355']

Reduced_cluster=['Cluster63','Cluster101','Cluster108']


for Clusters in Cluster_classical_dNdS:
  if len(Classical_dNdS_table.loc[Classical_dNdS_table['Cluster'].eq(Clusters)])==0:
    if Clusters in Cluster_with_corrected_topology:
      Aln="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed"
      Tree="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_constrain.nwk"
      #Modify leaf names 
      filename = Tree
      subprocess.run("cp "+ filename +" "+ filename+"2",shell=True)
      subprocess.call(["sed","-i",r"s@___@_+_@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@|@@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@/@@g",filename+"2"])
      filename = Aln
      subprocess.run("cp "+ filename +" "+ filename+"2",shell=True)
      subprocess.call(["sed","-i",r"s@|@@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@\[@_@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@]@_@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@(CDC48)@_CDC48_@g",filename+"2"])
      Aln="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed2"
      Tree="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_constrain.nwk2"
    if Clusters in Cluster_with_good_topology :
      Aln="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed"
      Tree="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+".aa.aln.treefile"
      #Modify leaf names 
      filename = Tree
      subprocess.run("cp "+ filename +" "+ filename+"2",shell=True)
      subprocess.call(["sed","-i",r"s@___@_+_@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@|@@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@/@@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@___|_@_+__@g",filename+"2"])
      filename = Aln
      subprocess.run("cp "+ filename +" "+ filename+"2",shell=True)
      subprocess.call(["sed","-i",r"s@|@@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@\[@_@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@]@_@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@(CDC48)@_CDC48_@g",filename+"2"])
      subprocess.call(["sed","-i",r"s@!@-@g",filename+"2"])
      Aln="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed2"
      #Aln_record = SeqIO.parse(Aln,"fasta")
      #SeqIO.write(Aln_record, "/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed3", "fasta")
      #os.rename("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed3","/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed3") 
      Tree="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+".aa.aln.treefile2"
    if Clusters in Reduced_cluster:
      Tree="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+".aa.aln.treefile_reduced"
      Aln="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/"+Clusters+"_NT.dna.trimmed_reduced"
    Output_path="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Classic_dNdS/"
    list_of_species_to_test=[]
    subtab=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk['Cluster'].eq(Clusters)]
    for index, row in subtab.iterrows():
      Scaffold_score=row['Scaffold_score']
      new_loci_name=row['target']
      if row['ORF_perc'] > 0.5:
        if row['Pseudogenized'] == "no":
          new_loci_name=new_loci_name+"_ORF_"+Scaffold_score
      if row['Pseudogenized'] == "yes":
        new_loci_name=new_loci_name+"_Pseudogenized_"+Scaffold_score
      if pd.isna(row['ORF_perc']):
          new_loci_name=new_loci_name+"_noORF_"+Scaffold_score
      new_loci_name=re.sub("\\(\\+\\)","_+_",new_loci_name)
      new_loci_name=re.sub("\\(-\\)","_-_",new_loci_name)
      new_loci_name=re.sub(":","_",new_loci_name)
      new_loci_name=re.sub("PoFV","PoEFV",new_loci_name)
      list_of_species_to_test.append(new_loci_name)
      list_of_species_to_test=list(dict.fromkeys(list_of_species_to_test))
    #list_of_species_to_test=["PoEFV_scaffold_10234_3639-4663_+__","PoEFV_scaffold_84473_912-1202_+__ORF_B"]
    tree = EvolTree(Tree,binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")
    R = tree.get_midpoint_outgroup()
    # and set it as tree outgroup
    tree.set_outgroup(R)
    new_tree=tree.write (format=0,outfile=Tree)
    tree = EvolTree(Tree,binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")
    print("\n")
    print("Here is the tree that will be tested for Cluster : ", Clusters)
    print(tree)
    print("\n")
    tree.workdir = Output_path
    tree.link_to_alignment(Aln)
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
     print(df_marks.shape[0]) 
     if df_marks.shape[0] < 5:
       df_marks=df_marks.loc[df_marks["Node_name"].str.contains("__")]
     else:
       df_marks = df_marks.iloc[1:]
    marks=list(df_marks['Node_id'].astype(int))
    tree.mark_tree(marks, ['#1'])
    print(tree.write ())
    # run neutral model 
    tree.run_model('b_neut'+'.'+str(Clusters),CodonFreq=4,estFreq=1,getSE=1,noisy=3,verbose=1)
    print("neutral model analysis done.")
    print("\n")
    #run free model 
    tree.run_model('b_free'+'.'+str(Clusters),CodonFreq=4,estFreq=1,getSE=1,noisy=3,verbose=1)
    print("free-ratio model analysis done.")
    print("\n")
    #Open the results 
    tree.link_to_evol_model(Output_path+'/b_neut'+'.'+str(Clusters)+'/out', 'b_neut')
    for model in  tree._models:
                    fb_model_neut=tree.get_evol_model(model)
    tree.link_to_evol_model(Output_path+'/b_free'+'.'+str(Clusters)+'/out', 'b_free')
    for model in  tree._models:
                    fb_model_free=tree.get_evol_model(model)              
    #Charging the free model output
    df_free = pd.read_csv(StringIO(fb_model_free.__str__()), skiprows=6,names=['marks','omega','node_ids','name'])
    df_free = df_free.applymap(lambda x: x.split(":")[1])
    #Remove white space
    df_free=df_free.applymap(str.strip).rename(columns=str.strip)
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
    print("dN/dS average is :",omega[0])
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
    with open(Output_path+'/b_free'+'.'+str(Clusters)+'/out', "r") as ifile:
      for line in ifile:
        if line.startswith("SEs for parameters:"):
            SE=next(ifile, ' ').strip()
            SE=re.split('\s+', SE)
            SE=SE[-1]
            #Insert the value into a new dataframe 
            Classical_dNdS_table=Classical_dNdS_table.append({"Cluster":Clusters,"dNdS":omega[0],"Pvalue_dNdS":pvalue,"SE_dNdS":SE,"SP1":"na","SP2":"na"},ignore_index=True)
    print("Cluster : ", Clusters , " dNdS_analysis_done!")
    Classical_dNdS_table.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Clusters/Classic_dNdS/dN_dS_table.tab",sep=";",index=False)

###

# Merge both analysis 
Pairwise_dNdS_table['dNdS_type']="pairwise"
Classical_dNdS_table['dNdS_type']="classic"
dNdS_table=Pairwise_dNdS_table.append(Classical_dNdS_table)


dNdS_table['SE_plus_dNdS']=dNdS_table['dNdS'].astype(float)+ dNdS_table['SE_dNdS'].astype(float)
del dNdS_table['SP1']
del dNdS_table['SP2']
# Add to the main table 

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk.merge(dNdS_table,on="Cluster",how="outer")

#########
# save8 #
#########

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS.loc[PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS['Scaffold_score'].eq("B"),"Scaffold_score"]="C"

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS.loc[~(PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS['scaf_name'].eq("scaffold_11747"))] # False positive 

PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat_euk_dNdS.m8",sep=";",index=False)


PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS= pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat_euk_dNdS.m8",sep=";")



######################
# Some informations about the EVEs #
######################
# Number clusters with PoEFV loci: 74 
# Number clusters with PoEFV loci paralogs : 37 
# Total number of PoEFV loci : 141 
# Total number of complet ORFs PoEFV loci : 92
# Total number of incomplet ORFs PoEFV loci : 42
# Total number of pseudogenized (with stop codon) PoEFV loci : 32
# Total number of pseudogenized and a complete ORF PoEFV loci : 10

## PoFV genome informations 
# Number of genes within PoFV genome : 128 
# PoFV % genome content AT 71.8%
# PoFV CD% genome content : 82.7% 
## PoEFV genome informations 
# PoEFV % genome content AT 70.2%
# PoEFV CD% genome content : 91.4% 

## dN/dS analysis informations 
# Number of clusters with possible dN/dS calculation: 35/74
  # Number of clusters with possible pairwise dN/dS calculation : 19/35 
  # Number of clusters with possible PAML dN/dS calculation : 16/35 
    # Number cluster with constrained topo : 11/16
    # Number cluster with non-constrained topo : 5/16
    # Number of loci with dN/dS significantly < 1: 29/35 
    

summary_tab=PoFV_vs_filamentous_ORFs_cov_HSPs_clusters_ORFs_repeat_euk_dNdS

summary_tab=summary_tab.merge(summary_tab.drop_duplicates(subset = ["target"]).groupby(['Cluster']).size().reset_index(name='Nb_loci'),on="Cluster")

#for index, row in summary_tab.iterrows():

summary_tab['Cluster_type']=np.nan
grouped = summary_tab.groupby(['Cluster'])

for group_name, group in grouped:
    type=np.nan
    original_group_length=len(group['target'].unique())
    Pseudogenized_length=len(group.loc[(group['Pseudogenized']=="yes") | (group['Pseudogenized']=="yes") & (group['ORF_perc'].ge(0.5))]['target'].unique())
    Complet_length=len(group.loc[(group['Pseudogenized']=="no") & (group['ORF_perc'].ge(0.5))]['target'].unique())
    if original_group_length == Pseudogenized_length:
      type="only_pseudo"
    elif original_group_length == Complet_length:
      type="only_complet"
    else:
      type="both"
    summary_tab.loc[summary_tab['Cluster'].eq(group_name),"Cluster_type"]=type
    


# COunt number of paralogues within ORFS


temp = summary_tab.drop_duplicates(subset = "Cluster")
temp=temp.loc[ temp['target'].str.contains("\\(")]
temp.sort_values(by=["query"], inplace = True)
#temp=temp.loc[temp['Nb_loci'].ge(2)]
#Distribution of number of paralogs per PoFV ORFs
temp.Nb_loci.value_counts()

summary_tab.loc[summary_tab['Cluster_type'].eq("only_pseudo") & summary_tab['Pvalue_dNdS'].lt(0.05),"Pvalue_dNdS"]=10
for i in ["A"]:
  print("# Number clusters with PoEFV loci: ", len(summary_tab.loc[summary_tab['target'].str.contains("\\(")]['Cluster'].unique())) 
  print("# Number clusters with PoEFV loci paralogs : ", len(summary_tab.loc[summary_tab['target'].str.contains("\\(") & summary_tab['Nb_loci'].ge(2) ]['Cluster'].unique()))  # 37
  print(" Number of loci corresponding to paralogs : ", len( summary_tab.loc[summary_tab['Nb_loci'].ge(2) & summary_tab['target'].str.contains("\\(")]['target'].unique())) # 104
  print("# Total number of PoEFV loci : ",len(summary_tab['target'].unique()) ) 
  print("# Total number of complet ORFs PoEFV loci : ", len(summary_tab.loc[(summary_tab['Pseudogenized']=="no") & (summary_tab['ORF_perc'].ge(0.5))]['target'].unique())) 
  print("# Total number of pseudogenized (with stop codon) PoEFV loci : ", len(summary_tab.loc[summary_tab['Pseudogenized']=="yes"]['target'].unique() )) 
  print("# Total number of pseudogenized and a complete ORF PoEFV loci : ", len(summary_tab.loc[(summary_tab['Pseudogenized']=="yes") & (summary_tab['ORF_perc'].ge(0.5))]['target'].unique())) 
  print("# Total number of clusters with only pseudogenized loci : ", len(summary_tab.loc[summary_tab['Cluster_type'].eq("only_pseudo")]['Cluster'].unique()))
  print("# Total number of clusters with only complete loci : ", len(summary_tab.loc[summary_tab['Cluster_type'].eq("only_complet")]['Cluster'].unique()))
  print("# Total number of clusters with both pseudogenized and complete loci : ", len(summary_tab.loc[summary_tab['Cluster_type'].eq("both")]['Cluster'].unique()))
  print("# Average percentage of identity : ", summary_tab.loc[summary_tab['query'].str.contains("PoFV")]['pident'].mean())
  #
  print("\n")
  print("## PoFV genome informations") 
  print("# Number of genes within PoFV genome : 128")
  print("# PoFV % genome content AT : 71.8%")
  print("# PoFV CD% genome content : 82.7%")
  print("## PoEFV genome informations")
  print("# PoEFV % genome content AT : ",np.mean(AT_scaffolds_content)*100,"%")
  print("# PoEFV CD% genome content : ",np.mean(CD_scaffolds_perc)*100,"%")
  #
  print("\n")
  print("## dN/dS analysis informations")
  print("# Number of clusters with possible dN/dS calculation : ",len(Cluster_with_2_paralogs)+len(Cluster_with_corrected_topology)+ len(Cluster_with_good_topology))
  print("# Number of loci with possible dN/dS calculation ; ", len(summary_tab.loc[summary_tab['Cluster'].isin(Cluster_with_2_paralogs+Cluster_with_corrected_topology+Cluster_with_good_topology)]['target'].unique()))
  print("  # Number of clusters with possible pairwise dN/dS calculation : ", len(Cluster_with_2_paralogs))
  print("  # Number of clusters with possible PAML dN/dS calculation : ", len(Cluster_with_corrected_topology)+ len(Cluster_with_good_topology))
  print("    # Number cluster with constrained topo :", len(Cluster_with_corrected_topology))
  print("    # Number cluster with non-constrained topo : ", len(Cluster_with_good_topology))
  print("    # Number of loci with dN/dS significantly < 1: ", len(summary_tab.loc[summary_tab['SE_plus_dNdS'].lt(1) & summary_tab['Pvalue_dNdS'].lt(0.03)]['Cluster'].unique()))




##################################################
### Code to display the loci in a plot within R ##
##################################################
import pandas as pd 
import numpy as np
import os,re 
from Bio import SeqIO 

PoEFV_tab=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/Filamentous_ORFs_vs_Porseoliae_result_cov_HSPs_NR_clusters_ORFs_repeat_euk_dNdS.m8",sep=";")

PoEFV_tab=PoEFV_tab.loc[PoEFV_tab['query'].str.contains("PoFV")]

PoEFV_tab['Complete']=np.nan

PoEFV_tab.loc[PoEFV_tab['ORF_perc'].ge(0.5) & PoEFV_tab['Pseudogenized'].eq("no"),"Complete"]="yes"
PoEFV_tab.loc[PoEFV_tab['Pseudogenized'].eq("yes"),"Complete"]="no"
PoEFV_tab.loc[PoEFV_tab['ORF_perc'].lt(0.5),"Complete"]="no"
PoEFV_tab.loc[PoEFV_tab['ORF_perc'].isna(),"Complete"]="no"

PoEFV_tab['Purifying_selection']="no"
PoEFV_tab.loc[summary_tab['SE_plus_dNdS'].lt(1) & PoEFV_tab['Pvalue_dNdS'].lt(0.03),"Purifying_selection"]="yes"

Porseoliae_assembly="/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/Platygaster_orseoliae_corrected2.fa"
Porseoliae_assembly=SeqIO.to_dict(SeqIO.parse(Porseoliae_assembly, "fasta"))

PoEFV_tab['scaf_length']=np.nan

for index, row in PoEFV_tab.drop_duplicates(subset = ["scaf_name"]).iterrows():
  PoEFV_tab.loc[PoEFV_tab['scaf_name'].eq(row['scaf_name']),"scaf_length"]=len(str(Porseoliae_assembly[row['scaf_name']].seq))

PoEFV_tab=PoEFV_tab.drop_duplicates(subset = ["target"])

PoEFV_tab.sort_values(by=["scaf_name","start"], inplace = True)

PoEFV_tab['PoFV_scaf_name']=PoEFV_tab['query'].str.replace("_orf.*","")
PoEFV_tab=PoEFV_tab[['target','scaf_name','start','end','strand','Pseudogenized','ORF_perc','query','Complete','Purifying_selection','scaf_length','Scaffold_score']]

PoEFV_tab['PoFV_ORFS']=PoEFV_tab['query'].str.replace("PoFV_","")
PoEFV_tab['scaf_name2']=PoEFV_tab['scaf_name']+' - '+PoEFV_tab['Scaffold_score']

# Create the tables 

PoFV_tab=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Genomes/PoFV/Predicted_orfs/Final_ORF_prediction_PoFV.gff",sep="\t",comment="#",header=None)

PoFV_tab.columns=['scaf_name','method','type','start','end','point1','strand','point2','informations']
PoFV_tab['PoFV_ORFS']=PoFV_tab['informations'].str.replace("Name=","")
PoFV_tab['PoFV_ORFS']=PoFV_tab['PoFV_ORFS'].str.replace(";.*","")
PoFV_tab['PoFV_ORFS']=PoFV_tab['PoFV_ORFS'].str.replace("PoFV_","")
PoFV_tab['PoFV_ORFS']=PoFV_tab['PoFV_ORFS'].str.replace(" .*","")
PoFV_tab['PoFV_ORFS']=PoFV_tab['PoFV_ORFS'].str.replace("_putative.*","")

# Add the ORF names
PoFV_gene_name_tab=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoFV_gene_names.txt",sep=";")
PoFV_gene_name_tab['PoFV_ORFS']=PoFV_gene_name_tab['PoFV_ORFS'].str.replace(" ","")

PoFV_tab=PoFV_tab.merge(PoFV_gene_name_tab,on="PoFV_ORFS",how="outer")
PoFV_tab.loc[PoFV_tab['Gene_names'].isna(),"Gene_names"]="Unknown"

for index, row in PoFV_tab.drop_duplicates(subset = ["scaf_name"]).iterrows():
  PoFV_tab.loc[PoFV_tab['scaf_name'].eq(row['scaf_name']),"scaf_length"]=len(str(Porseoliae_assembly[row['scaf_name']].seq))

#
count_PoFV_tab=PoFV_tab.drop_duplicates(subset = "PoFV_ORFS")
count_PoFV_tab.sort_values(by=["PoFV_ORFS"], inplace = True)
count_PoFV_tab['ORF_number'] = np.arange(count_PoFV_tab.shape[0])
count_PoFV_tab['ORF_number'] =count_PoFV_tab['ORF_number'] +1

# Number the ORFs 
PoFV_tab=PoFV_tab.merge(count_PoFV_tab[['PoFV_ORFS','ORF_number']],on="PoFV_ORFS")
PoFV_tab.sort_values(by=["PoFV_ORFS"], inplace = True)

PoFV_tab.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoFV_table_for_GGGenes_plot.txt",sep=";")

PoEFV_tab=PoEFV_tab.merge(count_PoFV_tab[['PoFV_ORFS','ORF_number']],on="PoFV_ORFS")
PoEFV_tab=PoEFV_tab.merge(PoFV_gene_name_tab,on="PoFV_ORFS",how="outer")
PoEFV_tab.loc[PoEFV_tab['Gene_names'].isna(),"Gene_names"]="Unknown"


PoEFV_tab.to_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/PoEFV_analysis/PoEFV_table_for_GGGenes_plot.txt",sep=";")




    
