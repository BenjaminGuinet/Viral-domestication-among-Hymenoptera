#!/usr/bin/python
import pandas as pd
import os.path
import time
import argparse



# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Add transcriptomic informations into MMseqs2 search tab.\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-b", "--blast", help="The blast file in .m8 format (check the column names and modify the colname 'target' wich corresponds to the Acc_number column in the mmseqs2 software")
parser.add_argument("-o", "--output_file", help="the output file name")
args = parser.parse_args()



blast_file=args.blast
output_file=args.output_file

#Example usage : python3 /beegfs/home/bguinet/these_scripts_2/Add_transcriptomic.py -b /beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono.m8 -o /beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM.m8
#list_of_names2=["Acromyrmex_echinatior","Anagyrus_pseudococci"]
#blast_file="/beegfs/data/bguinet/these/Clustering/mmseqs2_search_analysis/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid.m8"
#out_file="/beegfs/data/bguinet/these/Clustering/mmseqs2_search_analysis/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_env.m8"

#--------------#
start_time = time.time()
#####
# Part in order to get the correct species names format
####
#num_species = sum(1 for line in open(Species_name_file)) # count number of lines in the file in order to now the number of decimals / by the first number given
#list_of_names1=[]
#for names in open(Species_name_file,"r"):
#	list_of_names1.append(names)
#list_of_names2=[]

#for names in list_of_names1:
#	list_of_names2.append(names.replace("\n", ""))

#count_row = len(list_of_names2)
#filecount = 1



blast_file= pd.read_csv(blast_file,sep=";",header=0)

blast_file['Species_name']=blast_file['Species_name'].str.replace(r'.*__','',)

list_of_names2=blast_file.loc[~blast_file['Species_name'].str.contains("HSPs")]['Species_name'].unique()

#Create an empty dataframe filled with 0
columns=['Gene_Id','Nb_Transcripto_Reads','TPM']
Transcriptomic_df = pd.DataFrame(columns=columns)


#for i in list_to_remove:
#	print(i)
#	list_of_names2.remove(i)
#The idea is to parse all species genomes and to get scaffolds with candidates that have a good probability to be from eucaryotic origin.


for names in list_of_names2:
  if os.path.exists("/beegfs/data/bguinet/these/Genomes/"+str(names)+"/Mapping_transcriptomic/mapping_"+str(names)+"_genes.out"):
    file_exists=1
  else: 
    file_exists=0
  if os.path.exists("/beegfs/data/bguinet/these/Genomes/"+str(names)+"/Mapping_transcriptomic/mapping_"+str(names)+"_ovaries_genes.out") :
    file_ovarie_exists=1
  else:
    file_ovarie_exists=0
  
  if file_exists == 1 and  file_ovarie_exists==1:
    Sub_transcriptomic_file=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+str(names)+"/Mapping_transcriptomic/mapping_"+str(names)+"_genes.out",sep="\t")
    Sub_transcriptomic_file=Sub_transcriptomic_file[['Gene_Id','Reads','TPM']]
    Sub_transcriptomic_file.columns=['Gene_Id','Nb_Transcripto_Reads_whole_body','TPM_whole_body']
    Sub_transcriptomic_ovarie_file=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+str(names)+"/Mapping_transcriptomic/mapping_"+str(names)+"_ovaries_genes.out",sep="\t")
    Sub_transcriptomic_ovarie_file=Sub_transcriptomic_ovarie_file[['Gene_Id','Reads','TPM']]
    Sub_transcriptomic_ovarie_file.columns=['Gene_Id','Nb_Transcripto_Reads_ovarie','TPM_ovarie']
    Sub_transcriptomic_all_file = pd.merge(Sub_transcriptomic_file,Sub_transcriptomic_ovarie_file, how='outer',left_on=['Gene_Id'],right_on=['Gene_Id']).dropna()
    Transcriptomic_df=Transcriptomic_df.append(Sub_transcriptomic_all_file, sort=False)
    Transcriptomic_df['TPM_all']=Transcriptomic_df['TPM_whole_body']+Transcriptomic_df['TPM_ovarie']
  if file_exists == 0 and  file_ovarie_exists==1:
    Sub_transcriptomic_ovarie_file=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+str(names)+"/Mapping_transcriptomic/mapping_"+str(names)+"_ovaries_genes.out",sep="\t")
    Sub_transcriptomic_ovarie_file=Sub_transcriptomic_ovarie_file[['Gene_Id','Reads','TPM']]
    Sub_transcriptomic_ovarie_file.columns=['Gene_Id','Nb_Transcripto_Reads_ovarie','TPM_ovarie']
    Transcriptomic_df=Transcriptomic_df.append(Sub_transcriptomic_ovarie_file, sort=False)
    Transcriptomic_df['TPM_all']=Transcriptomic_df['TPM_ovarie']
  if file_exists == 1 and  file_ovarie_exists==0:
    Sub_transcriptomic_file=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+str(names)+"/Mapping_transcriptomic/mapping_"+str(names)+"_genes.out",sep="\t")
    Sub_transcriptomic_file=Sub_transcriptomic_file[['Gene_Id','Reads','TPM']]
    Sub_transcriptomic_file.columns=['Gene_Id','Nb_Transcripto_Reads_whole_body','TPM_whole_body']
    Transcriptomic_df=Transcriptomic_df.append(Sub_transcriptomic_file, sort=False)
    Transcriptomic_df['TPM_all']=Transcriptomic_df['TPM_whole_body']
  if file_exists == 0 and  file_ovarie_exists==0:
    continue
    
#Then we get the Transcriptomic_df and will have do merge its content to the blats file 

Number_query_start=len(blast_file['query'].unique())
Number_target_start=len(blast_file['target'].unique())
Blast_and_transcripto = pd.merge(Transcriptomic_df,blast_file, left_on=["Gene_Id"],right_on=['query'],how='outer')
Blast_and_transcripto.drop(['Gene_Id'],axis=1,inplace=True)

Blast_and_transcripto = Blast_and_transcripto[Blast_and_transcripto['query'].notna()]

print("Number of query remaining : ", len(Blast_and_transcripto['query'].unique()),'/',Number_query_start)
print("Number of target remaining : ", len(Blast_and_transcripto['target'].unique()),'/',Number_target_start)
Blast_and_transcripto.to_csv(output_file, sep=';',index=False)


print("Output file have been written to :", output_file)
