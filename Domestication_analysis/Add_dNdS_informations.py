from multipy.fdr import tst
import os 
import pandas as pd
import argparse
import numpy as np
# Print out a message when the program is initiated.
print('--------------------------------------------------------------------------\n')
print('                        Add dN/dS informations and calculate FDR pvalues .\n')
print('--------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow to detect and contact Hsp sequences within an alignment')
parser.add_argument("-b", "--blast", help="The blast file with all informations")
parser.add_argument("-d", "--dNdS_directory",help="The directory where to find all the .out dN/dS results")
parser.add_argument("-o", "--Output_file",help="The desired output file name")
args = parser.parse_args()

# Variable that stores fasta sequences
blast_table=args.blast
Output_file=args.Output_file
dNdS_directory=args.dNdS_directory

#Example usage 
#python3 Add_dNdS_informations.py -d /beegfs/data/bguinet/these/dNdS_analysis_filtred/dNdS_results/ -b /beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM.m8 -o /beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS.m8

#Load the blast table 
blast_table=pd.read_csv(blast_table,sep=";")

#Add dNdS results 

#Load the dNdS informations
dNdS_full_tab = pd.DataFrame(columns=("Clustername", "Event", "Mean_dNdS",  "Pvalue_dNdS","SE_dNdS"))

for filename in os.listdir(dNdS_directory):
	if filename.endswith(".out"):
		dNdS_tab= pd.read_csv(dNdS_directory+filename, sep="\t",index_col=None)
		dNdS_full_tab=dNdS_full_tab.append(dNdS_tab, ignore_index=True)

Full_tab_with_dNdS = pd.merge(blast_table, dNdS_full_tab, left_on=["Clustername","Event"],right_on=['Clustername',"Event"],how='outer')
Full_tab_with_dNdS = Full_tab_with_dNdS.drop_duplicates()
    
#for i in list_nb:
#	Full_tab_with_dNdS_pvalue_dNdS['FDR_pvalue_dNdS']=tst(Full_tab_with_dNdS_pvalue_dNdS['Pvalue_dNdS'], q=i)
#	print(i,len(Full_tab_with_dNdS_pvalue_dNdS.loc[Full_tab_with_dNdS_pvalue_dNdS['FDR_pvalue_dNdS'].astype(str)=="False"]))
#Here we take 0.05 because no coude 

Full_tab_with_dNdS_pvalue_dNdS= Full_tab_with_dNdS.drop_duplicates(subset=['Clustername','Event'], keep='first')
Full_tab_with_dNdS_pvalue_dNdS=Full_tab_with_dNdS_pvalue_dNdS[Full_tab_with_dNdS_pvalue_dNdS['Pvalue_dNdS'].notna()]
Full_tab_with_dNdS_pvalue_dNdS['FDR_pvalue_dNdS']=tst(Full_tab_with_dNdS_pvalue_dNdS['Pvalue_dNdS'], q=0.05)
Full_tab_with_dNdS_pvalue_dNdS=Full_tab_with_dNdS_pvalue_dNdS[['Event','Clustername','FDR_pvalue_dNdS']]

Full_tab_with_dNdS=Full_tab_with_dNdS.merge(Full_tab_with_dNdS_pvalue_dNdS, left_on=['Clustername','Event'], right_on=['Clustername','Event'],how='outer')


#Add number of loci within the same scaffold 
Full_tab_with_dNdS['Nb_same_loci_on_scaff'] = (Full_tab_with_dNdS.groupby(['Scaff_name']).Clustername.transform('nunique').mask(lambda x: x == 1, 0))

#Fill HSP values 
Full_tab_with_dNdS[['new1','new2']] = np.sort(Full_tab_with_dNdS[['query','HSP_name']].fillna('').values,1)
up = Full_tab_with_dNdS.groupby(['Clustername','new1','new2'])[['Event', 'Nloc', 'Nsp', 'Nsp_MRCA', 'Tips4dNdS.array', 'Boot', 'Nsp_losses', 'Mean_dNdS', 'Pvalue_dNdS', 'SE_dNdS', 'FDR_pvalue_dNdS', 'Nb_same_loci_on_scaff','TPM_all']].apply(lambda x : x.ffill().bfill())
Full_tab_with_dNdS.update(up)
Full_tab_with_dNdS.drop(['new1','new2'],1,inplace=True)


#Fill NaN Event 

s = Full_tab_with_dNdS[Full_tab_with_dNdS['Event'].isnull()].groupby(['Clustername', 'query']).ngroup()

to_fill = (s - s.groupby(Full_tab_with_dNdS['Clustername']).transform('min') + 1
           + Full_tab_with_dNdS.groupby('Clustername')['Event'].transform('max'))

Full_tab_with_dNdS['Event'] = Full_tab_with_dNdS['Event'].fillna(to_fill, downcast='infer')

#Fill Clustername Event with only NaN 

sizes = Full_tab_with_dNdS.groupby(['Clustername']).size()
Full_tab_with_dNdS['Size']=Full_tab_with_dNdS['Clustername'].map(sizes)
m = Full_tab_with_dNdS['Event'].isna()

Full_tab_with_dNdS.loc[m, 'Event_new']  = Full_tab_with_dNdS.loc[m].groupby(['Clustername'])['Size'].transform(lambda x:list(range(1,len(x)+1)))
mapping = Full_tab_with_dNdS.loc[m].groupby(['Clustername', 'query'])['Event_new'].first().copy()
Full_tab_with_dNdS.set_index(['Clustername', 'query'], inplace=True)
m = Full_tab_with_dNdS['Event'].isna()
Full_tab_with_dNdS.loc[m,'Event'] = Full_tab_with_dNdS.loc[m].index.map(mapping)
Full_tab_with_dNdS.reset_index(inplace=True)
Full_tab_with_dNdS.drop(columns=['Size', 'Event_new'], inplace=True)
Full_tab_with_dNdS['Event']=Full_tab_with_dNdS['Event'].astype(int)


Full_tab_with_dNdS.to_csv(Output_file,index=False,sep=";")
print("Output file written to : ", Output_file)

