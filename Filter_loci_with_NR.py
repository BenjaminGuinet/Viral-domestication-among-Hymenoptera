
import pandas as pd 
import argparse
import os
import io 

# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Filter loci with NR results                               .\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-l", "--blast_file", help="The blast file in .m8 format")
parser.add_argument("-o", "--out_file", help="The name of the output filtred file")
parser.add_argument("-Nr", "Nr_results", help="The Nr table result file")
args = parser.parse_args()

#Usage example Filter_loci_with_NR.py -b /beegfs/data/bguinet/these/Results/Candidate_viral_loci.m8 -Nr /beegfs/data/bguinet/these/Results/Candidate_viral_loci_Nr.m8 -o /beegfs/data/bguinet/these/Results/Candidate_viral_loci_filtred.m8
Blast_file= args.loci_file
out_file=args.out_file

blast = pd.read_csv(Blast_file,";")

tab_Nr= pd.read_csv("/beegfs/data/bguinet/these/Clustering4/Candidate_loci_filtred_end_reciprocal_result_filtred_taxid.m8",sep=";")

#If only viruses hits keep all clusters 
tab_Nr = tab_Nr.loc[~(tab_Nr['evalue'].gt(0.0000001) & tab_Nr['superkingdom'].isin(['Eukaryota', 'Archaea', 'Bacteria']))]

#Sort evalues 
tab_Nr=tab_Nr.sort_values(by='evalue', ascending=True)
NaN_tab_Nr = tab_Nr.loc[tab_Nr['species'].isna()]
tab_Nr = tab_Nr.drop_duplicates(subset=['query', 'species'], keep='first')
tab_Nr = tab_Nr.append(NaN_tab_Nr)
tab_Nr.superkingdom = tab_Nr.superkingdom.fillna('Unknown')

#Create table with Nr match informations 
#Take into account only non-hymenoptera matches 
tab_Nr['count_Hymenoptera'] = tab_Nr['order'].eq('Hymenoptera').groupby(tab_Nr['query']).transform('sum')
tab_Nr= tab_Nr.loc[~tab_Nr['order'].str.contains("Hymenoptera",na=False)]
tab_Nr = tab_Nr.join(pd.crosstab(tab_Nr['query'], tab_Nr['superkingdom']).add_prefix('Count_'), on='query')
tab_Nr['Count_Bacteria_Eukaryota_Archeae'] = tab_Nr['Count_Eukaryota'] + tab_Nr['Count_Bacteria']+ tab_Nr['Count_Archaea']
tab_Nr['Dif_Bacteria_Eukaryota_Archeae_Viruses'] = tab_Nr['Count_Bacteria_Eukaryota_Archeae'] - tab_Nr['Count_Viruses']

#Remove loci if to much bacteria hits compared to viruses (>20 +) and if less then 30 viruses hits 
tab_Nr_remove1=tab_Nr.loc[tab_Nr['Count_Bacteria_Eukaryota_Archeae'].ge(10)]
tab_Nr_remove1=tab_Nr_remove1.loc[~tab_Nr_remove1['Dif_Bacteria_Eukaryota_Archeae_Viruses'].lt(20)]
tab_Nr_remove1=tab_Nr_remove1.loc[~tab_Nr_remove1['Count_Viruses'].gt(30)]

#If less then 3 virus hit and more then 5 euk hit remove 
tab_Nr_remove2=tab_Nr.loc[tab_Nr['Count_Viruses'].lt(3) & tab_Nr['Count_Bacteria_Eukaryota_Archeae'].ge(5)]
tab_Nr_remove2=tab_Nr_remove2.loc[~tab_Nr_remove2['Count_Bacteria_Eukaryota_Archeae'].lt(15)]

#If more than 100 bacteria hits, remove 
tab_Nr_remove3=tab_Nr.loc[tab_Nr['Count_Bacteria'].gt(100)]

#Create a table 'tab_Nr_remove' composed of false positves loci 
tab_Nr_remove=tab_Nr_remove1.append(tab_Nr_remove2)
tab_Nr_remove=tab_Nr_remove.append(tab_Nr_remove3)

#Keep unknown loci 
tab_unknown = tab_Nr.loc[~tab_Nr['query'].isin(tab_Nr_remove['query'])]
tab_unknown = tab_unknown.drop_duplicates(subset = 'query', keep = 'first')

list_unknown= list(tab_unknown.loc[tab_unknown['Count_Viruses'].lt(3) & tab_unknown['Count_Unknown'].gt(1)]['query'].unique())

#Keep one representing loci to remove 
tab_Nr= tab_Nr.drop_duplicates(subset = 'query', keep = 'first')

#Remove from the original blast file 
blast =blast.loc[~blast['query'].isin(tab_Nr['query'])]

blast.to_csv(output_file,sep=";",indexe=false)


