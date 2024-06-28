from bioservices import UniProt
import pandas as pd
import sys 
import time
import argparse
import os
import re
# Print out a message when the program is initiated.
print('--------------------------------------------------------------------------------------------------------\n')
print('                        Add genomic ontology informations and Interproscan results into Clustering tab.\n')
print('--------------------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add ontology and Interproscan result')
parser.add_argument("-c", "--cluster_file", help="The file with clusters")
parser.add_argument("-I", "--interproscan_file", help="The Interproscan file in tsv format")
parser.add_argument("-o", "--out_file",help="The ouptut desired file name ")
args = parser.parse_args()

#Get the uniprot IDS

# Usage exemple : python3 /beegfs/home/bguinet/these_scripts_2/Add_interproscan_and_ontology.py -I ALL_Predicted_and_known_ORFs_intrproscan.tsv -c ALL_Predicted_and_known_ORFs_cluster.tab -o ALL_Predicted_and_known_ORFs_cluster_ontology.tab

List_viruses_name=["_AcMNPV","_LdMNPV","_CpV","_NeseNPV","_CuniNPV","_HzNV-1","_HzNV-2","_GbNV","_OrNV","_ToNV","_DiNV","_DmNV_kal","_DmNV_tom","_DmNV_esp","_DmNV_mau","_PmNV","_HgNV","_DhNV","_GpSGHV","_MdSGHV","_LbFV","_AmFV"]

start_time = time.time()

cluster= args.cluster_file
out_file=args.out_file
interproscan_file=args.interproscan_file

cluster="ALL_Predicted_and_known_ORFs_cluster.tab"
u = UniProt(verbose=False)
cluster=pd.read_csv(cluster,header=0,sep="\t")

cluster2=cluster.loc[~(cluster['Names'].str.contains("_\\+_") | ~(~cluster['Names'].str.contains("_\\-_")))]

cluster2['Names2']= cluster2['Names']
cluster2['Names2'] = cluster2['Names2'].str.replace('|'.join(List_viruses_name), '',regex=True)

out_dir=file  = os.path.dirname(os.path.realpath(out_file))

count_row = cluster2.shape[0]
filecount = 1

with open(out_dir+"/ID_uniprot.txt", "a") as file:
	for query in cluster2['Names2']:
		print(query)
		if  '_CP' in query:
			query=re.sub('_CP.*',"",query)
		df2=u.search(query,frmt="tab",columns=("lineage-id,id,organism,genes,protein names,go, go(biological process), go(molecular function),go(cellular component), go-id,reviewed"))
		try:
			if len(df2) != 0:
				df2=df2.replace('Taxonomic lineage IDs','query\tTaxonomic lineage IDs')#in order to add a column query to merge the two informations after
				df2=df2.replace('\n','\n'+query+'\t',1)
				matches = df2.count('\n')
				if matches > 2 :
					df2="\n".join(df2.split("\n", 2)[:2])
				print(df2,file=file)
			if len(df2) ==0:
				print("no value")
				df2='query\tTaxonomic lineage IDs\tEntry\tOrganism\tGene names\tProtein names\tGene ontology (GO)\tGene ontology (biological process)\tGene ontology (molecular function)\tGene ontology (cellular component)\tGene ontology IDs\tStatus\n'+query+'\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n'
				print(df2,file=file)
			filled_len = int(round(50 * filecount / float(count_row -1)))
			percents = round(100.0 * filecount / float(count_row ), 1)
			bar = '=' * filled_len + '-' * ((50) - filled_len)
			sys.stdout.write('[%s] %s%s Get Uniprot ids...%s\r' % (bar, percents, '%', ''))
			sys.stdout.flush()  # As suggested by Roid', right_index=True)
			filecount += 1
		except:
			print(query, " did not work")
file.close()

Uniprot_ids=pd.read_csv(out_dir+"/ID_uniprot.txt",header=0,sep="\t")
Uniprot_ids.drop(Uniprot_ids[Uniprot_ids['query'] == 'query'].index, inplace= True) 

#Allow to get a proper dataframe with only one row for colnames. 
cluster['Names2']= cluster['Names']
cluster['Names2'] = cluster['Names2'].str.replace('|'.join(List_viruses_name), '',regex=True)

df_merged = pd.merge(cluster, Uniprot_ids, left_on=['Names2'],right_on=['query'],how='outer')

df_merged=df_merged.drop(['Unnamed: 0'],axis=1)
df_merged=df_merged.drop_duplicates(keep='first')


#Transform and add InterproScan informations 

interproscan_file="ALL_Predicted_and_known_ORFs_intrproscan.tsv"
Interproscan=pd.read_csv(interproscan_file,sep="\t")

Interproscan.columns=['Prot_Acc','MD5','Seq_length','Analysis','Signature_acc','Signature_description','Start','End','Score','Status','Date','Interpro_acc','Interpro_description','GO']

Interproscan=Interproscan[['Prot_Acc','Analysis','Signature_acc','Signature_description','Score','Status','Interpro_acc','Interpro_description','GO']]

Interproscan=Interproscan.drop_duplicates(subset=['Prot_Acc','Analysis'], keep="first")

#Stack and unstack then collapse multiindex columns into a one level column
s=Interproscan.set_index(['Prot_Acc','Analysis']).stack().unstack('Prot_Acc').T
s.columns = [f'{a}_{b}' for a, b in s.columns]

df_merged_interpro = pd.merge(df_merged, s, left_on=['Names'],right_on=['Prot_Acc'],how='outer')

df_merged_interpro.to_csv(out_file, sep=';',na_rep='NaN', index=False)

print("\n")
print("Process done for ", count_row ,"sequences")
print("Output written at: ",out_file)
print("--- %s seconds ---" % (time.time() - start_time))

