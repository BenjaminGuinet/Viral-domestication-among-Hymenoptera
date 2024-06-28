###Get uniprot IDs#######

from bioservices import UniProt
import pandas as pd
import sys 
import time
import argparse
import os
import re
# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Add genomic ontology informations into MMseqs2 search tab.\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-b", "--blast_file", help="The blast file in .m8 format (check the column names and modify the colname 'target' wich corresponds to the Acc_number column in the mmseqs2 software")
parser.add_argument("-o", "--out_file",help="The ouptut path where to create the oufile")
args = parser.parse_args()

#Get the uniprot IDS

# Usage exemple python3 /beegfs/home/bguinet/these_scripts_2/Add_genomic_ontology.py -b /beegfs/data/bguinet/these/Clustering/mmseqs2_search_analysis/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid.m8 -o /beegfs/data/bguinet/these/Clustering/mmseqs2_search_analysis/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology.m8


start_time = time.time()

blast= args.blast_file
out_file=args.out_file

u = UniProt(verbose=False)
blast=pd.read_csv(blast,header=0,sep=";")

#blast.drop(["Unnamed: 0",], axis=1,inplace=True)
blast2=blast.drop_duplicates(subset='target', keep="last")
out_dir=file  = os.path.dirname(os.path.realpath(out_file))

count_row = blast2.shape[0]
filecount = 1

with open(out_dir+"/ID_uniprot.txt", "a") as file:
	for query in blast2['target']:
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
df_merged = pd.merge(blast, Uniprot_ids, left_on=['target'],right_on=['query'],how='outer')
print(list(df_merged))
#df_merged.drop(['query_bis','query_y'], axis=1,inplace=True)
#df_merged.drop(['Unnamed: 0', 'Unnamed:_0', 'Unnamed:_0.1','query_y'],axis=1))
df_merged= df_merged.rename(columns={'query_x': 'query'})
df_merged= df_merged.dropna(subset = ['target'])
df_merged=df_merged.drop_duplicates(keep='first')

print("Number of queries with ontology  assigned : ", len(df_merged[df_merged['Gene ontology IDs'].notna()]['query'].unique()))

try:
	df_merged.drop(['Unnamed: 0'], axis = 1, inplace = True) 
except:
	print("")
df_merged.to_csv(out_file, sep=';',na_rep='NaN', index=False)

print("\n")
print("Process done for ", count_row ,"sequences")
print("Output written at: ",out_file)
print("--- %s seconds ---" % (time.time() - start_time))

