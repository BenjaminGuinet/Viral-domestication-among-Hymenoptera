import pandas as pd
import argparse
import os
import io 

# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Filter clusters       .\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-c", "--clustering_file", help="The blast file in .m8 format")
parser.add_argument("-o", "--out_file", help="The name of the outpit reciprocal file")
args = parser.parse_args()

Clustering_file= args.clustering_file
out_file=args.out_file
#Clustering_file="/beegfs/data/bguinet/these/Clustering3/Candidate_viral_loci_and_viral_protein_clu.tab"

df = pd.read_csv(Clustering_file,"\t")
df=df.drop(columns=['Unnamed: 0'])
df.columns = ["Clusternames", "Names"]
"""
I keep a cluster if: - There is at least one seqname with + or - in its name AND at least one other seqname.

#I remove if :
#-There is only seqname with + or - in their name. 
#-There is only other seqname.

"""
list1 = df.loc[df.Names.str.contains(":")].Clusternames.unique()
list2 = df.loc[~df.Names.str.contains(":")].Clusternames.unique()
final_list = set(list1).intersection(set(list2))
df=df.loc[df.Clusternames.isin(final_list)]

Nb_clusters=len(df['Clusternames'].unique())
#Remove cluster where there is more than 50 different species 

df['sp_names']=df["Names"].str.rsplit(":", 1).str[-1] 
subdf=df.loc[df['Names'].str.contains(":")]
subdf['sp_count'] = subdf.groupby('Clusternames')['sp_names'].unique()
count_df=subdf.groupby('Clusternames').sp_names.nunique()
count_df=pd.DataFrame({'Clusternames':count_df.index, 'nb_sp':count_df.values})

df['virus_names']=df["Names"].str.rsplit(":", 1).str[-1] 
subdfvirus=df.loc[~df['Names'].str.contains(":")]
subdfvirus['virus_count'] = subdfvirus.groupby('Clusternames')['virus_names'].unique()
count_dfvirus=subdfvirus.groupby('Clusternames').virus_names.nunique()
count_dfvirus=pd.DataFrame({'Clusternames':count_dfvirus.index, 'nb_virus':count_dfvirus.values})

new_df=pd.merge(count_df,df,on='Clusternames',how="outer")
new_df=pd.merge(count_dfvirus,new_df,on='Clusternames',how="outer")

df_huge_Cluster=new_df.loc[new_df['nb_sp'].ge(50) & new_df['nb_virus'].ge(5)]
new_df=new_df.loc[new_df['nb_sp'].lt(50)]

print("Number of cluster removed : ",Nb_clusters-len(new_df['Clusternames'].unique()),"/", Nb_clusters)
print("Number remaining clusters :", len(new_df['Clusternames'].unique()))
print("Number of cluster with more than 50 species :", len(df_huge_Cluster['Clusternames'].unique()))
new_df=new_df[['Clusternames','Names']]

new_df.to_csv(out_file, sep='\t')
df_huge_Cluster.to_csv("/beegfs/data/bguinet/these/Clustering3/Filtred_Candidate_viral_loci_and_viral_protein_clu_huge.tab",sep="\t",index=False)
