import pandas as pd 
import os 
import re 
import argparse


# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Filter TEs  results                                       .\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow  to filter TEs research and count number of overlapping TEs within sequences')
parser.add_argument("-b", "--Blast_table", help="The name of the Blast file with all filtred candidate loci")
parser.add_argument("-TE", "TE_dir", help="The TE directory with .tab blast result extension files ")
parser.add_argument("-o", "output_file", help="The output desired file name ")
args = parser.parse_args()


#Example usage : python3 Filter_TEs.py -b /beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS_synteny_ORFs2.m8 -TE /beegfs/data/bguinet/these/Repeat_env_analysis2/ /beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS_synteny_ORFs2_TE.m8

Blast_file= args.loci_file
out_file=args.out_file
Blast_table=args.Blast_table
TE_table_dir=args.TE_dir
Overlapping_tab_output= args.Overlapping_tab_output
seqoverlap=args.seqoverlap

Blast_file="/beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS_synteny_ORFs2.m8"
out_file="/beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS_synteny_ORFs2_TE.m8"
TE_table_dir = "/beegfs/data/bguinet/these/Repeat_env_analysis2/"

seqoverlap=10
###Check overlapping TEs withint candidat homologous sequences 

import pandas as pd 
#Open the TE results 
TE_table_dir="/beegfs/data/bguinet/these/Repeat_env_analysis2/"
#Gather all .tab files within a directory 
TE_table= pd.DataFrame(columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen','tcov','strand','Newqstart','Newqend','qmatchlen'])
count_table=0
for file in os.listdir(TE_table_dir):
	if file.endswith(".tab"):
		if "TE_" in file:
			continue
		else:
			subTEtab=pd.read_csv(TE_table_dir+file,sep=";")
			TE_table=TE_table.append(subTEtab)
			count_table+=1
print(count_table ," TE tables opened")

TE_table.drop(['qstart','qend'], axis=1, inplace=True)
TE_table.columns=['Scaffold_species','target','pident','alnlen','mismatch','gapopen','tstart','tend','evalue','bits','qlen','tlen','tcov','strand','qstart','qend','qmatchlen']
#TE_table=TE_table.loc[TE_table['tcov'].gt(0.20)]
TE_table=TE_table.loc[TE_table['evalue'].lt(0.0000000001)]



if run_overlapping_TE == "yes"
#######

#Open the table with all homologous filtred results 
Blast_table=Blast_file
Blast_table=pd.read_csv(Blast_table,sep=";")

#Run the overlapping analysis
print("Overlapping analysis running ...")

#Remove duplicate queries 
Blast_table=Blast_table.drop_duplicates(subset=['query'], keep='first')
Blast_table['Scaffold_species']=Blast_table['Scaffold_name']+":"+Blast_table['Species_name']

#Dataframe to fill : 
Overlapping_TEs_within_homologous_viral_sequences_tab = pd.DataFrame(columns=['Clustername','Query','Qstart', 'Qend','Qevalue','cov_depth_BUSCO','cov_depth_scaffold','Overlapping_TE','TE_start','TE_end','TE_evalue','Overlap_percentage','TE_orientation'])

#The following code compare TE positions within scaffolds compared to homologous viral sequence positions.
#The parameter -seqoverlap define the minimum overlapping % of a TE with an homologous viral sequence to count as overlap  

#Remove Homologous with HSPs and without any TEs 
Blast_table=Blast_table.loc[~Blast_table['query'].str.contains("-HSPs")]
Blast_table=Blast_table.loc[Blast_table['Scaffold_species'].isin(TE_table['Scaffold_species'].unique())]

Nb_queries =len(Blast_table['Scaffold_species'].unique())
count=0
for Species in Blast_table['Scaffold_species'].unique():
	print(count,"/",Nb_queries)
	count+=1
	#Extract subtab for each species 
	subBlast_table=Blast_table.loc[Blast_table['Scaffold_species'].str.contains(Species)]
	#Loop over subBlast_table coordinates 
	for index, row in subBlast_table.iterrows():
		subTE=TE_table.loc[TE_table['Scaffold_species'].str.contains(Species)]
		if len(subTE)>0:
			Clustername=row['Clustername']
			cov_depth_BUSCO=row['cov_depth_BUSCO']
			cov_depth_scaffold=row['cov_depth_candidat']
			qevalue=row['evalue']
			Species= row['Scaffold_species']
			Query= row['query']
			start_end=re.sub(".*:","", re.sub("\\(.*","",row['query']))
			Blast_table_start = int(re.sub("-.*","",start_end))
			Blast_table_end = int(re.sub(".*-","",start_end))
			#
			#Full overlapping inside 
			if len(subTE.loc[subTE['qstart'].ge(Blast_table_start) & subTE['qend'].le(Blast_table_end) & subTE['qend'].ge(Blast_table_start) & subTE['qstart'].le(Blast_table_end)]) >0:
				subTE=subTE.loc[subTE['qstart'].ge(Blast_table_start) & subTE['qend'].le(Blast_table_end) & subTE['qend'].ge(Blast_table_start) & subTE['qstart'].le(Blast_table_end)]
				print("---- Full overlapping inside ----")
				for index2, row2 in subTE.iterrows():
					TE_start=row2['qstart']
					TE_end=row2['qend']
					TE_evalue=row2['evalue']
					TE=row2['target']
					Overlap_percentage = ((TE_end - TE_start) / (Blast_table_end - Blast_table_start))*100
					print(Blast_table_start,' - ',Blast_table_end)
					print(TE_start,' - ',TE_end)
					print("Overlapp : " ,Overlap_percentage)
					if Overlap_percentage >= int(seqoverlap):
						Overlapping_TEs_within_homologous_viral_sequences_tab = Overlapping_TEs_within_homologous_viral_sequences_tab.append({'Clustername':Clustername,'Query': Query, 'Qstart':Blast_table_start,'Qend':Blast_table_end,'Qevalue':qevalue , 'cov_depth_BUSCO':cov_depth_BUSCO,'cov_depth_scaffold':cov_depth_scaffold,'Overlapping_TE':TE, 'TE_start':TE_start,'TE_end':TE_end, 'TE_evalue':TE_evalue,'Overlap_percentage':Overlap_percentage,'TE_orientation':'inside'},ignore_index=True)
			#
			#Full overlapping outside 
			if len(subTE.loc[subTE['qstart'].le(Blast_table_start) & subTE['qend'].ge(Blast_table_end) &  subTE['qend'].ge(Blast_table_start) &  subTE['qstart'].le(Blast_table_start)] ) >0:
				subTE=subTE.loc[subTE['qstart'].le(Blast_table_start) & subTE['qend'].ge(Blast_table_end) &  subTE['qend'].ge(Blast_table_start) &  subTE['qstart'].le(Blast_table_start)] 
				print("---- Full overlapping outside ----")
				for index2, row2 in subTE.iterrows():
					TE_start=row2['qstart']
					TE_end=row2['qend']
					TE_evalue=row2['evalue']
					TE=row2['target']
					Overlap_percentage = 100
					print(Blast_table_start,' - ',Blast_table_end)
					print(TE_start,' - ',TE_end)
					print("Overlapp : " ,Overlap_percentage)
					Overlapping_TEs_within_homologous_viral_sequences_tab = Overlapping_TEs_within_homologous_viral_sequences_tab.append({'Clustername':Clustername,'Query': Query, 'Qstart':Blast_table_start,'Qend':Blast_table_end,'Qevalue':qevalue , 'cov_depth_BUSCO':cov_depth_BUSCO,'cov_depth_scaffold':cov_depth_scaffold,'Overlapping_TE':TE, 'TE_start':TE_start,'TE_end':TE_end, 'TE_evalue':TE_evalue,'Overlap_percentage':Overlap_percentage,'TE_orientation':'full'},ignore_index=True)
			#
			#Partial overlapping right
			if len(subTE.loc[subTE['qstart'].ge(Blast_table_start) & subTE['qend'].ge(Blast_table_start) & subTE['qstart'].le(Blast_table_end) & subTE['qend'].ge(Blast_table_end)]) >0:
				subTE=subTE.loc[subTE['qstart'].ge(Blast_table_start) & subTE['qend'].ge(Blast_table_start) & subTE['qstart'].le(Blast_table_end) & subTE['qend'].ge(Blast_table_end)]
				print("---- Partial overlapping right  ----")
				for index2, row2 in subTE.iterrows():
					TE_start=row2['qstart']
					TE_end=row2['qend']
					TE_evalue=row2['evalue']
					TE=row2['target']
					Overlap_percentage =  ((Blast_table_end - TE_start) / (Blast_table_end - Blast_table_start)) *100
					print(Blast_table_start,' - ',Blast_table_end)
					print(TE_start,' - ',TE_end)
					print("Overlapp : " ,Overlap_percentage)
					if Overlap_percentage >= int(seqoverlap):
						Overlapping_TEs_within_homologous_viral_sequences_tab = Overlapping_TEs_within_homologous_viral_sequences_tab.append({'Clustername':Clustername,'Query': Query, 'Qstart':Blast_table_start,'Qend':Blast_table_end,'Qevalue':qevalue , 'cov_depth_BUSCO':cov_depth_BUSCO,'cov_depth_scaffold':cov_depth_scaffold,'Overlapping_TE':TE, 'TE_start':TE_start,'TE_end':TE_end, 'TE_evalue':TE_evalue,'Overlap_percentage':Overlap_percentage,'TE_orientation':'right'},ignore_index=True)
			#
			#Partial overlapping left 
			if len(subTE.loc[subTE['qstart'].le(Blast_table_start) & subTE['qend'].le(Blast_table_end) & subTE['qend'].ge(Blast_table_start) & subTE['qstart'].le(Blast_table_end)]) >0:
				subTE=subTE.loc[subTE['qstart'].le(Blast_table_start) & subTE['qend'].le(Blast_table_end) & subTE['qend'].ge(Blast_table_start) & subTE['qstart'].le(Blast_table_end)]
				print("---- Partial overlapping left  ----")
				for index2, row2 in subTE.iterrows():
					TE_start=row2['qstart']
					TE_end=row2['qend']
					TE_evalue=row2['evalue']
					TE=row2['target']
					Overlap_percentage = ((TE_end - Blast_table_start) / (Blast_table_end - Blast_table_start))*100
					print(Blast_table_start,' - ',Blast_table_end)
					print(TE_start,' - ',TE_end)
					print("Overlap : " ,Overlap_percentage)
					if Overlap_percentage >= int(seqoverlap):
						Overlapping_TEs_within_homologous_viral_sequences_tab = Overlapping_TEs_within_homologous_viral_sequences_tab.append({'Clustername':Clustername,'Query': Query, 'Qstart':Blast_table_start,'Qend':Blast_table_end,'Qevalue':qevalue , 'cov_depth_BUSCO':cov_depth_BUSCO,'cov_depth_scaffold':cov_depth_scaffold,'Overlapping_TE':TE, 'TE_start':TE_start,'TE_end':TE_end, 'TE_evalue':TE_evalue,'Overlap_percentage':Overlap_percentage,'TE_orientation':'left'},ignore_index=True)
			#No overlap 
			if len(Overlapping_TEs_within_homologous_viral_sequences_tab.loc[Overlapping_TEs_within_homologous_viral_sequences_tab['Query']==Query]) ==0:
				TE="None"
				Overlapping_TE="None"
				TE_start="None"
				TE_end="None"
				TE_evalue="None"
				Overlap_percentage="None"
				Overlapping_TEs_within_homologous_viral_sequences_tab = Overlapping_TEs_within_homologous_viral_sequences_tab.append({'Clustername':Clustername,'Query': Query, 'Qstart':Blast_table_start,'Qend':Blast_table_end,'Qevalue':qevalue , 'cov_depth_BUSCO':cov_depth_BUSCO,'cov_depth_scaffold':cov_depth_scaffold,'Overlapping_TE':TE, 'TE_start':TE_start,'TE_end':TE_end, 'TE_evalue':TE_evalue,'Overlap_percentage':Overlap_percentage,'TE_orientation':'None'},ignore_index=True)
		else:
			Clustername=row['Clustername']
			qevalue=row['evalue']
			cov_depth_BUSCO=row['cov_depth_BUSCO']
			cov_depth_scaffold=row['cov_depth_candidat']
			Query=row['query']
			TE="None"
			Overlapping_TE="None"
			TE_start="None"
			TE_end="None"
			TE_evalue="None"
			Overlap_percentage="None"
			Overlapping_TEs_within_homologous_viral_sequences_tab = Overlapping_TEs_within_homologous_viral_sequences_tab.append({'Clustername':Clustername,'Query': Query, 'Qstart':Blast_table_start,'Qend':Blast_table_end,'Qevalue':qevalue , 'cov_depth_BUSCO':cov_depth_BUSCO,'cov_depth_scaffold':cov_depth_scaffold,'Overlapping_TE':TE, 'TE_start':TE_start,'TE_end':TE_end, 'TE_evalue':TE_evalue,'Overlap_percentage':Overlap_percentage,'TE_orientation':'None'},ignore_index=True)

print("\n Overlapping process done  ")

print("Representated clusters : ")

#Save the the Overlapping_TEs_within_homologous_viral_sequences_tab
Overlapping_TEs=Overlapping_TEs_within_homologous_viral_sequences_tab.loc[~Overlapping_TEs_within_homologous_viral_sequences_tab['Overlap_percentage'].str.contains("None", na=False)]

list_to_remove=['Cluster18904', 'blast4058', 'Cluster14435',
'Cluster18974', 'Cluster15629', 'Cluster13425', 'Cluster12627',
'blast5888', 'Cluster1638', 'blast4202', 'Cluster348',
'blast4072', 'Cluster5044', 'blast3899',
'blast6052', 'Cluster6294', 'Cluster7451', 'blast2889',
'Cluster17424', 'Cluster3482', 'Cluster11858', 'Cluster3712',
'Cluster8698', 'Cluster15850', 'blast7673',
'Cluster15568', 'Cluster3966', 'Cluster11808', 'Cluster13439',
'Cluster7069', 'Cluster4987', 'Cluster1861', 'Cluster127',
'Cluster12543', 'Cluster17317', 'blast4467', 'blast6440',
'Cluster9418', 'Cluster3354', 'Cluster9138', 'blast2922',
'Cluster11771', 'Cluster16155', 'blast4397', 'blast1227',
'Cluster683', 'blast811', 'blast4070', 'Cluster3493','Cluster11916','Cluster17599','Cluster12361','Cluster17920']

Filter_Overlapping_TEs = Overlapping_TEs.loc[(Overlapping_TEs['TE_orientation']=="inside") & (Overlapping_TEs['Overlap_percentage'].gt(20)) | (Overlapping_TEs['TE_orientation']=="full") | (Overlapping_TEs['TE_orientation']=="right") & (Overlapping_TEs['Overlap_percentage'].gt(50)) | (Overlapping_TEs['TE_orientation']=="left") & (Overlapping_TEs['Overlap_percentage'].gt(50))]
Filter_Overlapping_TEs.pivot_table(index=['Clustername'], aggfunc='size')

Filter_Overlapping_TEs['Species_name']=Filter_Overlapping_TEs['Query'].str.replace('.*:','')

Filter_Overlapping_TEs=Filter_Overlapping_TEs.loc[~Filter_Overlapping_TEs['Clustername'].isin(list_to_remove)]


Overlapping_TEs_within_homologous_viral_sequences_tab.to_csv(TE_table_dir+"Overlapping_TE_table.tab",sep=";",index=False)
# 665/1087 loci presenting overlapping TEs on Cluster4185 (86/94 species )


if run_count_TE == "yes"
###########
#Count TE in scaffolds containing homologous insertions candidates

#Count number of repeat within scaffolds 

Blast_table="/beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS_synteny_ORFs2.m8"
Blast_table=pd.read_csv(Blast_table,sep=";")
Blast_table=Blast_table.drop_duplicates(subset=['query'], keep='first')
Blast_table['Scaffold_species']=Blast_table['Scaffold_name']+":"+Blast_table['Species_name']
 
TE_table_count=TE_table.groupby(['Scaffold_species']).size().reset_index(name='count')

list_to_remove=['Cluster18904', 'blast4058', 'Cluster14435',
'Cluster18974', 'Cluster15629', 'Cluster13425', 'Cluster12627',
'blast5888', 'Cluster1638', 'blast4202', 'Cluster348',
'blast4072', 'Cluster5044', 'blast3899',
'blast6052', 'Cluster6294', 'Cluster7451', 'blast2889',
'Cluster17424', 'Cluster3482', 'Cluster11858', 'Cluster3712',
'Cluster8698', 'Cluster15850', 'blast7673',
'Cluster15568', 'Cluster3966', 'Cluster11808', 'Cluster13439',
'Cluster7069', 'Cluster4987', 'Cluster1861', 'Cluster127',
'Cluster12543', 'Cluster17317', 'blast4467', 'blast6440',
'Cluster9418', 'Cluster3354', 'Cluster9138', 'blast2922',
'Cluster11771', 'Cluster16155', 'blast4397', 'blast1227',
'Cluster683', 'blast811', 'blast4070', 'Cluster3493','Cluster11916','Cluster17599','Cluster12361','Cluster17920']

TE_merged=TE_table_count.merge(Blast_table,on="Scaffold_species",how="right")
TE_merged['count_repeat1']= TE_merged['count_repeat']
TE_merged['count_repeat']= TE_merged['count']

#Changed scaffold scores 

TE_merged['count_eucaryote_plus_repeat']=TE_merged['count_eucaryote']+TE_merged['count_repeat']
TE_merged.loc[TE_merged['count_eucaryote'].eq(0) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("True"),'Scaffold_score2' ]='X'
TE_merged.loc[TE_merged['count_repeat'].eq(0) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("True"),'Scaffold_score2']='X'

TE_merged.loc[TE_merged['count_eucaryote'].gt(0) & TE_merged['count_eucaryote'].lt(5) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("True"),'Scaffold_score2']='F'
TE_merged.loc[TE_merged['count_repeat'].gt(0) & TE_merged['count_repeat'].lt(2) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("True"),'Scaffold_score2']='F'

TE_merged.loc[TE_merged['count_eucaryote'].eq(0) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("No_cov"),'Scaffold_score2']='E'
TE_merged.loc[TE_merged['count_repeat'].eq(0) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("No_cov"),'Scaffold_score2']='E'

TE_merged.loc[TE_merged['count_eucaryote_plus_repeat'].ge(5) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("True"),'Scaffold_score2']='D'
TE_merged.loc[TE_merged['count_eucaryote'].ge(5) & TE_merged ['FDR_pvalue_cov'].astype(str).str.contains("True"),'Scaffold_score2']='D'
TE_merged.loc[TE_merged['count_repeat'].ge(0) & TE_merged ['FDR_pvalue_cov'].astype(str).str.contains("True"),'Scaffold_score2']='D'

TE_merged.loc[TE_merged['count_eucaryote'].eq(0) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("False"),'Scaffold_score2']='C'
TE_merged.loc[TE_merged['count_repeat'].eq(0) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("False"),'Scaffold_score2']='C'

TE_merged.loc[TE_merged['count_eucaryote'].ge(1) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("No_cov"),'Scaffold_score2']='B'
TE_merged.loc[TE_merged['count_repeat'].ge(1) & TE_merged ['FDR_pvalue_cov'].astype(str).str.contains("No_cov"),'Scaffold_score2']='B'

TE_merged.loc[TE_merged['count_eucaryote'].ge(1) & TE_merged['FDR_pvalue_cov'].astype(str).str.contains("False"),'Scaffold_score2']='A'
TE_merged.loc[TE_merged['count_repeat'].ge(1) & TE_merged ['FDR_pvalue_cov'].astype(str).str.contains("False"),'Scaffold_score2']='A'

TE_merged['Scaffold_score']=TE_merged['Scaffold_score2']

#TE_merged.to_csv("/beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS_synteny_ORFs3.m8",sep=";",index=False)


############Find TE repeat around EVEs,  BUSCOs and random regions  #####

if run_flanking_TE == "yes" 

#This file was generated under R  and contains only A,B,C,D scaffolds 
Blast_table=pd.read_csv("/beegfs/data/bguinet/these/Clustering4/Env_table2.txt",sep=",")
Blast_table=Blast_table.loc[~Blast_table['New_query_bis2'].str.contains("HSP")]
Blast_table['start_end']=Blast_table['New_query_bis2'].str.replace('\\(.*','')
Blast_table['start_end']=Blast_table['start_end'].str.replace('.*:','')
Blast_table['Qstart']=Blast_table['start_end'].str.replace('-.*','')
Blast_table['Qend']=Blast_table['start_end'].str.replace('.*-','')
Blast_table['Qend']=Blast_table['Qend'].astype(int)
Blast_table['Qstart']=Blast_table['Qstart'].astype(int)
Blast_table['Scaffold_name']=Blast_table['New_query_bis2'].str.replace(':.*','')
Blast_table['Scaffold_species']=Blast_table['Scaffold_name']+":"+subBlast_table['Species_name']
Blast_table['query']=Blast_table['New_query_bis2']
m = Blast_table['New_query_bis2'].str.contains('\\(-\\)')
Blast_table['Qstrand'] = np.where(m, '-', '+')
subBlast_table=Blast_table


############Find TE repeat around EVEs #####

TE_distance_table_EVE = pd.DataFrame(columns=['Clustername','Query','Qstart', 'Qend','Qevalue','TE_Distance','TE_start','TE_end','TE_evalue','TE_direction'])

Nb_queries =len(Blast_table['Scaffold_species'].unique())
count=0

TE_distance_table_EVEs = pd.DataFrame(columns=['Species','query','Scaffold_species','Qstart','Qend','TE_start','TE_end','TE_evalue','Distance','Qstrand','TE_name'])

#Subset the TEs within EVEs scaffolds 
subTE=TE_table.loc[TE_table['Scaffold_species'].isin(subBlast_table['Scaffold_species'])]
subTE.rename({'qend':'TE_end','qstart':'TE_start','evalue':'TE_evalue','strand':'TE_strand','target':'TE_name'},axis = 1, inplace = True)
#subBlast_table['Qstart'] = subBlast_table['start']
#subBlast_table['Qend'] = subBlast_table['end']


subBlast_table=subBlast_table[~subBlast_table['query'].str.contains("HSPs")]

#Extract closest TEs for each EVEs
#subBlast_table.rename({'Newend':'Qend','Newstart':'Qstart'},axis = 1, inplace = True)
new_tab = subBlast_table.merge(subTE, on='Scaffold_species')
new_tab['diff'] = ((new_tab['Qstart'] - new_tab['Qend']) - (new_tab['TE_start'] - new_tab['TE_end'])).abs()
out = new_tab[new_tab['diff'] == new_tab.groupby('Scaffold_species')['diff'].transform('min')]
cond=[out['TE_start'] > out['Qend'],out['TE_start'] < out['Qend']]
values=[out['TE_start']- out['Qend'],out['Qstart'] - out['TE_end']]
out['Distance']=np.select(cond,values)

#Save it within the global dataframe 
out=out[['query','Scaffold_species','Qstart','Qend','Qstrand','TE_start','TE_end','TE_evalue','TE_strand','TE_name','Distance','genomic_structure','Clustername','Event','Scaffold_score']]

TE_distance_table_EVEs=TE_distance_table_EVEs.append(out)
TE_distance_table_EVEs=TE_distance_table_EVEs.loc[~TE_distance_table_EVEs['Distance'].lt(0)]

TE_distance_table_EVEs=TE_distance_table_EVEs.sort_values(by='Distance', ascending=True)
#Keep shortest per query
TE_distance_table_EVEs = TE_distance_table_EVEs.drop_duplicates(subset=['query'], keep='first')
TE_distance_table_EVEs['Status']=2
TE_distance_table_EVEs['Species'] = TE_distance_table_EVEs['Scaffold_species'].str.replace('.*:','')

#Add censored EVE data 

Censored_EVEs= subBlast_table.loc[~subBlast_table['query'].isin(TE_distance_table_EVEs['query'])][['query','Scaffold_species','Qstart','Qend','Qstrand','Scaffold_length_x','genomic_structure','Clustername','Event','Scaffold_score']]
Censored_EVEs.columns=['query','Scaffold_species','Qstart','Qend','Qstrand','Scaffold_length','genomic_structure','Clustername','Event','Scaffold_score']
Censored_EVEs['Status']=1
Censored_EVEs['Species']=Censored_EVEs['Scaffold_species'].str.replace('.*:','')
Censored_EVEs['Distance']=Censored_EVEs['Scaffold_length']
TE_distance_table_EVEs=TE_distance_table_EVEs.append(Censored_EVEs)
TE_distance_table_EVEs['Type']="EVEs"

#Keep only one representant sequence within events 
TE_distance_table_EVEs=TE_distance_table_EVEs.groupby(['Clustername', 'Event']).sample()
del TE_distance_table_EVEs['Clustername']
del TE_distance_table_EVEs['Event']

TE_distance_table_EVEs.to_csv("/beegfs/data/bguinet/these/Repeat_env_analysis2/TE_distance_table_EVEs.tab",sep=";",index=False)


############Find TE repeat around random  regions################

#This id done using a snakemake job in order to found all TEs within hymenoptera genomes

#Gather all files into one unique dataframe 

from Bio import SeqIO 
Genome_path="/beegfs/data/bguinet/these/Genomes/"
import os
import random
import numpy as np
tab_names=pd.read_csv("/beegfs/data/bguinet/these/Species_genome_names.txt",header=None)
list_names=tab_names[0].unique()

def Assembly_to_choose(sp_names):
      sp=sp_names
      sp=str(sp)
      if os.path.exists(Genome_path+sp+"/"+sp+"_correctedbis.fa"):
       Genome=sp+"_correctedbis.fa"
      elif os.path.exists(Genome_path+sp+"/"+sp+"_corrected2.fa"):
       Genome=sp+"_corrected2.fa"
      elif os.path.exists(Genome_path+sp+"/"+sp+"_corrected.fa"):
       Genome=sp+"_corrected.fa"
      else:
       Genome=sp+".fa"
      return(str(Genome_path+sp+"/"+Genome))

#For each species, sample 1000 random positions in scaffolds, and look for transposable elements around those regions. 
Genome_informations=pd.DataFrame(columns=["Species_name","Scaffold_name","Scaffold_length"])
TE_distance_table_random = pd.DataFrame(columns=['Query','Qstart', 'Qend','Qstrand','TE_Distance','TE_start','TE_end','TE_evalue','TE_direction','TE_name','TE_strand'])
#list_names=["Lasius_niger"]
for sp in list_names:
	print("Getting genomic informations on : ",sp)
	#Load the species's genome 
	Assembly=Assembly_to_choose(sp)
	record= SeqIO.to_dict(SeqIO.parse(Assembly, "fasta")) 
	df=pd.DataFrame(record.items(), columns=['Scaffold_name', 'Scaffold_length'])
	count_df=df.Scaffold_length.str.len()
	count_df=count_df.to_frame()
	df['Scaffold_length']=count_df['Scaffold_length']
	#keep scaffolds with at least 10000pb
	df_20000=df.loc[df['Scaffold_length'].ge(20000)]
	df_10000=df.loc[df['Scaffold_length'].ge(10000)]
	df_5000=df.loc[df['Scaffold_length'].ge(5000)]
	df_1=df.loc[df['Scaffold_length'].ge(1)]
	random_scaffolds_20000=random.choices(df_20000['Scaffold_name'].unique(), k = 250)
	random_scaffolds_10000=random.choices(df_10000['Scaffold_name'].unique(), k = 250)
	random_scaffolds_5000=random.choices(df_5000['Scaffold_name'].unique(), k = 250)
	random_scaffolds_1=random.choices(df_1['Scaffold_name'].unique(), k = 250)
	random_scaffolds=list(itertools.chain(random_scaffolds_1,random_scaffolds_5000,random_scaffolds_10000,random_scaffolds_20000))
	#df.loc[df['Scaffold_length'].isin(random_scaffolds)]
	#Open Repeat BlastX results
	Repeat_tab=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+sp+"/run_repeat/Assembly_"+sp+"_search_on_RepeatPeps_result.m8",sep="\t",header=None)
	Repeat_tab.columns=["query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","qlen","tlen","qlen","qcov","tcov"]
	Repeat_tab=Repeat_tab.loc[Repeat_tab['evalue'].lt(0.0000000001)]
	Repeat_tab=Repeat_tab.loc[~Repeat_tab['target'].str.contains("Unknown")]
	Repeat_tab=df.merge(Repeat_tab,left_on="Scaffold_name",right_on="query")
	sLength=Repeat_tab.shape[0]
	#Change coordinates
	Repeat_tab['strand']=np.where(Repeat_tab["qstart"]>Repeat_tab["qend"],'-','+')
	#Repeat_tab[['qstart', 'qend']] = np.sort(Repeat_tab[['qstart', 'qend']].values, axis=1)
	m = Repeat_tab['strand'].eq('-')
	Repeat_tab['Newqstart'] = np.where(m, Repeat_tab['Scaffold_length'].sub(Repeat_tab['qstart']), Repeat_tab['qstart'])
	Repeat_tab['Newqend'] = np.where(m, Repeat_tab['Scaffold_length'].sub(Repeat_tab['qend']), Repeat_tab['qend'])
	Repeat_tab.drop(['qstart','qend'], axis=1, inplace=True)
	Repeat_tab.rename(columns={'Newqstart': 'qstart','Newqend': 'qend'},inplace=True, errors='raise')
	#Remove overlapping matches
	Repeat_tab=Repeat_tab.loc[Repeat_tab['evalue'].lt(0.0000000001)]
	Repeat_tab = Repeat_tab.sort_values(['query', 'qend', 'qstart'])
	c1 = Repeat_tab['query'].shift() != Repeat_tab['query']
	c2 = Repeat_tab['qend'].shift() - Repeat_tab['qstart'] < 0
	Repeat_tab['overlap'] = (c1 | c2).cumsum()
	#Finally, we get the row with the maximum sum in each group using groupby.
	Repeat_tab['qmatchlen'] = Repeat_tab['qend'].astype(int) - Repeat_tab['qstart'].astype(int)
	Repeat_tab=Repeat_tab.sort_values(['evalue'], ascending=True).groupby('overlap').first()
	Random_position_table=pd.DataFrame(columns=['Query','Scaffold_name','Qstart','Qend','strand','Scaffold_length'])
	i=1
	for scaffolds in random_scaffolds:
		len_scaffold=df[df['Scaffold_name']==scaffolds]['Scaffold_length'].iloc[0]
		#Get random position within the scaffold
		random_position=random.randint(1, len_scaffold) 
		random_name="random_"+str(i)
		Random_position_table=Random_position_table.append({"Query": random_name,"Scaffold_name":  scaffolds,"Qstart":random_position,"Qend":random_position+1,"Qstrand":"+","Scaffold_length":len_scaffold}, ignore_index=True)
		i+=1
	print("All  random positions loaded\n")
	sub_Repeat_tab=Repeat_tab.loc[Repeat_tab['Scaffold_name'].isin(random_scaffolds)]
	sub_Repeat_tab.rename({'qend':'TE_end','qstart':'TE_start','evalue':'TE_evalue','target':'TE_name','strand':'TE_strand','qlen':'Scaffold_length'},axis = 1, inplace = True)
	sub_Repeat_tab=sub_Repeat_tab[['Scaffold_name','TE_start','TE_end','TE_evalue','TE_name','TE_strand','Scaffold_length']]
	#Extract closest TEs for each EVEs
	new_tab = Random_position_table.merge(sub_Repeat_tab, on='Scaffold_name')
	new_tab['diff'] = ((new_tab['Qstart'] - new_tab['Qend']) - (new_tab['TE_start'] - new_tab['TE_end'])).abs()
	out_random = new_tab[new_tab['diff'] == new_tab.groupby('Scaffold_name')['diff'].transform('min')]
	cond=[out_random['TE_start'] > out_random['Qend'],out_random['TE_start'] < out_random['Qend']]
	values=[out_random['TE_start']- out_random['Qend'],out_random['Qstart'] - out_random['TE_end']]
	out_random['Distance']=np.select(cond,values)
	print("Distance to TEs around random position done")
	#Save it within the global dataframe 
	#out_random.columns=['Query','Scaffold_name','Qstart','Qend','strand','TE_start','TE_end','TE_evalue','diff','Distance','TE_name']
	out_random=out_random.loc[~out_random['Distance'].lt(0)]
	out_random=out_random.sort_values(by='Distance', ascending=True)
	#Keep shortest per query
	out_random = out_random.drop_duplicates(subset=['Query'], keep='first')
	out_random['Status']=2
	out_random['Species'] = sp
	#Add censored random data 
	Censored_random= Random_position_table.loc[~Random_position_table['Scaffold_name'].isin(new_tab['Scaffold_name'])][['Query','Scaffold_name','Qstart','Qend','Scaffold_length']]
	Censored_random['Status']=1
	Censored_random['Species']=sp
	Censored_random['Distance']=0
	del out_random["Scaffold_length_y"]
	out_random.rename({'Scaffold_length_x':'Scaffold_length'},axis = 1, inplace = True)
	TE_distance_table_random=TE_distance_table_random.append(out_random)
	TE_distance_table_random=TE_distance_table_random.append(Censored_random)
TE_distance_table_random['Scaffold_species']=TE_distance_table_random['Scaffold_name']+":"+TE_distance_table_random['Species']
TE_distance_table_random['Genomic_structure']="Random"
TE_distance_table_random.drop(['strand'], axis = 1, inplace = True)
TE_distance_table_random.columns=['query', 'Qstart', 'Qend', 'Qstrand','TE_Distance', 'TE_start', 'TE_end', 'TE_evalue', 'TE_direction', 'TE_name', 'TE_strand','Scaffold_name', 'Scaffoldlength','diff', 'Distance', 'Status', 'Species','Scaffold_species', 'Genomic_structure']
#Generate random Qstrand 
del TE_distance_table_random['Qstrand']
TE_distance_table_random['Qstrand'] = TE_distance_table_random['Qstrand']=np.random.choice(['+','-'],len(TE_distance_table_random))

cond=[TE_distance_table_random['TE_end'] < TE_distance_table_random['Qstart'],TE_distance_table_random['TE_start'] > TE_distance_table_random['Qend'] ,TE_distance_table_random['TE_start'].isna() ]
values=["left","right","NA"]
TE_distance_table_random['Direction']=np.select(cond,values)
values=[-TE_distance_table_random['Distance'],TE_distance_table_random['Distance'],TE_distance_table_random['Distance']]
TE_distance_table_random['Distance']=np.select(cond,values)
TE_distance_table_random.to_csv("/beegfs/data/bguinet/these/Repeat_env_analysis2/TE_distance_table_random.tab",sep=";",index=False)

############Find TE repeat around BUSCOs #####

#Load BUSCO files 
Species_name_file= "/beegfs/data/bguinet/these/Species_genome_names.txt"

list_of_names=[]
for names in open(Species_name_file,"r"):
        list_of_names.append(names.replace("\n", ""))

TE_distance_table_BUSCO = pd.DataFrame(columns=['Species','Scaffold_species','Qstart','Qend','Qstrand','Scaffold_length','TE_start','TE_end','TE_evalue','TE_strand','TE_name','Distance','Status'])


def Assembly_to_choose(sp_names):
      sp=sp_names
      sp=str(sp)
      if os.path.exists(Genome_path+sp+"/"+sp+"_correctedbis.fa"):
       Genome=sp+"_correctedbis.fa"
      elif os.path.exists(Genome_path+sp+"/"+sp+"_corrected2.fa"):
       Genome=sp+"_corrected2.fa"
      elif os.path.exists(Genome_path+sp+"/"+sp+"_corrected.fa"):
       Genome=sp+"_corrected.fa"
      else:
       Genome=sp+".fa"
      return(str(Genome_path+sp+"/"+Genome))

Genome_path="/beegfs/data/bguinet/these/Genomes/"
from Bio import SeqIO 
#list_of_names=['Lasius_niger']
for sp in list_of_names:
	Species=sp
	print("Running TE research on : ",Species)
	#Open BUSCO table 
	if sp in ['Bombus_impatiens', 'Cremastinae_A', 'Ganaspis_brasiliensis']:
		print("1")
		BUSCO_table=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+sp+"/run_busco/run_BUSCO_v3/full_table_"+sp+"_BUSCO_v3_corrected_GC_cov.tsv",sep="\t")
		BUSCO_table.rename({"scaf_length":"Scaffold_length"},axis = 1, inplace = True)
		BUSCO_table=BUSCO_table.loc[BUSCO_table['Status'].str.contains("Complete",na=False)]
		BUSCO_table['Scaffold_species'] = BUSCO_table['scaf_name']+":"+sp
		BUSCO_table['Qstrand'] =  "+" 
		BUSCO_table=BUSCO_table.loc[~BUSCO_table['Scaffold_species'].isna()]
		BUSCO_table['query']=BUSCO_table['Scaffold_species']+":"+BUSCO_table['Busco_id']
		BUSCO_table=BUSCO_table.sample(n=1000)
		Assembly=Assembly_to_choose(sp)
		record= SeqIO.to_dict(SeqIO.parse(Assembly, "fasta")) 
		Scaf_length_table = pd.DataFrame(columns=['Contig','Scaffold_length'])
		for scaf in BUSCO_table['scaf_name'].unique():
			scaf_len=len(str(record[scaf].seq))
			Scaf_length_table = Scaf_length_table.append({'scaf_name':scaf,'Scaffold_length': scaf_len},ignore_index=True)
			#print(Scaf_length_table)
		del BUSCO_table['Scaffold_length']
		BUSCO_table=BUSCO_table.merge(Scaf_length_table,on="scaf_name")
 		print(BUSCO_table)
		BUSCO_table=BUSCO_table.sample(n=1000)
	else:
		try:
			print("2")
			BUSCO_table=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+sp+"/run_busco/run_BUSCO_v3/full_table_"+sp+"_BUSCO_v3_corrected_GC_cov.tsv",sep="\t")
			BUSCO_table.rename({"scaf_length":"Scaffold_length"},axis = 1, inplace = True)
			BUSCO_table=BUSCO_table.loc[BUSCO_table['Status'].str.contains("Complete",na=False)]
			BUSCO_table['Scaffold_species'] = BUSCO_table['scaf_name']+":"+sp
			BUSCO_table=BUSCO_table.loc[~BUSCO_table['Scaffold_species'].isna()]
			BUSCO_table['query']=BUSCO_table['Scaffold_species']+":"+BUSCO_table['Busco_id']
			BUSCO_table['Qstrand'] =  "+" 
			BUSCO_table=BUSCO_table.sample(n=1000)
		except:
				print("3")
				BUSCO_table=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+sp+"/run_busco/run_BUSCO_v3/full_table_"+sp+"_BUSCO_v3_newcoord.tsv",sep=";")
				#BUSCO_table.rename({"scaf_length":"Scaffold_length"},axis = 1, inplace = True)
				BUSCO_table=BUSCO_table.loc[BUSCO_table['Status'].str.contains("Complete",na=False)]
				BUSCO_table['Scaffold_species'] = BUSCO_table['Contig']+":"+sp
				BUSCO_table['query']=BUSCO_table['Scaffold_species']+":"+BUSCO_table['Contig']
				BUSCO_table=BUSCO_table.loc[~BUSCO_table['Scaffold_species'].isna()]
				BUSCO_table=BUSCO_table.sample(n=1000)
				BUSCO_table['Qstrand'] =  "+" 
				Assembly=Assembly_to_choose(sp)
				record= SeqIO.to_dict(SeqIO.parse(Assembly, "fasta")) 
				Scaf_length_table = pd.DataFrame(columns=['Contig','Scaffold_length'])
				for scaf in BUSCO_table['New_query'].unique():
					scaf_len=len(str(record[scaf].seq))
					Scaf_length_table = Scaf_length_table.append({'Contig':scaf,'Scaffold_length': scaf_len},ignore_index=True)
					#print(Scaf_length_table)
				BUSCO_table=BUSCO_table.merge(Scaf_length_table,on="Contig")
	#Subset the TEs within BUSCO scaffolds 
		BUSCO_table['query']=BUSCO_table['Scaffold_species']+":"+BUSCO_table['Busco_id']
	subTE=TE_table.loc[TE_table['Scaffold_species'].str.contains(sp)]
	subTE.rename({"qlen":"Scaffold_length"},axis = 1, inplace = True)
	subTE.rename({'qend':'TE_end','qstart':'TE_start','evalue':'TE_evalue','strand':'TE_strand','target':'TE_name'},axis = 1, inplace = True)
	#Extract closest TEs for each BUSCOs
	BUSCO_table.rename({'Newend':'Qend','Newstart':'Qstart'},axis = 1, inplace = True)
	new_tab = BUSCO_table.merge(subTE, on='Scaffold_species')
	new_tab['diff'] = ((new_tab['Qstart'] - new_tab['Qend']) - (new_tab['TE_start'] - new_tab['TE_end'])).abs()
	out = new_tab[new_tab['diff'] == new_tab.groupby('Scaffold_species')['diff'].transform('min')]
	cond=[out['TE_start'] > out['Qend'],out['TE_start'] < out['Qend']]
	values=[out['TE_start']- out['Qend'],out['Qstart'] - out['TE_end']]
	out['Distance']=np.select(cond,values)
	#Save it within the global dataframe 
	out=out[['query','Scaffold_species','Qstart','Qend','Qstrand','TE_start','TE_end','TE_evalue','Distance','TE_name','TE_strand','Scaffold_length_x']]
	out.rename({"Scaffold_length_x":"Scaffold_length"},axis = 1, inplace = True)
	out=out.sort_values(by='Distance', ascending=True)
	out['Status']=2
	######Add censored BUSCO data
	Censored_BUSCO=BUSCO_table.loc[~BUSCO_table['query'].isin(out['query'])][['query','Scaffold_species','Qstart','Qend','Scaffold_length','Qstrand']]
	Censored_BUSCO.columns=['query','Scaffold_species','Qstart','Qend','Scaffold_length','Qstrand']
	Censored_BUSCO['Status']=1
	out=out.append(Censored_BUSCO)
	out['Species']=sp
	TE_distance_table_BUSCO=TE_distance_table_BUSCO.append(out)

#Remove TEs inside BUSCOs
TE_distance_table_BUSCO=TE_distance_table_BUSCO.loc[~TE_distance_table_BUSCO['Distance'].lt(0)]

#Filter results 
TE_distance_table_BUSCO=TE_distance_table_BUSCO.sort_values(by='Distance', ascending=True)
#Keep shortest per query
TE_distance_table_BUSCO = TE_distance_table_BUSCO.drop_duplicates(subset=['query'], keep='first')
TE_distance_table_BUSCO['Type']="BUSCO"
TE_distance_table_BUSCO.to_csv("/beegfs/data/bguinet/these/Repeat_env_analysis2/TE_distance_table_BUSCO.tab",sep=";",index=False)

#######
###Merge both BUSCO and EVEs flanking TE results :
###

TE_distance_table_BUSCO['genomic_structure'] ="Eukaryote"
TE_distance_table_BUSCO['Scaffold_score']="A"
TE_distance_table_BUSCO_and_EVE = TE_distance_table_BUSCO.append(TE_distance_table_EVEs)

cond=[TE_distance_table_BUSCO_and_EVE['TE_end'] < TE_distance_table_BUSCO_and_EVE['Qstart'],TE_distance_table_BUSCO_and_EVE['TE_start'] > TE_distance_table_BUSCO_and_EVE['Qend'] ,TE_distance_table_BUSCO_and_EVE['TE_start'].isna() ]
values=["left","right","NA"]
TE_distance_table_BUSCO_and_EVE['Direction']=np.select(cond,values)
values=[-TE_distance_table_BUSCO_and_EVE['Distance'],TE_distance_table_BUSCO_and_EVE['Distance'],TE_distance_table_BUSCO_and_EVE['Distance']]
TE_distance_table_BUSCO_and_EVE['Distance']=np.select(cond,values)

TE_distance_table_BUSCO_and_EVE.to_csv("/beegfs/data/bguinet/these/Repeat_env_analysis2/TE_distance_table_BUSCO_and_EVE.tab",sep=";",index=False)
########
