#!/usr/bin/python
import pandas as pd
import numpy as np
import sys
import argparse
import Bio
from Bio import SeqIO
from statistics import mean
import random
import time
import re
import os


# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Add genomic environment informations into MMseqs2 search tab.\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-b", "--blast", help="The blast file in .m8 format (check the column names and modify the colname 'target' wich corresponds to the Acc_number column in the mmseqs2 software")
parser.add_argument("-o", "--out",help="The ouptut path where to create the oufile")
parser.add_argument("-s", "--species_name_file", help="introduce the txt file with species names")
parser.add_argument("-nB","--Generate_new_busco_cordinate", help="yes if you need to generate new coordinates for BUSCO tsv result (because MEC changed the coordinates of scaffolds)")
args = parser.parse_args()


#Usage : python3 Add_genomic_env.py -i /beegfs/home/bguinet/M2_script/short_file_species_name.txt -b /beegfs/home/bguinet/M2_script/Out_file.m8 -o /beegfs/home/bguinet/M2_script/Out_file_env.m8

Species_name_file=args.species_name_file
blast_file=args.blast
out_file=args.out

if args.Generate_new_busco_cordinate:
	Generate_new_busco_cordinate="yes"
else:
	Generate_new_busco_cordinate="no"

#--------------#
start_time = time.time()
#####
# Part in order to get the correct species names format
####
num_species = sum(1 for line in open(Species_name_file)) # count number of lines in the file in order to now the number of decimals / by the first number given
list_of_names1=[]
for names in open(Species_name_file,"r"):
	list_of_names1.append(names)
list_of_names2=[]

for names in list_of_names1:
	list_of_names2.append(names.replace("\n", ""))

count_row = len(list_of_names2)
filecount = 1
#blast_file="/beegfs/data/bguinet/these/Clustering2/Viralprot_vs_Viral_loci_result_all_match_and_cluster_without_duplicated_taxid_ontology_env_transcripto_mono_dNdS.m8"
#blast_file="/beegfs/data/bguinet/these/Clustering/mmseqs2_search_analysis/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid.m8"
#blast_file="/beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER.m8"
blast_file= pd.read_csv(blast_file,sep=";",header=0)
#del blast_file['Unnamed: 0']

#Create an empty dataframe filled with 0

"""
if os.path.exists("/beegfs/data/bguinet/these/Clustering3/All_informations_names_df_save"):
	All_informations_names_df=pd.read_csv("/beegfs/data/bguinet/these/Clustering3/All_informations_names_df", sep='\t')
	List_already_ran_species=list(All_informations_names_df["Scaffold"].str.split(":", n = 1, expand = True)[1].unique())
	list_of_names2=[item for item in list_of_names2 if item not in List_already_ran_species]
	columns=['Scaffold','cov_depth_BUSCO','cov_depth_candidat','pvalue_cov', 'GC_content_BUSCO', 'GC_content_scaffold','pvalue_gc','Number_viral_loci', 'Number_busco_loci','Scaffold_length']
	All_informations_names_df = pd.DataFrame(columns=columns)
	All_informations_names_df = All_informations_names_df.fillna(0)
else:
	columns=['Scaffold','cov_depth_BUSCO','cov_depth_candidat','pvalue_cov', 'GC_content_BUSCO', 'GC_content_scaffold','pvalue_gc','Number_viral_loci', 'Number_busco_loci','Scaffold_length']
	All_informations_names_df = pd.DataFrame(columns=columns)
	All_informations_names_df = All_informations_names_df.fillna(0)
"""
columns=['Scaffold','cov_depth_BUSCO','cov_depth_candidat','pvalue_cov', 'GC_content_BUSCO', 'GC_content_scaffold','pvalue_gc','Number_viral_loci', 'Number_busco_loci','Scaffold_length']
All_informations_names_df = pd.DataFrame(columns=columns)
All_informations_names_df = All_informations_names_df.fillna(0)


#_______________________________________
#Functions 
def read_file(filename):
 data = {}
 current_key = None
 with open(filename) as f:
	 for line in f:
	 	if line.startswith('>'):
	 		current_key = line.strip()[1:]
	 		data[current_key] = []
	 	elif current_key:
	 		data[current_key].append([int(i) for i in line.split(' ')])
 return data
def prepare_data(data):
 records = []
 for key in data:
		if len(data[key]) > 1:
		 # case where a key have multiple lines (we add a suffix)
		 for i, d in enumerate(data[key]):
			 records.append({
			 'scaffold': '{}_{}'.format(key, i),
			 'start': d[0],
			 'end': d[1]})
		elif len(data[key]) == 1:
	 		# otherwise no suffix needed
	 		records.append({
	 		'scaffold': key,
	 		'start': data[key][0][0],
	 		'end': data[key][0][1]})
 return records
#_______________________________________


def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

#list_of_names2=["Leptopilina_clavipes","Cotesia_vestalis","Platygaster_orseoliae","Megastigmus_dorsalis"]
#list_of_names2=["Wasmannia_auropunctata","Tryphoninae_A","Venturia_canescens"]
#list_of_names2=["Megastigmus_dorsalis","Megastigmus_dorsalis"]
for names in list_of_names2:
		if names=="Platygaster_orseoliae":
			Correct_conta="yes"
		else:
			Correct_conta="no"
		d=[]
		print("--------------------------------------------------")
		print("-  BEGENNING of the process for species : ",names," -")
		print("--------------------------------------------------")
		print("\n")
		print("Charging the candidat Viral and Busco scaffolds")
		if os.path.exists("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected2.fa"):
			fasta=1
			for seq_record in SeqIO.parse(("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected2.fa"), "fasta"):
				d.append({'Scaffold':seq_record.id, 'length':len(seq_record)})
			len_count=pd.DataFrame(d)
		elif os.path.exists("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa"):
			fasta=2
			for seq_record in SeqIO.parse(("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa"), "fasta"):
				d.append({'Scaffold':seq_record.id, 'length':len(seq_record)})
			len_count=pd.DataFrame(d)
		elif os.path.exists("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_bis.fa"):
			fasta=3
			for seq_record in SeqIO.parse(("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_bis.fa"), "fasta"):
				d.append({'Scaffold':seq_record.id, 'length':len(seq_record)})
			len_count=pd.DataFrame(d)
		print("\n")
		#Now we will add the number of virus and busco sequences into the scaffolds
		#Virus informations
		#scaff_count_viral=  pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+names+"/run_mmseqs2_V/result_mmseqs2_summary_V.m8", header=0, sep = " ")
		#print("Opening :",blast_file)
		scaff_count_viral= blast_file
		scaff_count_viral=scaff_count_viral[['Clustername','query','qstart','qend']]
		#Keep only scaffold of the target species
		scaff_count_viral=scaff_count_viral[scaff_count_viral['query'].str.contains(names,na=False)]
		scaff_count_viral['seqnames']=scaff_count_viral['query'].str.replace(r':.*', '')
		scaff_count_viral.seqnames = scaff_count_viral.seqnames.astype(str)
		scaff_count_viral=scaff_count_viral.drop_duplicates(subset='query', keep="first")
		scaff_count_viral_bis=pd.merge(len_count, scaff_count_viral, left_on=['Scaffold'],right_on=['seqnames'],how='outer')
		scaff_count_viral=scaff_count_viral_bis
    		#Now we will add the number of virus and busco sequences into the scaffolds
    		#Virus informations
       		#scaff_count_viral=  pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+names+"/run_mmseqs2_V/result_mmseqs2_summary_V.m8", header=0, sep = " ")
	        #print("Opening :",blast_file)
		#Keep only the scaffolds with candidate loci 
		scaff_count_viral=scaff_count_viral[scaff_count_viral['seqnames'].notna()]
		scaff_count_viral['Number_viral_loci'] = scaff_count_viral.groupby('seqnames')['seqnames'].transform('count')
		scaff_count_viral=scaff_count_viral.drop_duplicates(subset='seqnames', keep="first")
		scaff_count_viral = scaff_count_viral[np.isfinite(scaff_count_viral['Number_viral_loci'])]
		#Augustus informations:
		print("Adding Augustus informations")
		print("Adding BUSCO informations")
		if Generate_new_busco_cordinate=="yes":
			scaff_count_busco=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+names+"/run_busco/run_BUSCO_v3/full_table_"+names+"_BUSCO_v3.tsv", comment="#", header=0,sep="\t")
			scaff_count_busco['Busco_id_contig']=scaff_count_busco['Busco_id']+"_"+scaff_count_busco['Contig']
			#Change the scaffold ids for those who where broken during the process
			if os.path.exists("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/intervals.txt"):
			                interval_file="TRUE"
			                New_coordinates_df = pd.DataFrame(columns=['Busco_id', 'Status', 'scaffold', 'start', 'end', 'Score', 'Length', 'start2', 'end2', 'Name2', 'strand', 'length', 'Newstart', 'Newend'])
			                data = read_file("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/intervals.txt")
			                records = prepare_data(data)
			                df2 = pd.DataFrame(records)
			                if len(df2.loc[df2['scaffold'].str.contains("[0-9]_")])!=0:
						#once we have the df we will add two new columns which corresponds to the new minus strand intervals :
			                        total = df2.groupby(df2.scaffold.str.extract('^([^\.]+)')[0])['end'].transform('max')
			                        df2['start2'] = total - df2['start']
			                        df2['end2'] = total - df2['end']
			                        df2.columns = ["scaffold",  "start_plus",     "end_plus",  "start_minus",  "end_minus"]
			                        df2['end_plus']=df2['end_plus']+1
			                        df2=df2[df2['scaffold'].str.contains("[1-9]_(?!.*_)")]
			                        df2['scaffold2']=df2['scaffold']
			                        df2[['scaffold2','sub_scaffold']] = df2.scaffold2.str.split(r'(?<=[0-9])_', expand=True)
			                        scaff_count_busco3= scaff_count_busco[scaff_count_busco['Contig'].isin(df2['scaffold2'])]
			                        grouped = scaff_count_busco3.groupby('Contig')
			                        for query, group in grouped:
			                                df1=scaff_count_busco3.loc[scaff_count_busco3['Contig']==query]
			                                if len(df1)!=0:
			                                        grouped2=df1.groupby('Busco_id')
			                                        length_query=len(grouped2)
			                                        for query2, group in grouped2:
			                                                df1bis=df1.loc[df1['Busco_id']==query2]
			                                                df1bis=df1bis[:1]
			                                                df1bis.columns=['Busco_id', 'Status', 'scaffold', 'start', 'end', 'Score', 'Length', 'Busco_id_contig']
			                                                Newstart = []
			                                                Newend = []
			                                                Newnames=[]
			                                                subdf2=df2.loc[df2['scaffold'].str.contains(df1bis['scaffold'].values[0])][['scaffold','start_plus','end_plus']]
			                                                if len(subdf2.loc[subdf2['start_plus'].le(int(df1bis['start'].values[0])) & subdf2['end_plus'].ge(int(df1bis['end'].values[0]))])!=0:
			                                                        subdf2=subdf2.loc[subdf2['start_plus'].le(int(df1bis['start'].values[0])) & subdf2['end_plus'].ge(int(df1bis['end'].values[0]))]
			                                                else:
			                                                         subdf2=subdf2.loc[subdf2['end_plus'].le(int(df1bis['end'].values[0]))]
			                                                Newstart.append(subdf2['start_plus'].iloc[0])
			                                                Newend.append(subdf2['end_plus'].iloc[0])
			                                                Newnames.append(subdf2['scaffold'].iloc[0])
			                                                df1bis['start2']=Newstart
			                                                df1bis['end2']=Newend
			                                                df1bis['Name2']=Newnames
									#Now we will get the new coordinates fo the candidate genes
			                                                df1bis["start"] = pd.to_numeric(df1bis["start"])
			                                                df1bis["end"] = pd.to_numeric(df1bis["end"])
			                                                df1bis["start2"] = pd.to_numeric(df1bis["start2"])
			                                                df1bis["end2"] = pd.to_numeric(df1bis["end2"])
									#Add length
			                                                df1bis['strand']="+"
			                                                df1bis['length']=np.where(df1bis['strand'].eq('-'), df1bis['end']-df1bis['start'],df1bis['end']-df1bis['start'])
			                                                Newstart = []
			                                                Newend = []
			                                                NS=df1bis['start'].values[0]-df1bis['start2'].values[0]
			                                                NE=NS+df1bis['length'].values[0]
			                                                Newstart.append(NS)
			                                                Newend.append(NE)
			                                                df1bis['Newstart']=Newstart
			                                                df1bis['Newend']=Newend
			                                                New_coordinates_df=New_coordinates_df.append(df1bis)
			                        New_coordinates_df=New_coordinates_df[['Busco_id_contig', 'Newstart', 'Newend', 'Name2']]
			                        New_coordinates_df.columns=['Busco_id_contig', 'Newstart', 'Newend', 'New_query']
			                        scaff_count_busco2  = pd.merge(scaff_count_busco, New_coordinates_df, on='Busco_id_contig',how="outer")
			                        #Fill NaN value (where scaffold did not break)
			                        scaff_count_busco2.New_query = scaff_count_busco2.New_query.fillna(scaff_count_busco2.Contig)
			                        scaff_count_busco2.Newstart = scaff_count_busco2.Newstart.fillna(scaff_count_busco2.Start)
			                        scaff_count_busco2.Newend =  scaff_count_busco2.Newend.fillna(scaff_count_busco2.End)
			                else:
			                        interval_file="FALSE"
			                        scaff_count_busco['New_query']=scaff_count_busco.Contig
			                        scaff_count_busco['Newstart']=scaff_count_busco['Start']
			                        scaff_count_busco['Newend']=scaff_count_busco['End']
			                        scaff_count_busco2=scaff_count_busco
			                        print(scaff_count_busco2)
			else:
			                interval_file="FALSE"
			                scaff_count_busco['New_query']=scaff_count_busco.Contig
			                scaff_count_busco['Newstart']=scaff_count_busco['Start']
			                scaff_count_busco['Newend']=scaff_count_busco['End']
			                scaff_count_busco2=scaff_count_busco
			print(names)
			scaff_count_busco2.to_csv("/beegfs/data/bguinet/these/Genomes/"+names+"/run_busco/run_BUSCO_v3/full_table_"+names+"_BUSCO_v3_newcoord.tsv",sep=";")
		else:
			scaff_count_busco2=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+names+"/run_busco/run_BUSCO_v3/full_table_"+names+"_BUSCO_v3_newcoord.tsv",sep=";")
		scaff_count_busco=scaff_count_busco2
		scaff_count_busco=scaff_count_busco[['Busco_id_contig','Busco_id','Length','Newstart','Newend','New_query']]
		scaff_count_busco.columns = ['Busco_id_contig','Busco_id','Length','start','end','scaffold']
		scaff_count_busco['Number_busco_loci'] = scaff_count_busco.groupby('scaffold')['scaffold'].transform('count')
		scaff_count_busco=scaff_count_busco.drop_duplicates(subset='scaffold', keep="first")
		scaff_count_busco_viral = pd.merge(scaff_count_viral, scaff_count_busco, left_on=['seqnames'],right_on=['scaffold'],how='outer')
		scaff_count_busco_viral=scaff_count_busco_viral.drop(['query','scaffold'], axis=1)
		scaff_count_busco_viral=scaff_count_busco_viral[scaff_count_busco_viral['seqnames'].notna()]
		#Add Gc and cov informations
		#Open the file containing the Gc & Cov informations for all scaffolds
		if os.path.isfile("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/cov_GC_new_mean.tab") ==True:
						scaff_count_GC_cov=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/cov_GC_new_mean.tab", sep = "\t")
						if 'Unnamed: 0' in scaff_count_GC_cov:
							scaff_count_GC_cov=scaff_count_GC_cov.drop(columns=['Unnamed: 0'])
						scaff_count_GC_cov.columns=["scaf_name","GC","Arithmetic_cov_depth","Scaffold_length","Median_cov_depth"] 
						No_mapping="False"
		elif os.path.isfile("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/cov_GC_new.tab") ==True:
						scaff_count_GC_cov=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/cov_GC_new.tab", sep = "\t")
						if 'Unnamed: 0' in scaff_count_GC_cov:
				   			scaff_count_GC_cov=scaff_count_GC_cov.drop(columns=['Unnamed: 0'])
						scaff_count_GC_cov.columns=["scaf_name","GC","Median_cov_depth","Scaffold_length"] 
						No_mapping="False"
		else:
						No_mapping="True"
						print(names,"does not have a mapping analysis")
						print("\n")
		if No_mapping=="False":
			#Get only scaffold cov for busco:
			print("scaff_count_GC_cov",scaff_count_GC_cov)
			print("scaff_count_busco",scaff_count_busco)
			df_cov_gc_scaff_busco = pd.merge(scaff_count_GC_cov, scaff_count_busco, left_on=['scaf_name'],right_on=['scaffold'],how='inner',)
			#Create a list where will be stored all coverage informations for scaffold that contain at least one BUSCO
			list_scaff_count_GC_cov_busco=[]
			#For each cov_dept, put it into the list
			for i in df_cov_gc_scaff_busco['Median_cov_depth']:
				list_scaff_count_GC_cov_busco.append(i)
			if Correct_conta=="yes":
				list_scaff_count_GC_cov_busco_high=[]
				list_scaff_count_GC_cov_busco_low=[]
				for i in df_cov_gc_scaff_busco['Median_cov_depth']:
					if i >= 21:
						list_scaff_count_GC_cov_busco_high.append(i)
					if i < 21:
						list_scaff_count_GC_cov_busco_low.append(i)
				print("Mean cov High:",mean(list_scaff_count_GC_cov_busco_high))
				print("Mean cov Low:",mean(list_scaff_count_GC_cov_busco_low))
			#print(list_scaff_count_GC_cov_busco, file = open("list_scaff_count_GC_cov_busco"+names, "w"))
			#From this list, we will decide if the scaffold containing a viral sequence as stastiticaly the same cov as an eucaryote scaffold:
			#For that we will parse each candidate scaffold an compare it to the theorique distribution:
			#So we create a dataframe that will contain informations of only candidates viral scaffolds:
			df_cov_gc_scaff_viral = pd.merge(scaff_count_GC_cov, scaff_count_viral, left_on=['scaf_name'],right_on=['seqnames'],how='inner')
			list_cov_scaff_viral=[]
			len_scaff_count_GC_cov_busco=len(list_scaff_count_GC_cov_busco)
			mean_list_scaff_count_GC_cov_busco=mean(list_scaff_count_GC_cov_busco)
			#import statistics
			#mean_list_scaff_count_GC_cov_busco= statistics.median(list_scaff_count_GC_cov_busco)
			print("----------------------------------------")
			print("Coverage Filter processing... ")
			print("----------------------------------------")
			print("\n")
			if Correct_conta=="yes":
				for index, row in df_cov_gc_scaff_viral.iterrows():
					if row['Median_cov_depth'] >= 21:
						mean_list_scaff_count_GC_cov_busco=mean(list_scaff_count_GC_cov_busco_high)
						len_scaff_count_GC_cov_busco_high=len(list_scaff_count_GC_cov_busco_high)
						if row['Median_cov_depth'] > mean_list_scaff_count_GC_cov_busco:
							pvalue=(sum(i > row['Median_cov_depth'] for i in list_scaff_count_GC_cov_busco_high)/len_scaff_count_GC_cov_busco_high)
						else:
							pvalue=(sum(i < row['Median_cov_depth'] for i in list_scaff_count_GC_cov_busco_high)/len_scaff_count_GC_cov_busco_high)
						list_cov_scaff_viral.append({'scaf_name':row['scaf_name'],'Scaffold_length':row['Scaffold_length'],'cov_depth_BUSCO':mean_list_scaff_count_GC_cov_busco,'pvalue_cov':pvalue,'cov_depth_candidat':row['Median_cov_depth']})
					if row['Median_cov_depth'] < 21:
						mean_list_scaff_count_GC_cov_busco=mean(list_scaff_count_GC_cov_busco_low)
						len_scaff_count_GC_cov_busco_low=len(list_scaff_count_GC_cov_busco_low)
						if row['Median_cov_depth'] > mean_list_scaff_count_GC_cov_busco:
							pvalue=(sum(i > row['Median_cov_depth'] for i in list_scaff_count_GC_cov_busco_low)/len_scaff_count_GC_cov_busco_low)
						else:
							pvalue=(sum(i < row['Median_cov_depth'] for i in list_scaff_count_GC_cov_busco_low)/len_scaff_count_GC_cov_busco_low)
						list_cov_scaff_viral.append({'scaf_name':row['scaf_name'],'Scaffold_length':row['Scaffold_length'],'cov_depth_BUSCO':mean_list_scaff_count_GC_cov_busco,'pvalue_cov':pvalue,'cov_depth_candidat':row['Median_cov_depth']})
			else:
				for index, row in df_cov_gc_scaff_viral.iterrows():
					if row['Median_cov_depth'] > mean_list_scaff_count_GC_cov_busco:
						pvalue=(sum(i > row['Median_cov_depth'] for i in list_scaff_count_GC_cov_busco)/len_scaff_count_GC_cov_busco)
					else:
						pvalue=(sum(i < row['Median_cov_depth'] for i in list_scaff_count_GC_cov_busco)/len_scaff_count_GC_cov_busco)
					list_cov_scaff_viral.append({'scaf_name':row['scaf_name'],'Scaffold_length':row['Scaffold_length'],'cov_depth_BUSCO':mean_list_scaff_count_GC_cov_busco,'pvalue_cov':pvalue,'cov_depth_candidat':row['Median_cov_depth']})
			print("Coverage informations processing done")
			print("Number of coverage scaffold for BUSCO distribution", len_scaff_count_GC_cov_busco)
			print("Number of candidate scaffolds analyzed: ", len(df_cov_gc_scaff_viral))
			print("\n")
			candidates_cov_scaff=pd.DataFrame(list_cov_scaff_viral)
		else: 
			print("no mapping detected") 
		print("----------------------------------------")
		print(" G+C Filter processsing... ")
		print("----------------------------------------")
		print("\n")
		#first we will open the chimera made by the concatenation of all the scaffols containing at leas one busco sequence.
		chimera= SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Genomes/"+names+"/run_busco/Chimera_Busco_scaff.fa", "fasta"))
		for record in chimera:
		    range_chimera= len(chimera['Chimera_seq'].seq)
		#In order to count the GC content of all scaffolds with a viral seq we have to download the genome of the species
		#Change all in upper cases
		if  fasta==1:
		    records = (rec.upper() for rec in SeqIO.parse("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected2.fa", "fasta"))
		elif  fasta==2:
		    records = (rec.upper() for rec in SeqIO.parse("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa", "fasta"))
		elif  fasta==3:
		    records = (rec.upper() for rec in SeqIO.parse("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_bis.fa", "fasta"))
		else:
		    records = (rec.upper() for rec in SeqIO.parse("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+".fa", "fasta"))
		#Convert to a dict object
		Genome_fasta_species= SeqIO.to_dict(records)
		#for record in Genome_fasta_species:
		count_row = scaff_count_viral.shape[0]
		filecount = 1
		#Create a liste that will be transformed into a tab in order to ad pvalue informations
		list_gc_scaff_viral=[]
		for index, row in scaff_count_viral.iterrows():
		        print(row['Scaffold'])
		    	#For each scaff name, check the content of G+C
		    	#length_sequence_viral_scaff=len(Genome_fasta_species[row['scaf_name']].seq)
		        length_sequence_viral_scaff=row['length']
		        GC_content_scaffold = (Genome_fasta_species[row['Scaffold']] .seq.count('G') + Genome_fasta_species[row['Scaffold']] .seq.count('C'))
		        GC_content_candidat = (Genome_fasta_species[row['Scaffold']] .seq[int(row['qstart']):int(row['qend'])].count('G')+Genome_fasta_species[row['Scaffold']].seq[int(row['qstart']):int(row['qend'])].count('C'))
		        GC_content_scaffold = GC_content_scaffold - GC_content_candidat #In order to not take into account the GC content of the viral tract which can change wether it has been transfered recently or not.
		        length_sequence_viral_scaff= length_sequence_viral_scaff - len(Genome_fasta_species[row['Scaffold']] .seq[int(row['qstart']):int(row['qend'])])
				# The tract to test changes and ecomes the size of the whole scaffold - the size of the viral tract.
		        if length_sequence_viral_scaff ==0 :
		                        length_sequence_viral_scaff=1
		        else:
		                        GC_content_scaffold  = int(GC_content_scaffold)  / int(length_sequence_viral_scaff)
				#Now we will sample several part in the chimera with the same scaffold length as the candidate and get a distribution of G+C values
		    	#Number of sampling
		        nb_interval=500
		    	#Limits from 0 to the length of the chimera sequence
		        limit_low= 0
		        limit_high=range_chimera
		        liste_GCcontent_sample=[]
		        for x in range(nb_interval):
		                        number = random.randint(limit_low,limit_high-length_sequence_viral_scaff)
		    					#Extract the sequence from the coordinates randomly generated across the chimera sequence
		                        random_sequence=chimera['Chimera_seq'].seq[number:number + length_sequence_viral_scaff]
		    					#count GC content of each sample took in chimera
		                        random_sequence=random_sequence.upper()
		                        GC_content_sample=(random_sequence.count('G') + random_sequence.count('C'))/length_sequence_viral_scaff
		                        liste_GCcontent_sample.append(GC_content_sample)
		    	#CAlculation of the probability for an observed value to higher or lower that expected
		    	#We calcul the % of value above or below the observed value. If the % is very low it means that few observed values show the same GC pattern
		    	#we the hypothesis that they are from eucaryotic origine.
		        if GC_content_scaffold > mean(liste_GCcontent_sample):
		                        pvalue_GC=(sum(i > GC_content_scaffold for i in liste_GCcontent_sample)/len(liste_GCcontent_sample))
		        else:
		                        pvalue_GC=(sum(i < GC_content_scaffold for i in liste_GCcontent_sample)/len(liste_GCcontent_sample))
        		list_gc_scaff_viral.append({'scaf_name':row['Scaffold'],'Scaffold_length':row['length'],'pvalue_gc':pvalue_GC,'GC_content_scaffold':GC_content_scaffold,'GC_content_BUSCO':mean(liste_GCcontent_sample)})
		candidates_gc_scaff=pd.DataFrame(list_gc_scaff_viral)
		#print(candidates_gc_scaff)
		print("G+C Filter processing done")
		print("Informations:")
		print("Number of iteration done: ", nb_interval)
		print("Number of candidate scaffolds analyzed: ", count_row)
		print("\n")
		#Merge the two dataframe dogether GC + cov information
		if No_mapping=="True":
			candidates_gc_cov_scaff= candidates_gc_scaff
		else:
			candidates_gc_cov_scaff = pd.merge(candidates_gc_scaff, candidates_cov_scaff, left_on=['scaf_name'],right_on=['scaf_name'],how='outer')
			#print(candidates_gc_cov_scaff)
			#candidates_gc_cov_scaff.drop(['cov_depth_y','Scaffold_length'],axis=1,inplace=True)
			###Merge quantity of viral and busco locus in scaffolds and pvalues of cov and gc
		Count_scaff_gc_cov=pd.merge(candidates_gc_cov_scaff,scaff_count_busco_viral,left_on=['scaf_name'],right_on=['Scaffold'],how='outer',)
		Count_scaff_gc_cov=Count_scaff_gc_cov.fillna("NA")
		    #Change the scaffold names in order to add all these informations into a huge dataframe, then this dataframe at this end will be merged with the blast dataframe
		Count_scaff_gc_cov['Scaffold']=Count_scaff_gc_cov['Scaffold'].astype(str) + ":" + names
		print("Count_scaff_gc_cov\n",Count_scaff_gc_cov)
		print(list(Count_scaff_gc_cov))
		#print("okok",Count_scaff_gc_cov)
		if No_mapping=="True":
			print("No mapping information detected")
		else:
			Count_scaff_gc_cov.drop(["Scaffold_length_x", "scaf_name", "Scaffold_length_y", "scaf_name"],axis=1,inplace=True)
			#Reorder the columns
			print(Count_scaff_gc_cov)
			Count_scaff_gc_cov= Count_scaff_gc_cov[['Scaffold','length','cov_depth_BUSCO','cov_depth_candidat','pvalue_cov', 'GC_content_BUSCO', 'GC_content_scaffold','pvalue_gc','Number_viral_loci', 'Number_busco_loci']]
		All_informations_names_df=pd.concat([All_informations_names_df,Count_scaff_gc_cov], axis=0, join='outer', ignore_index=False,keys=None, levels=None, names=None, verify_integrity=False,copy=True)
		print("G+C and Coverage informations processing done for species :", names )
		print("\n")
		All_informations_names_df['Number_viral_loci'] = All_informations_names_df['Number_viral_loci'].fillna(0)
		All_informations_names_df['Number_busco_loci'] = All_informations_names_df['Number_busco_loci'].fillna(0)
		All_informations_names_df.to_csv("/beegfs/data/bguinet/these/Clustering3/All_informations_names_df_save", sep='\t')

All_informations_names_df.to_csv("/beegfs/data/bguinet/these/Clustering3/All_informations_names_df_huge", sep='\t')

###Merge quantity, gc and cov with blast informations.
print("Merging of G+C and Coverage informations with the MMseqs analysis table")
print("Merging of G+C and Coverage informations with the MMseqs analysis table")
#blast_file= pd.read_csv("/beegfs/data/bguinet/these/Clustering2/Viralprot_vs_Viral_loci_result_all_match_and_cluster_without_duplicated_taxid_ontology.m8",sep=";",header=0)
#blast_file.drop(['Unnamed: 0','Unnamed:_0','Unnamed:_0.1'],axis=1,inplace=True)
###Merge quantity, gc and cov with blast informations.
print("Merging of G+C and Coverage informations with the MMseqs analysis table")
#To do so we will have to transforme the names
#All_informations_names_df=  pd.read_csv("/beegfs/data/bguinet/these/Clustering2/All_informations_names_df", header=0, sep = "\t")
#All_informations_names_df=All_informations_names_df[['GC_content_BUSCO', 'GC_content_scaffold','pvalue_gc','cov_depth_BUSCO', 'cov_depth_candidat','pvalue_cov','scaffnames','scaf_name']]
#All_informations_names_df.drop(["Unnamed: 0","Scaffold_length","qend", "qstart","Clustername"],axis=1,inplace=True)
#Load the table with eucaryotic genes within scaffolds
#
#Eucaryotic_scaff_table=pd.read_csv("/beegfs/data/bguinet/these/Eucaryotic_env_analysis/Table_eucaryote_virus_final.txt",sep=";")
#Eucaryotic_scaff_table['Scaff_and_species']=Eucaryotic_scaff_table['Scaff_name']+":"+Eucaryotic_scaff_table['Species_name']

#All_informations_names_df = pd.merge(Eucaryotic_scaff_table,All_informations_names_df, left_on=["Scaff_and_species"],right_on=['scaffnames'],how='outer')
#blast_file["query_bis"]=blast_file["query"].replace(':[^:]+:',":",regex=True) #Allow to get a proper forme of the query names
#blast_file=blast_file.drop(['count_eucaryote','count_virus','cov_depth_BUSCO','cov_depth_candidat','pvalue_cov','pvalue_gc','GC_content_BUSCO','GC_content_scaffold'], axis=1)
#blast_with_GC_cov_inf = pd.merge(blast_file, All_informations_names_df, left_on=["query_bis"],right_on=['scaffnames'],how='outer')
#blast_with_GC_cov_inf = blast_with_GC_cov_inf[blast_with_GC_cov_inf['query'].notna()]
#blast_with_GC_cov_inf= blast_with_GC_cov_inf.sort_values(by='Clustername', ascending=True)
#blast_with_GC_cov_inf.drop(['scaffnames','query_y','query_bis'],axis=1,inplace=True)
#blast_with_GC_cov_inf.rename(columns={'length': 'Scaf_length'}, inplace=True)
#blast_with_GC_cov_inf.dropna(subset = ['Clustername'],inplace=True)
#blast_with_GC_cov_inf.to_csv("/beegfs/data/bguinet/these/Clustering2/Viralprot_vs_Viral_loci_result_all_match_and_cluster_without_duplicated_taxid_ontology_env.m8", sep=';')
"""

All_informations_names_df= pd.read_csv("/beegfs/data/bguinet/these/Clustering3/All_informations_names_df",sep="\t")

blast_file=pd.read_csv("/beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat.m8",sep=";")
import numpy as np
blast_file["query_bis"]=blast_file["query"].replace(':[^:]+:',":",regex=True)
blast_with_GC_cov_inf = pd.merge(blast_file, All_informations_names_df, left_on=["query_bis"],right_on=['Scaffold'],how='outer')
blast_with_GC_cov_inf = blast_with_GC_cov_inf[blast_with_GC_cov_inf['query'].notna()]
blast_with_GC_cov_inf['Number_busco_loci'] = blast_with_GC_cov_inf['Number_busco_loci'].replace(np.nan, 0)
blast_with_GC_cov_inf['Number_viral_loci'] = blast_with_GC_cov_inf['Number_viral_loci'].replace(np.nan, 0)
del blast_with_GC_cov_inf['Scaffold_length']

blast_with_GC_cov_inf=blast_with_GC_cov_inf.rename(columns={"length": "Scaffold_length"})

#Proceed FDR analysis 


from multipy.fdr import tst

#Check wich FDR to take 
#list_nb=[0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01]

#for i in list_nb:
#	blast_with_GC_cov_inf['FDR_pvalue_cov']=tst(blast_with_GC_cov_inf['pvalue_cov_median'], q=i)
#	print(i,len(blast_with_GC_cov_inf.loc[blast_with_GC_cov_inf['FDR_pvalue_cov'].astype(str)=="False"]))
#Here we take 0.05


#-Pvalue_cov
FDR_with_cov_inf= blast_with_GC_cov_inf.drop_duplicates(subset=['Scaffold'], keep='first')
#FDR_with_cov_inf=FDR_with_cov_inf[~FDR_with_cov_inf['pvalue_cov'].str.contains("No_cov")]
FDR_with_cov_inf['cov_depth_BUSCO'] = FDR_with_cov_inf['cov_depth_BUSCO'].str.replace('No_cov','-1')
FDR_with_cov_inf['pvalue_cov'] = FDR_with_cov_inf['pvalue_cov'].str.replace('No_cov','-1')
FDR_with_cov_inf['pvalue_cov']=FDR_with_cov_inf['pvalue_cov'].astype(float)
FDR_with_cov_inf['cov_depth_BUSCO']=FDR_with_cov_inf['cov_depth_BUSCO'].astype(float)
FDR_with_cov_inf=FDR_with_cov_inf.loc[FDR_with_cov_inf['cov_depth_BUSCO']>0]
FDR_with_cov_inf=FDR_with_cov_inf[['Scaffold','pvalue_cov']]
FDR_with_cov_inf['FDR_pvalue_cov']=tst(FDR_with_cov_inf['pvalue_cov'], q=0.05)
FDR_with_cov_inf=FDR_with_cov_inf[['Scaffold','FDR_pvalue_cov']]
#for i in list_nb:
#	Full_tab_with_dNdS_pvalue_GC['FDR_pvalue_gc']=tst(Full_tab_with_dNdS_pvalue_GC['pvalue_gc'], q=i)
#	print(i,len(Full_tab_with_dNdS_pvalue_GC.loc[Full_tab_with_dNdS_pvalue_GC['FDR_pvalue_gc'].astype(str)=="False"]))
#Here we take 0.05 because no coude 

FDR_with_GC_inf= blast_with_GC_cov_inf.drop_duplicates(subset=['Scaffold'], keep='first')
FDR_with_GC_inf=FDR_with_GC_inf[FDR_with_GC_inf['pvalue_gc'].notna()]
FDR_with_GC_inf=FDR_with_GC_inf[['Scaffold','pvalue_gc']]
FDR_with_GC_inf['FDR_pvalue_gc']=tst(FDR_with_GC_inf['pvalue_gc'], q=0.05)
FDR_with_GC_inf=FDR_with_GC_inf[['Scaffold','FDR_pvalue_gc']]

#Merge cov and gc FDR
blast_with_GC_cov_inf=blast_with_GC_cov_inf.merge(FDR_with_GC_inf,left_on='Scaffold', right_on='Scaffold',how='outer')
blast_with_GC_cov_inf=blast_with_GC_cov_inf.merge(FDR_with_cov_inf,left_on='Scaffold', right_on='Scaffold',how='outer')


#Replace NA by no_cov
blast_with_GC_cov_inf.fillna({'cov_depth_BUSCO': "No_cov", 'cov_depth_candidat': "No_cov",'pvalue_cov':"No_cov", 'FDR_pvalue_cov':"No_cov"}, inplace=True)

blast_with_GC_cov_inf['cov_depth_candidat'] = blast_with_GC_cov_inf['cov_depth_candidat'].str.replace('No_cov','-1')
blast_with_GC_cov_inf['cov_depth_BUSCO'] = blast_with_GC_cov_inf['cov_depth_BUSCO'].str.replace('No_cov','-1')
blast_with_GC_cov_inf['diff_cov']= blast_with_GC_cov_inf['cov_depth_candidat'].astype(float)-blast_with_GC_cov_inf['cov_depth_BUSCO'].astype(float)
blast_with_GC_cov_inf['diff_cov']=abs(blast_with_GC_cov_inf['diff_cov'])
blast_with_GC_cov_inf.loc[blast_with_GC_cov_inf.diff_cov.lt(10) & blast_with_GC_cov_inf.diff_cov.gt(0) ,'FDR_pvalue_cov']='False'
blast_with_GC_cov_inf['cov_depth_candidat'] = blast_with_GC_cov_inf['cov_depth_candidat'].str.replace('-1','No_cov')
blast_with_GC_cov_inf['cov_depth_BUSCO'] = blast_with_GC_cov_inf['cov_depth_BUSCO'].str.replace('-1','No_cov')

blast_with_GC_cov_inf.to_csv("/beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env.m8",index=False,sep=";")
"""
print("\n")
print("Process done")
print("Output written at: ",out_file)
print("--- %s seconds ---" % (time.time() - start_time))



