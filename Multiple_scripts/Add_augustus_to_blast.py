import argparse
import pandas as pd
import numpy as np
# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Merge Augustus results.\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-s", "--species_name_file", help="introduce the txt file with species names")
parser.add_argument("-b", "--blast_file", help="the blast file")
parser.add_argument("-o", "--out_file", help="The blast output file with augustus informations")
parser.add_argument("-a", "--augustus_file", help="The augustus file (reduced)")
parser.add_argument("-r", "--repeat_file", help="The repeat file (reduced)")
args = parser.parse_args()

Species_name_file=args.species_name_file
Blast_file=args.blast_file
Augustus_file=args.augustus_file
Out_file=args.out_file
Repeat_file=args.repeat_file
#Example usage python3 Add_augustus_to_blast.py -s /beegfs/data/bguinet/these/Species_genome_names_without_outgroup.txt -b /beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER.m8 -o /beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus.m8 -a /beegfs/data/bguinet/these/Genomes/Augustus_search_analysis/Augustus_Uniprot_tax_result_nr2_sub.m8 -r /beegfs/data/bguinet/these/Repeat_env_analysis/All_scaffold_repeat_search_result_dna_strand_V_reduced.m8


#Here we take the 5 first best hits for each Augustus prediction, then if at leats one of them is annoted as Arthropoda, we conclude that this augustus prediction is indeed a good gene prediction. 

Tab=pd.read_csv(Augustus_file,sep="\t",header=None)
Tab.columns = ['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','taxid','taxlineage','taxname']
Tab=Tab.loc[Tab['taxlineage'].str.contains('Arthropoda') | Tab['taxlineage'].str.contains('Insecta') ]
Tab=Tab.drop_duplicates(subset ="query")

Tab['tcov']=(Tab['alnlen']*100)/Tab['tlen']
Tab=Tab.drop(Tab[(Tab.evalue > 0.00001)].index)
Tab=Tab.drop(Tab[(Tab.tcov < 10)].index)
Tab=Tab.drop(Tab[(Tab.tcov < 30) & (Tab.pident < 50)].index)

#Count number of Augustus insect/arthopods gene within each scaffold containing candidate EVEs 
Blast=pd.read_csv(Blast_file,sep=";")

Blast=Blast[['query']]
Blast=pd.concat([Blast[['query']], Blast['query'].str.split(':', expand=True)], axis=1)
Blast = Blast.drop_duplicates(subset = "query", keep="first")

list_of_names1=[]
for names in open(Species_name_file,"r"):
        list_of_names1.append(names)
list_of_names2=[]
for names in list_of_names1:
        list_of_names2.append(names.replace("\n", ""))

      
Table_eucaryote_final = pd.DataFrame() #creates a new dataframe that's empty

Table_augustus_eucaryote_cordinates = pd.DataFrame() #creates a new dataframe that's empty

      
for Species_names in list_of_names2:
  Sub_blast=Blast.loc[Blast[2]==Species_names]
  List_scaff_analyse=[]
  for index, row in Sub_blast.iterrows():
      List_scaff_analyse.append(row[0])
  List_scaff_analyse=list(dict.fromkeys(List_scaff_analyse))
  #Open the augustus table in order to get the gene names 
  Augustus=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+str(Species_names)+"/run_augustus/run_augustus_new.out",comment="#",sep="\t",header=None)
  Augustus['Gene_name'] = Augustus[8].str.extract('(g\d+)')
  Augustus['Complet_name'] = Augustus['Gene_name']+":"+Augustus[0]+":"+Augustus[3].astype(str)+"-"+Augustus[4].astype(str)+"("+Augustus[6]+"):"+Species_names
  Augustus=Augustus.drop_duplicates(subset='Gene_name', keep="first")
  #Augustus.to_csv("/beegfs/data/bguinet/these/Eucaryotic_env_analysis/Augustus_eucaryotic_virus_genes_with_coordinates.tab",sep=";",index=False)
  Augustus=Augustus[[0,'Gene_name','Complet_name']]
  Augustus.columns=['Scaff_name','Gene_name','Complet_name_augustus']
  Augustus['Gene_species_name']=Augustus['Gene_name']+":"+str(Species_names)
  #Merge Augustus and tab in order to add scaffold informations 
  Table2_2 = pd.merge(Tab,Augustus, how='left',left_on=['query'],right_on=['Gene_species_name'])
  Table2_2= Table2_2[Table2_2['Scaff_name'].notna()]
  Table_augustus_eucaryote_cordinates=Table_augustus_eucaryote_cordinates.append(Table2_2)
  Table_augustus_eucaryote_cordinates.drop_duplicates(subset ="Complet_name_augustus",keep = 'first',inplace=True) 
  #print(Species_names)
  #if Species_names=="Venturia_canescens":
  # Table_augustus_eucaryote_virus_cordinates.to_csv("/beegfs/data/bguinet/these/Eucaryotic_env_analysis/Table_augustus_eucaryote_cordinates",sep=";")
  Table3=pd.merge(Table2_2,Augustus, how='outer',left_on=['query'],right_on=['Gene_species_name'])
  Table3=Table3.loc[Table3['query'].str.contains(Species_names,na=False)]
  Table3_eucaryote=Table3.loc[Table3['taxlineage'].str.contains("Arthropoda|Insecta",na=False)]
  Table3_eucaryote.drop_duplicates(subset =["target","pident","alnlen","mismatch","Scaff_name_y"],keep = 'first',inplace=True)  
  #Table3_virus=Table3.loc[Table3['kingdom'].str.contains('Viru',na=False)]
  #Count number of eucaryote loci within each scaffold 
  #Table3_virus=Table3_virus.groupby(['Scaff_name_y']).size().reset_index(name='count')
  #Table3_virus.columns=['Scaff_name','count_virus']
  Table3_eucaryote=Table3_eucaryote.groupby(['Scaff_name_y']).size().reset_index(name='count')
  Table3_eucaryote.columns=['Scaff_name','count_eucaryote'] 
  Table3_eucaryote['Species_name']=Species_names
  Table3_eucaryote=Table3_eucaryote.groupby(["Scaff_name","Species_name"]).agg({"count_eucaryote":"sum"}).reset_index()
  Table_eucaryote_final = Table3_eucaryote.append(Table_eucaryote_final, ignore_index = True) #
  print("Augustus insect gene counting done for : ", Species_names)
  
Table_eucaryote_final.to_csv("/beegfs/data/bguinet/these/Eucaryotic_env_analysis/Table_eucaryote_final.txt",sep=";")

#Open Blast file 
Blast=pd.read_csv(Blast_file,sep=";")
Blast[['query','query1','query2','query3']]=pd.concat([Blast[['query']], Blast['query'].str.split(':', expand=True)], axis=1)
Blast['Scaff_name']=Blast['query1']+":"+Blast['query3']
Blast.drop(['query1','query2','query3'], axis = 1, inplace = True) 
Blast_augustus_tab=pd.merge(Table_eucaryote_final,Blast,left_on=['Scaff_name'],right_on=['Scaff_name'],how='outer')


#Adding Repeat to the blast dataframe 

#Load reapeat file 

Repeat_file=pd.read_csv(Repeat_file,sep=";")
#Count number of repeat in each scaffolds 
Repeat_file=Repeat_file.loc[Repeat_file['evalue'].lt(0.00000005)]
Repeat_file=Repeat_file.loc[Repeat_file['pident'].gt(0.59)]
Repeat_file=Repeat_file.loc[Repeat_file['tcov'].gt(0.20)]
#Repeat_file=Repeat_file.drop(Repeat_file[(Repeat_file['alnlen'].lt(150)) & (Repeat_file['tcov'].lt(0.5))].index)
#Repeat_file=Repeat_file.loc[~Repeat_file['alnlen'].lt(150) & Repeat_file['tcov'].lt(0.5)]
#Repeat_file=Repeat_file.loc[Repeat_file['tcov'].gt(0.20) | Repeat_file['evalue'].lt(0.0000000005)]

Repeat_file=Repeat_file.groupby(['query']).size().reset_index(name='count')

Repeat_file.columns=['Scaff_name','count_repeat']
Blast_augustus_repeat_file=pd.merge(Repeat_file,Blast_augustus_tab,left_on=['Scaff_name'],right_on=['Scaff_name'],how='outer')

Blast_augustus_repeat_file['count_repeat'] = Blast_augustus_repeat_file['count_repeat'].replace(np.nan, 0)
Blast_augustus_repeat_file['count_eucaryote'] = Blast_augustus_repeat_file['count_eucaryote'].replace(np.nan, 0)
Blast_augustus_repeat_file.to_csv(Out_file,sep=";",index=False)

print("All Augustus eucaryote genes have been written to : ",Out_file)

