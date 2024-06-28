import pandas as pd 
import argparse
import os
import io 
from Bio import SeqIO

# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Filter loci with NR results                               .\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informations to blast file')
parser.add_argument("-bl", "--blast", help="blast file strand summary") 
parser.add_argument("-Nr", "--Nr_blast", help="The Nr table result file after adding taxid info") 
parser.add_argument("-fa", "--fasta", help="Fasta file of candidate loci") 
parser.add_argument("-o", "--out_file", help="Nr blast file filtered") 
parser.add_argument("-ofa", "--out_fasta", help="Fasta file of candidate loci filtered")
args = parser.parse_args()

#Open mmseqs2 NR results file 
NR_blast=pd.read_csv(arg.Nr_blast,sep=";")

#If only viruses hits keep all clusters 
NR_blast = NR_blast.loc[~(NR_blast['evalue'].gt(0.0000001) & NR_blast['superkingdom'].isin(['Eukaryota', 'Archaea', 'Bacteria']))]

#Sort evalues 
NR_blast=NR_blast.sort_values(by='evalue', ascending=True)
NaN_NR_blast = NR_blast.loc[NR_blast['species'].isna()]
NR_blast = NR_blast.drop_duplicates(subset=['query', 'species'], keep='first')
NR_blast = NR_blast.append(NaN_NR_blast)
NR_blast.superkingdom = NR_blast.superkingdom.fillna('Unknown')

#Create table with Nr match informations 
#Take into account only non-hymenoptera matches 
NR_blast['count_Diptera'] = NR_blast['order'].eq('Diptera').groupby(NR_blast['query']).transform('sum')
NR_blast= NR_blast.loc[~NR_blast['order'].str.contains("Diptera",na=False)]
NR_blast = NR_blast.join(pd.crosstab(NR_blast['query'], NR_blast['superkingdom']).add_prefix('Count_'), on='query')
NR_blast['Count_Bacteria_Eukaryota_Archeae'] = NR_blast['Count_Eukaryota'] + NR_blast['Count_Bacteria']+ NR_blast['Count_Archaea']
NR_blast['Dif_Bacteria_Eukaryota_Archeae_Viruses'] = NR_blast['Count_Bacteria_Eukaryota_Archeae'] - NR_blast['Count_Viruses']

#Remove loci if to much bacteria hits compared to viruses (>20 +) and if less then 30 viruses hits 
NR_blast_remove1=NR_blast.loc[NR_blast['Count_Bacteria_Eukaryota_Archeae'].ge(10)]
NR_blast_remove1=NR_blast_remove1.loc[~NR_blast_remove1['Dif_Bacteria_Eukaryota_Archeae_Viruses'].lt(20)]
NR_blast_remove1=NR_blast_remove1.loc[~NR_blast_remove1['Count_Viruses'].gt(30)]

#If less then 3 virus hit and more then 5 euk hit remove 
NR_blast_remove2=NR_blast.loc[NR_blast['Count_Viruses'].lt(3) & NR_blast['Count_Bacteria_Eukaryota_Archeae'].ge(5)]
NR_blast_remove2=NR_blast_remove2.loc[~NR_blast_remove2['Count_Bacteria_Eukaryota_Archeae'].lt(15)]

#If more than 100 bacteria hits, remove 
NR_blast_remove3=NR_blast.loc[NR_blast['Count_Bacteria'].gt(100)]

#Create a table 'NR_blast_remove' composed of false positves loci 
NR_blast_remove=NR_blast_remove1.append(NR_blast_remove2)
NR_blast_remove=NR_blast_remove.append(NR_blast_remove3)

#Keep unknown loci 
tab_unknown = NR_blast.loc[~NR_blast['query'].isin(NR_blast_remove['query'])]
tab_unknown = tab_unknown.drop_duplicates(subset = 'query', keep = 'first')

list_unknown= list(tab_unknown.loc[tab_unknown['Count_Viruses'].lt(3) & tab_unknown['Count_Unknown'].gt(1)]['query'].unique())

#Keep one representing loci to remove 
NR_blast= NR_blast.drop_duplicates(subset = 'query', keep = 'first')

out_file=args.out_file
NR_blast.to_csv(out_file,sep=";",Index=FALSE)

#Filter fasta file of candidate loci using filtered Blast
Fasta_file=open(args.fasta, "r")
record_dict_aa = SeqIO.to_dict(SeqIO.parse(Fasta_file_aa, "fasta"))
record_dict_dna = SeqIO.to_dict(SeqIO.parse(Fasta_file_dna, "fasta"))

#Save AA format
with open(args.out_fasta_aa,"w") as output:
    for i in Blast_filtered['query'].unique():
        print('>',record_dict_aa['query'].id,sep="",file=output)
        print(record_dict['query'].seq,file=output)

#Save DNA format
with open(args.out_fasta_dna,"w") as output:
    for i in Blast_filtered['query'].unique():
        print('>',record_dict_dna['query'].id,sep="",file=output)
        print(record_dict['query'].seq,file=output)
