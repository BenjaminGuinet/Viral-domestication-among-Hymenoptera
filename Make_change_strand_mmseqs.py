#!/usr/bin/python
import pandas as pd
import numpy as np
import sys 
import argparse
import os

# Print out a message when the program is initiated.
print('----------------------------------------------------------------\n')
print('                        Make change strand of mmseqs2 output.\n')
print('----------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add strand information and change coordinate position in the blast output file')
parser.add_argument("-b", "--blast_file", help="introduce the blast file in .m8 format")
parser.add_argument("-o", "--out",help="The ouptut path where to create the strand.m8 file")
parser.add_argument("-t", "--species_type", help="introduce the name of db species types, ex: Virus or Hymenoptera")
args = parser.parse_args()

"""
python3 Make_change_strand_mmseqs2.py -b /beegfs/data/bguinet/these/Genomes/Leptopilina_boulardi/run_mmseqs2_V/result_mmseqs2.m8 -o /beegfs/data/bguinet/these/Genomes/Leptopilina_boulardi/run_mmseqs2_V -t virus; done
done
"""


# Variable that stores fasta sequences
blast_file=args.blast_file
Output_path=args.out
Type_species=args.species_type
#Allow to get a letter to know where to take and save the files
if Type_species == "virus":
    Type_species = "V"
if Type_species == "Virus":
    Type_species = "V"
if Type_species == "Hymenoptera":
    Type_species = "H"
if Type_species == "hymenoptera":
    Type_species = "H"

if Type_species == "H":

	blast=pd.read_table(blast_file,header=None,comment='#')
	blast.columns = ["query", "target", "pidentity", "alnlength", "mismatches", "gapopens", "qstart", "qend", "sstart", "ssend", "evalue,","bit_score"]
	sLength=blast.shape[0]
	#blast['sstart'] = blast['sstart'].apply(lambda x: x - 1)#to correct the issue 
	#blast['ssend'] = blast['ssend'].apply(lambda x: x + 2)#because the coordinate are from codons
	blast['strand']=np.where(blast["sstart"]>blast["ssend"],'-','+')
	blast[['sstart', 'ssend']] = np.sort(blast[['sstart', 'ssend']].values, axis=1)
	blast_file=blast_file.replace('.m8','')
	blast.to_csv(blast_file+"_strand_H.m8", sep='\t')
	print("New file with strand modifications written at: ", blast_file+"_strand_H.m8")

if Type_species == "V":
       
	blast=pd.read_table(blast_file,header=None)
	blast.columns = ["query", "target", "pident", "alnlen", "mismatch", "gapopen","qstart", "qend","qlen", "tstart", "tend","tlen", "evalue", "bits"]
	sLength=blast.shape[0]
	
	blast['strand']=np.where(blast["qstart"]>blast["qend"],'-','+')
	blast[['qstart', 'qend']] = np.sort(blast[['qstart', 'qend']].values, axis=1)
	m = blast['strand'].eq('-')
	#blast[['qstart','qend']]=np.where(blast['strand'].eq('-')[:,None],np.column_stack((blast['qstart']+1,blast['qend']+0)),blast[['qstart','qend']].values) #In order to mathc the 0-based method of samtools 
        #blast[['qstart','qend']]=np.where(blast['strand'].eq('+')[:,None],np.column_stack((blast['qstart']-1,blast['qend']+0)),blast[['qstart','qend']].values)
	blast['Newqstart'] = np.where(m, blast['qlen'].sub(blast['qend']), blast['qstart'])
	blast['Newqend'] = np.where(m, blast['qlen'].sub(blast['qstart']), blast['qend'])
#	blast['Newqstart'] = blast['sstart'].apply(lambda x: x - 1)
	blast_file=blast_file.replace('.m8','')
	blast.to_csv(blast_file+"_strand_V.m8", sep='\t')
	print("New file with strand modifications written at: ", blast_file+"_strand_V.m8")

if Type_species == "V2":
       
	blast=pd.read_table(blast_file,header=None)
	blast.columns = ["query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits","tlen"]
	sLength=blast.shape[0]
	#blast['qstart'] = blast['qstart'].apply(lambda x: x - 1)#to correct the issue 
	#blast['qend'] = blast['qend'].apply(lambda x: x + 2)#because the coordinate are from codons
	blast['strand']=np.where(blast["qstart"]>blast["qend"],'-','+')
	blast[['qstart', 'qend']] = np.sort(blast[['qstart', 'qend']].values, axis=1)
	m = blast['strand'].eq('-')
	blast['Newqstart'] = np.where(m, blast['qlen'].sub(blast['qend']), blast['qstart'])
	blast['Newqend'] = np.where(m, blast['qlen'].sub(blast['qstart']), blast['qend'])
	blast_file=blast_file.replace('.m8','')
	blast.to_csv(blast_file+"_strand_V.m8", sep='\t')
	print("New file with strand modifications written at: ", blast_file+"_strand_V.m8")

if Type_species == "Repeat":
	blast=pd.read_table(blast_file,header=None)
	blast.columns = ['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen', 'tcov']
	sLength=blast.shape[0]
	blast['strand']=np.where(blast["qstart"]>blast["qend"],'-','+')
	blast[['qstart', 'qend']] = np.sort(blast[['qstart', 'qend']].values, axis=1)
	m = blast['strand'].eq('-')
	blast['Newqstart'] = np.where(m, blast['qlen'].sub(blast['qend']), blast['qstart'])
	blast['Newqend'] = np.where(m, blast['qlen'].sub(blast['qstart']), blast['qend'])
	blast_file=blast_file.replace('.m8','')
	blast.to_csv(blast_file+"_strand_V.m8", sep='\t')
	print("New file with strand modifications written at: ", blast_file+"_strand_V.m8")
	
