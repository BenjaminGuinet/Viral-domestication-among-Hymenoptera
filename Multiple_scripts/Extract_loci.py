import re
import textwrap
import os 
import pandas as pd
import pybedtools
from pybedtools import BedTool
from Bio import SeqIO 
import numpy as np 
import sys
sys.tracebacklimit = None
import argparse

# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Add genomic environment informations into MMseqs2 search tab.\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-l", "--loci_file", help="The loci summary file")
parser.add_argument("-s", "--species_name", help="Only used when -c not equal to no")
parser.add_argument("-o","--outfile", help="loci output file")
parser.add_argument("-b","--bedfile_output", help="The bedfile output")
parser.add_argument("-fasta","--fasta_genome", help="The fasta genome file")
parser.add_argument("-c","--correct_scaf_name", help="Put no if you do not need to change coordinatesand scaf name because of MEC")
args = parser.parse_args()

outfile=args.outfile
bedfile_output=args.bedfile_output

try:
	fasta_genome=args.fasta
	correct_scaf_name=args.correct_scaf_name
except:
	print("ok")
	correct_scaf_name="yes"
Loci_file=args.loci_file
names=args.species_name

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
    	 		'end': d[1],
    	 		'frag':'{}'.format(i)})
    		elif len(data[key]) == 1:
    	 		# otherwise no suffix needed
    	 		records.append({
    	 		'scaffold': key,
    	 		'start': data[key][0][0],
    	 		'end': data[key][0][1],
    	 		'frag':"no"})
     return records
#Bedfile=pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+names+"/Recover_fasta_loci.bed",header=None,sep="\t")


if correct_scaf_name =="no":
    Bedfile=pd.read_csv(Loci_file,sep=" ")
    Bedfile.columns=['scaffold','start','end','width','Info','strand']
    print("Generating the bed file...")
    Bedfile['Species']="Species"
    Bedfile['Info']="0"
    Bedfile=Bedfile[['scaffold', 'start','end','Species','Info','strand']]
    Bedfile[['start','end']]=np.where(Bedfile['strand'].eq('-')[:,None],np.column_stack((Bedfile['start']-1,Bedfile['end']+0)),Bedfile[['start','end']].values)
    Bedfile[['start','end']]=np.where(Bedfile['strand'].eq('+')[:,None],np.column_stack((Bedfile['start']-1,Bedfile['end']+0)),Bedfile[['start','end']].values)
    Bedfile.to_csv(bedfile_output,header=None,sep="\t",index=False)
    print("Bedfile written to : ",bedfile_output )
    print("Extracting the fasta sequences...")
    print("Nothing has been corrected")
    a = pybedtools.example_bedtool(bedfile_output)
    fasta = pybedtools.example_filename(fasta_genome)
    try:
          b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
    except Exception:
          pass
    #b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
    if os.path.isfile(outfile):
        print("File create succesfully!")
else:
    Bedfile=pd.read_csv(Loci_file,sep=" ")
    Bedfile.columns=['scaffold','start','end','width','Info','strand']
    print("Analysis the species :",names)
    #Functions 
    Interval2='ok'
    #Change the scaffold ids for those who where broken during the process
    if os.path.exists("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/intervals.txt"):
                    Interval='TRUE'
                    New_coordinates_df = pd.DataFrame(columns=['Busco_id', 'Status', 'scaffold', 'start', 'end', 'Score', 'Length', 'start2', 'end2', 'Name2', 'strand', 'length', 'Newstart', 'Newend'])
                    data = read_file("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/intervals.txt")
                    records = prepare_data(data)
                    df2 = pd.DataFrame(records)
                    #once we have the df we will add two new columns which corresponds to the new minus strand intervals :
                    total = df2.groupby(df2.scaffold.str.extract('^([^\.]+)')[0])['end'].transform('max')
                    df2['start2'] = total - df2['start']
                    df2['end2'] = total - df2['end']
                    df2.columns = ["scaffold",  "start_plus","end_plus",'frag',  "start_minus",  "end_minus"]
                    df2['end_plus']=df2['end_plus']+1
                    df2=df2[~df2['frag'].str.contains("no")]
                    if df2.empty:
                        Interval2='FALSE'
                        Interval='FALSE'
                        Bedfile=Bedfile.drop_duplicates(['scaffold','start','end','strand'], keep="first") 
                        Bedfile.to_csv(bedfile_output,header=None,sep="\t",index=False)
                    else:
                        del df2['frag']
                        df2['scaffold2']=df2['scaffold']
                        df2[['scaffold2','sub_scaffold']] = df2.scaffold2.str.split(r'(?<=[0-9])_', expand=True)
                        scaff_count_busco3= Bedfile[Bedfile['scaffold'].isin(df2['scaffold2'])]
                        grouped = scaff_count_busco3.groupby('scaffold')
                        for query, group in grouped:
                                        df1=scaff_count_busco3.loc[scaff_count_busco3['scaffold']==query]
                                        if len(df1)!=0:
                                                grouped2=df1.groupby('start')
                                                length_query=len(grouped2)
                                                for query2, group in grouped2:
                                                        df1bis=df1.loc[df1['start']==query2]
                                                        df1bis=df1bis[:1]
                                                        #df1bis.columns=['Busco_id', 'Status', 'scaffold', 'start', 'end', 'Score', 'Length', 'Busco_id_contig']
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
                                                        #df1bis['strand']="+"
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
                        New_coordinates_df=New_coordinates_df[['scaffold','Name2', 'Newstart', 'Newend','width','Info','strand']] 
                        #Bedfile[['start','end']]=np.where(Bedfile['strand'].eq('-')[:,None],np.column_stack((Bedfile['start']-1,Bedfile['end']+0)),Bedfile[['start','end']].values)
                        #Bedfile[['start','end']]=np.where(Bedfile['strand'].eq('+')[:,None],np.column_stack((Bedfile['start']-1,Bedfile['end']+0)),Bedfile[['start','end']].values)
                        New_coordinates_df[['Newstart','Newend']]=np.where(New_coordinates_df['Newstart'].eq(0)[:,None],np.column_stack((New_coordinates_df['Newstart']+3,New_coordinates_df['Newend']+0)),New_coordinates_df[['Newstart','Newend']].values)
                        New_coordinates_df2  = pd.merge(Bedfile, New_coordinates_df, left_on=['scaffold'],right_on=['scaffold'],how="outer")
                                        #Fill NaN value (where scaffold did not break)
                        New_coordinates_df2.Name2 = New_coordinates_df2.Name2.fillna(New_coordinates_df2.scaffold)
                        New_coordinates_df2.Newstart = New_coordinates_df2.Newstart.fillna(New_coordinates_df2.start)
                        New_coordinates_df2.Newend =  New_coordinates_df2.Newend.fillna(New_coordinates_df2.end)
                        New_coordinates_df2.strand_y =  New_coordinates_df2.strand_y.fillna(New_coordinates_df2.strand_x)
                        New_coordinates_df2=New_coordinates_df2[['Name2', 'Newstart', 'Newend','strand_y']] 
                        New_coordinates_df2['Newstart']=New_coordinates_df2['Newstart'].astype(int)
                        New_coordinates_df2['Newend']=New_coordinates_df2['Newend'].astype(int)
                        New_coordinates_df2['Species']="Species"
                        New_coordinates_df2['Info']=0
                        New_coordinates_df2=New_coordinates_df2[['Name2', 'Newstart','Newend','Species','Info','strand_y']]
                        New_coordinates_df2[['Newstart','Newend']]=np.where(New_coordinates_df2['strand_y'].eq('-')[:,None],np.column_stack((New_coordinates_df2['Newstart']-1,New_coordinates_df2['Newend']+0)),New_coordinates_df2[['Newstart','Newend']].values)
                        New_coordinates_df2[['Newstart','Newend']]=np.where(New_coordinates_df2['strand_y'].eq('+')[:,None],np.column_stack((New_coordinates_df2['Newstart']-1,New_coordinates_df2['Newend']+0)),New_coordinates_df2[['Newstart','Newend']].values)
                        New_coordinates_df2=New_coordinates_df2.drop_duplicates(['Name2','Newstart','Newend','strand_y'], keep="first")
                        New_coordinates_df2.to_csv(bedfile_output,header=None,sep="\t",index=False)
    else:
        Interval='FALSE'
        Bedfile.columns=['scaffold','start','end','width','Info','strand']
        print("Generating the bed file...")
        Bedfile['Species']="Species"
        Bedfile['Info']="0"
        Bedfile=Bedfile[['scaffold', 'start','end','Species','Info','strand']]
        Bedfile[['start','end']]=np.where(Bedfile['strand'].eq('-')[:,None],np.column_stack((Bedfile['start']-1,Bedfile['end']+0)),Bedfile[['start','end']].values)
        Bedfile[['start','end']]=np.where(Bedfile['strand'].eq('+')[:,None],np.column_stack((Bedfile['start']-1,Bedfile['end']+0)),Bedfile[['start','end']].values)
        Bedfile=Bedfile.drop_duplicates(['scaffold','start','end','strand'], keep="first") 
        Bedfile.to_csv(bedfile_output,header=None,sep="\t",index=False)
    #Load the fasta file 
    if Interval == 'TRUE':
        from Bio import SeqIO
        d=[]
        if os.path.exists("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa"):
            for seq_record in SeqIO.parse(("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa"), "fasta"):
                d.append({'Scaffold':seq_record.id, 'length':len(seq_record)})
            len_count=pd.DataFrame(d)
        else:
            for seq_record in SeqIO.parse(("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+".fa"), "fasta"):
                d.append({'Scaffold':seq_record.id, 'length':len(seq_record)})
            len_count=pd.DataFrame(d)
        scaff_bis=pd.merge(len_count, New_coordinates_df2, left_on=['Scaffold'],right_on=['Name2'],how='outer')
        if scaff_bis['Scaffold'].isnull().values.any() ==True:
            New_corrected='TRUE'
            print("There have been an issue betwenn the fasta genome file and coordinates, trying to generate a nex fasta file from the coordinate dataframe...")
            data = read_file("/beegfs/data/bguinet/these/Genomes/"+names+"/Mapping/intervals.txt")
            records = prepare_data(data)
            df2 = pd.DataFrame(records)
            chrs = {}
            for seq in SeqIO.parse("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+".fa", "fasta"):
                chrs[seq.id] = seq.seq
            List_scaf = list()
            for i in chrs.keys():
                List_scaf.append(i)    
            with open("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected2.fa", 'w') as out_fasta_corrected_file:
                for index, row in df2.iterrows():
                    if re.search(r"[0-9]_",row['scaffold']):
                        New_scaff_name=row['scaffold']
                        Old_scaff_name=New_scaff_name.rsplit('_', 1)[0]
                        print(">",New_scaff_name,sep='',file=out_fasta_corrected_file)
                        print(textwrap.fill(str(chrs[Old_scaff_name][int(row['start']):int(row['end'])]),width=80),file=out_fasta_corrected_file)
                    else:
                        New_scaff_name=row['scaffold']
                        print(">",New_scaff_name,sep='',file=out_fasta_corrected_file)
                        print(textwrap.fill(str(chrs[New_scaff_name][int(row['start']):int(row['end'])]),width=80),file=out_fasta_corrected_file)
            #Reopen the new corrected fasta file 
            d=[]
            for seq_record in SeqIO.parse(("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected2.fa"), "fasta"):
                        d.append({'Newscaffold':seq_record.id, 'length':len(seq_record)})
            len_count=pd.DataFrame(d)
    # pybedtools version
    print("Printing loci to the species directory")
    #a = pybedtools.example_bedtool("/beegfs/data/bguinet/these/Genomes/"+names+"/Recover_fasta_loci_newcoordinate.bed")
    #fasta = pybedtools.example_filename("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa")
    #b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
    if os.path.isfile("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected2.fa"):
        if Interval=='TRUE':    
            a = pybedtools.example_bedtool(bedfile_output)
            fasta = pybedtools.example_filename("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected2.fa")
            try:
                    b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
            except Exception:
                    pass
    elif  os.path.isfile("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa"):
        if Interval=='TRUE': 
            print("Using the first corrected file")
            a = pybedtools.example_bedtool(bedfile_output)
            fasta = pybedtools.example_filename("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa")
            try:
                    b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
            except Exception:
                    pass
            #b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
        if Interval2=='FALSE': 
            print("ok")
            a = pybedtools.example_bedtool(bedfile_output)
            fasta = pybedtools.example_filename("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_corrected.fa")
            try:
                    b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
            except Exception:
                    pass
    else : 
        print("Nothing has been corrected but we need to create a new assembly file without informations")
        record_dict = SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+".fa", "fasta"))
        with open("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_bis.fa", 'w') as out_fasta_corrected_file:
              SeqIO.write(record_dict.values(), output_handle, 'fasta')
        a = pybedtools.example_bedtool(bedfile_output)
        fasta = pybedtools.example_filename("/beegfs/data/bguinet/these/Genomes/"+names+"/"+names+"_bis.fa")
        try:
              b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
        except Exception:
              pass
        #b = a.sequence(fi=fasta,s=True,fo=outfile,fullHeader=True)
    if os.path.isfile(outfile):
    	print("File create succesfully!")
