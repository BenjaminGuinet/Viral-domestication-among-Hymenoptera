import pandas as pd 
import os
list_of_names1=[]
for names in open("/beegfs/data/bguinet/these/Species_genome_names.txt","r"):
        list_of_names1.append(names)
list_of_names2=[]
for names in list_of_names1:
        list_of_names2.append(names.replace("\n", ""))


table=pd.read_csv("/beegfs/data/bguinet/these/Transciptomic_analysis/Table_RNA_reads3.txt",sep=";")

Genome_size_table = pd.DataFrame(columns=("Species_name","Genome_size"))
for species in list_of_names2:
        sub_table=table.loc[table['ScientificName'].str.contains(species)]
        try:
            table2= pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+str(species)+"/Genome_assembly_statistics/report.tsv",header=None,sep="\t")
            Genome_size_table =Genome_size_table.append({"Species_names":table2[1][0],"Genome_size":table2[1][7]}, ignore_index=True)
        except:
        		print(row['ScientificName'], " not in data")


#Whole Body

table_whole=table.loc[table['Tissue'].str.contains('whole', na=False) | table['Tissue'].str.contains('Whole', na=False)]
table_whole=table_whole.loc[~table_whole['Tissue'].str.contains('head', na=False)]
table_whole=table_whole.loc[~table_whole['Tissue'].str.contains('brain', na=False)]
table_whole=table_whole.loc[~table_whole['LibrarySource'].str.contains('META', na=False)]

table=table_whole

new_df = pd.merge(table, Genome_size_table,  how='right', left_on=['ScientificName'], right_on = ['Species_names'])
new_df= new_df.drop_duplicates()

new_df["Genome_size"] = pd.to_numeric(new_df["Genome_size"])
new_df['Nb_reads_coverage_X']= new_df['Genome_size']*300
new_df['Nb_reads_coverage_X'] = pd.to_numeric(new_df['Nb_reads_coverage_X'])
new_df = new_df[new_df['Run'].notna()]
new_df['Perc_coverage']= (new_df['bases'].astype(int)*300)/new_df['Nb_reads_coverage_X'].astype(int)

new_df=new_df.sort_values(['Species_names', 'Perc_coverage'], ascending=[True, False])
new_df['Sex'] = pd.Categorical(new_df['Sex'], ['female','pooled male and female', 'male and female','male'])
new_df=new_df.sort_values("Sex")
new_df['cumsum'] = new_df['Perc_coverage'].groupby(new_df['Species_names']).cumsum()
new_df['Perc_cumulative'] = new_df.groupby(['Species_names'])['Perc_coverage'].apply(lambda x: x.cumsum())

#Remove SRA that have already reached the coverage of X.
new_df=new_df[new_df.Perc_cumulative.gt(300).groupby(new_df.Species_names).cumsum()<=1].copy()




for species in list_of_names2:
		sub_table=new_df.loc[new_df['ScientificName'].str.contains(species)]
		if sub_table.empty:
			print("No SRA GENOMIC data")
		else: 
			with open("/beegfs/data/bguinet/these/Transciptomic_analysis2/Transcriptomic_jobs2/Mapping_job_"+str(species)+"_Transcriptomic.sh",'w') as w:
				w.write("#!/bin/bash\n")
				w.write("#SBATCH --mem 10G\n")
				w.write("#SBATCH -t 72:00:00\n")
				w.write("#SBATCH --cpus-per-task=10\n")
				w.write("#SBATCH -e /beegfs/data/bguinet/these/Transciptomic_analysis2/Transcriptomic_jobs2/Mapping_job_"+str(species)+"_Transcriptomic.error\n")
				w.write("#SBATCH -o /beegfs/data/bguinet/these/Transciptomic_analysis2/Transcriptomic_jobs2/Mapping_job_"+str(species)+"_Transcriptomic.out\n")
				w.write("#SBATCH -J Mapping_job_"+str(species)+"\n")
				w.write("\n")				
				w.write("#Download reads\n")			
				List_SRA_number=[]
				for sra in sub_table['Run']:
						List_SRA_number.append(sra)
				for sra in List_SRA_number:
						print(str(species)+"_"+str(sra))
						#Download all reads 
						w.write("\n")
						w.write("rm -rf /beegfs/data/bguinet/sra_reads/sra/"+str(sra)+"*\n")
						w.write("/beegfs/data/bguinet/TOOLS/sratoolkit.2.11.0-ubuntu64/bin/prefetch --max-size 100000000 "+str(sra)+" && /beegfs/data/bguinet/TOOLS/sratoolkit.2.11.0-ubuntu64/bin/fasterq-dump -t /beegfs/data/bguinet/these/Genomes/"+str(species)+"/Mapping/Transcriptomic_reads/ --threads 10 -f "+str(sra)+" -O /beegfs/data/bguinet/these/Genomes/"+str(species)+"/Mapping/Transcriptomic_reads/\n")
						w.write("/beegfs/data/bguinet/Bguinet_conda/bin/pigz --best /beegfs/data/bguinet/these/Genomes/"+str(species)+"/Mapping/Transcriptomic_reads/"+str(sra)+"*\n")
						List_SRA_number=[]
						List_SRA_number1=[]
						List_SRA_number2=[]
						List_SRA_single=[]
						for sra in sub_table['Run']:
							if sub_table.loc[sub_table['Run'].str.contains(sra)]['LibraryLayout'].iloc[0]=="SINGLE":
								List_SRA_single.append(sra+".fastq.gz")
							elif sub_table.loc[sub_table['Run'].str.contains(sra)]['LibraryLayout'].iloc[0]=="PAIRED":
								List_SRA_number.append(sra)
						for i in List_SRA_number:
							i1=i+"_1.fastq.gz"
							List_SRA_number1.append(i1)
							i2=i+"_2.fastq.gz"
							List_SRA_number2.append(i2)
				w.write("date;hostname;pwd\n")
				w.write("###############################################################\n")
				w.write("#       INPUT for species :"+str(species)+"\n")
				w.write("###############################################################\n")
				w.write("# ASSEMBLY FILE :\n")
				w.write("mkdir /beegfs/data/bguinet/these/Genomes/"+str(species)+"/Mapping/Transcriptomic_reads \n")
				if os.path.exists("/beegfs/data/bguinet/these/Genomes/"+species+"/"+species+"_corrected2.fa"):
					w.write("FILE=/beegfs/data/bguinet/these/Genomes/"+str(species)+"/"+str(species)+"_corrected2.fa\n")                
				elif os.path.exists("/beegfs/data/bguinet/these/Genomes/"+species+"/"+species+"_corrected.fa"):
					w.write("FILE=/beegfs/data/bguinet/these/Genomes/"+str(species)+"/"+str(species)+"_corrected.fa\n")
				elif os.path.exists("/beegfs/data/bguinet/these/Genomes/"+species+"/"+species+"_bis.fa"):
					w.write("FILE=/beegfs/data/bguinet/these/Genomes/"+str(species)+"/"+str(species)+"_bis.fa\n")            
				w.write("OUT=/beegfs/data/bguinet/these/Genomes/"+str(species)+"/Mapping/Transcriptomic_reads\n")
				w.write("# READS to map onto assembly\n")
				w.write("\n")
				w.write("###############################################################\n")
				w.write("#       SOFTWARES                                             #\n")
				w.write("##############################################################\n")
				w.write("\n")
				w.write("# MAPPING HISAT2\n")
				w.write("HISAT2=/beegfs/data/bguinet/Bguinet_conda/bin\n")
				w.write("#bedtools\n")
				w.write("BEDTOOLS=/beegfs/data/bguinet/Bguinet_conda/bin\n")
				w.write("# samtools\n")
				w.write("SAMTOOLS=/beegfs/data/soft/samtools-1.9/bin\n")
				w.write("\n")
				w.write("###############################################################\n")
				w.write("#        MAP SHORT READS TO INFER COVERAGE DEPTH              \n")
				w.write("###############################################################\n")
				w.write("cd $OUT\n")
				w.write("\n")
				w.write("# build index\n")
				w.write("$HISAT2/hisat2-build $FILE mapping_index -p 10 \n")
				w.write("# mapping\n")
				if len(List_SRA_single)>=1 and len(List_SRA_number1) >= 1:
					w.write("$HISAT2/hisat2 --dta -U "+','.join(str(item) for item in List_SRA_single)+" -k 1 -p 10 -q -x mapping_index -1 "+','.join(str(item) for item in List_SRA_number1)+" -2 "+','.join(str(item) for item in List_SRA_number2)+"  | $SAMTOOLS/samtools view -o mapping_"+str(species)+".bam 2> stats_mapping.txt\n")
				if len(List_SRA_single)>=1 and len(List_SRA_number1) == 0 :
					w.write("$HISAT2/hisat2 --dta -U "+','.join(str(item) for item in List_SRA_single)+" -k 1 -p 10 -q -x mapping_index | $SAMTOOLS/samtools view -o mapping_"+str(species)+".bam 2> stats_mapping.txt\n")
				if len(List_SRA_single)==0 and len(List_SRA_number1) >= 1 :
					w.write("$HISAT2/hisat2 --dta  -k 1 -p 10 -q -x mapping_index -1 "+','.join(str(item) for item in List_SRA_number1)+" -2 "+','.join(str(item) for item in List_SRA_number2)+"  | $SAMTOOLS/samtools view -o mapping_"+str(species)+".bam 2> stats_mapping.txt\n")
				w.write("\n")
				w.write("# TPM calculation \n")
				w.write("\n")
				w.write("#First create the GFF file \n")
				w.write("python3 /beegfs/home/bguinet/these_scripts_2/Make_GFF_file.py -s "+str(species)+"\n")
				w.write("\n")
				w.write("#Second run the TPMcalculator algorithm\n")
				w.write("/beegfs/data/bguinet/anaconda3/envs/my_python3.6_env/bin/TPMCalculator -g /beegfs/data/bguinet/these/Genomes/"+str(species)+"/Mapping/Transcriptomic_reads/GFF_candidates.gff -b /beegfs/data/bguinet/these/Genomes/"+str(species)+"/Mapping/Transcriptomic_reads/mapping_"+str(species)+".bam -e -a\n")
				w.write("\n")
