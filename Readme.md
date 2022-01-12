### Downloading viral databases

**NCBI viruses proteins database :** https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=Viruses,%20taxid:10239

* The download involved 07/10/19 **1,471,031 viral sequences** in their protein forms.* *


### Purging of contaminating sequences

Since among the potential candidates, we expected a certain proportion of laboratory or natural contaminants, it was important to exclude as many as possible. 
Thus, we eliminated the bacterium-specific phage family sequences from the dataset. Finally, historically, 
the DNA circles within polydnavirus viral particles have been annotated as being of viral origin and named "polydnavirus",
when in fact they are virulence factors of eukaryotic origin. The only genes of viral origin known to date are those involved in the production of viral particles. 
Therefore, we chose to eliminate these proteins from the dataset to avoid false positives.




--------

## Search for BUSCO genes for each genome.
 This BUSCO step is important in multiple following processes (check assembly quality, build a species phylogeny and build mutliple normal distributions such as 
the expected coverage or G+C% distrubtion along the genomes.

* Snakemake rule: **Busco_analysis**

* Script used : ***Busco v3***

Important files created : 

- **/beegfs/data/bguinet/these/Genomes/{species_names}/run_busco/run_BUSCO_v3/short_summary_{species_names}_BUSCO_v3.txt** (BUSCO summary statistics)
- **/beegfs/data/bguinet/these/Genomes/{species_name}/run_busco/run_BUSCO_v3/compiled_busco_aa** (AA fasta BUSCO sequences)

---------


## Search for viral sequence homology for each genome with MMseqs2 search

![Image description](mmseqs_step.png)

Mmseqs2 query = genome ; db : virus protein sequences

* Snakemake rule: **Homology_analysis**

* Script used : ***Mmseqs2 search***

Important files created : 

- **/beegfs/data/bguinet/these/Genomes/{species_name}/run_mmseqs2/result_mmseqs2.m8** (Mmseqs results)
- **/beegfs/data/bguinet/these/Genomes/{species_name}/run_mmseqs2_V/result_mmseqs2_summary_V.m8** (Candidate oci summary)

*Note : On Lbbe slurm workers, this programm does not run on pbil-27. 

----------

## Characterization of viral and eukaryotic loci 

As a result, we will have a number of viral sequences with homology to certain locus found in Hymenoptera genomes. Among these numerous hits, some candidates may be false positives corresponding to eukaryotic genes present in viruses following a horizontal gene transfer from a eukaryotic to a viral genome, these phenomena being well known in viruses with large genomes (Herniou *et al.,* 2003). In order to eliminate some of these sequences, we use the principle of sequence overlap between viral hits and BUSCO genes (a priori of eukaryotic origin). Thus, the presence of a viral hit at the same location in the genome as that of a BUSCO gene probably indicates a viral capture event of this eukaryotic gene by a viral genome. Using the R code ```Overlapping_sequences.R```, we thus define 3 types of locus, each being defined by overlapping several MMseqs2 hits: (1) viral locus sharing no overlap with a BUSCO gene, (2) Hymenoptera BUSCO loci sharing no overlap with a viral hit and finally (3) a viral/hymenoptera BUSCO locus having an overlap of at least 10% of the size of the viral hit. Thereafter, only the non-overlapping viral loci (category 1) are retained since they are good candidates for endogenization. Finally, during the MMseq2 search, a single viral protein may present homologies at several locations within a scaffold, corresponding to several HSPs (High Scoring Pairs). We concatenate these HSPs in the alignment part of the pipeline (see [HSPs analysis](#HSPs-analysis)).



#### 1) Adding strands and changing the direction of the coordinates in the file result_mmseqs2.m8

![Image description](mmseqs2_output_step)


#### 2) For refseq mmseqs2 table
```
cat /beegfs/data/bguinet/these/Species_genome_names.txt | while read line; do python3 Make_change_strand_mmseqs2.py -b /beegfs/data/bguinet/these/Genomes/${line}/run_mmseqs2_V/result_mmseqs2.m8 -o /beegfs/data/bguinet/these/Genomes/${line}/run_mmseqs2_V -t virus; done
```

#### 3) For BUSCO tblastn tab (because busco already did a tblastn research, we will take it from the source)
```
To delete >  cat /beegfs/data/bguinet/these/Species_genome_names.txt | while read line; do python3 Make_change_strand_mmseqs2.py -b /beegfs/data/bguinet/these/Genomes/${line}/run_busco/run_BUSCO_v3/blast_output/tblastn_${line}_BUSCO_v3.tsv -o /beegfs/data/bguinet/these/Genomes/${line}/run_busco/run_BUSCO_v3/blast_output/ -t hymenoptera; done


to keep 
import pandas as pd 

list_of_names1=[]
for names in open("/beegfs/data/bguinet/these/Genomes/Species_genome_names_test.txt","r"):
	list_of_names1.append(names)

list_of_names2=[]

for names in list_of_names1:
	list_of_names2.append(names.replace("\n", ""))
  
for sp_name in list_of_names2:
  print(sp_name)
  df1 = pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+sp_name+"/run_busco/run_BUSCO_v3/full_table_"+sp_name+"_BUSCO_v3.tsv",comment='#',sep="\t")
  df1['Strand']='NA'
  row=0
  for Busco_id, Contig in zip (df1['Busco_id'],df1['Contig']):
    if df1.at[row,'Status']=="Missing":
      row+=1
      continue
    else:
      try: 
        df2 = pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+sp_name+"/run_busco/run_BUSCO_v3/augustus_output/predicted_genes/"+Busco_id+".out.1",comment='#',sep="\t",header=None)
      except: 
        try:
            df2 = pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+sp_name+"/run_busco/run_BUSCO_v3/augustus_output/predicted_genes/"+Busco_id+".out.2",comment='#',sep="\t",header=None)
        except:
            df2 = pd.read_csv("/beegfs/data/bguinet/these/Genomes/"+sp_name+"/run_busco/run_BUSCO_v3/augustus_output/predicted_genes/"+Busco_id+".out.3",comment='#',sep="\t",header=None)
      Strand=df2[6][0]
      df1.at[row, 'Strand'] = Strand
      row+=1
      print(row," / ",len(df1['Busco_id']))
  df1.to_csv("/beegfs/data/bguinet/these/Genomes/"+sp_name+"/run_busco/run_BUSCO_v3/full_table_"+sp_name+"_BUSCO_v3_strand.tsv",header=True,sep="\t")
    
    
    
    
```

Generate a file : **result_mmseqs2_strand_V.m8** de type :
```
"query", "tlen", "target", "pident", "alnlen", "mismatch", "gapopen","qstart", "qend", "tstart", "tend", "evalue", "bits", "strand"
```

![Image description](Overlapping_step.png)

#### 4) Switch to R, in order to launch R on the server use */beegfs/data/soft/R-3.5.2/bin/R
```
/beegfs/data/soft/R-3.5.2/bin/Rscript Overlapping_sequences_BUSCO_Viral_loci.R /beegfs/data/bguinet/these/Species_genome_names.txt /beegfs/data/bguinet/these/Genomes/ #Permet d'executer le fichier . R
```

Creation of :

**Matches_Apis_mellifera_strand_V.m8** : un fichier avec les brains d'affichés ainsi que les coordonnées changées de sens
**Matches_",i,"_summary_V.txt** : a file of type "seqnames" "start" "end" "width" "strand" "type" in which the HSPs are grouped together and have a new coordinate.



Recovery of all loci according to their coordinates in the genomes of each species in a single file. **All_fasta_viral_loci.fna**

![Image description](All_non_overlapping_candidate_step.png)

```
bash Recover_loci_sequences2.sh /beegfs/data/bguinet/these/Genomes/ /beegfs/home/bguinet/these_scripts_2/ /beegfs/data/bguinet/these/Species_genome_names.txt
```

### Filtering of loci using NR database. 

We will remove loci matching with Bacterial of Insects protein sequences fellowing some criterias.


* Snakemake file : **Homologous_snakemake**
* Snakemake rule: **Filter_loci**

* Script used : ***mmsesq createdb, search, convertalis,Filter_loci_with_NR.py***

Important file created :  **/beegfs/data/bguinet/these/Results/Candidate_viral_loci_filtred.m8** (used in Filter candidate)


## Gene clustering analysis 

We will perform now clustering with the idea of gathering candidate loci and viral proteins in the same clusters by sequence homologies. 
For this we proceed to 5 steps, the main part of the pipeline code is written in the file **Snakefile_clustering_analysis**. 

### Create fasta files with virus and filtred loci to feed the clustering method 

- We will concatenate the filtred queries loci with the proteins from the target (virus proteins database)

```
cat /beegfs/data/bguinet/these/Viral_sequences_loci/All_fasta_viral_loci_filtred.aa /beegfs/data/bguinet/these/NCBI_protein_viruses/All_viral_protein_sequences_without_contamination_controls.fa >> /beegfs/data/bguinet/these/Results/Candidate_viral_loci_and_viral_protein.aa
```

### Connected component Clustering

Note: We use the following parameters: 
      #--cluster-mode 1 : Connected component clustering that uses the transitive connection to cover more distant homologs.
      #--cov-mode 0 : Bidirectional coverage, where only sequences with overlapping sequence lengths greater than 30% of the longer of the two sequences are clustered, 
      #(in this case always the viral sequence since the loci are defined by viral hits). 
      #-evalue 0.0001 : To eliminate false positives. 
      # --cluster-reassign : During the cascade clustering of Mmseqs2, as the representative of a cluster can change at each iteration, 
      #it can happen that some members that were already close to a cluster do not fulfill the clustering criteria anymore. We therefore correct this by reassigning 
      #these sequences

* Snakemake file: **/Snakefile_Clustering_analysis**
* Snakemake rule: **Clustering**

* Script used : ***Mmseqs2 cluster***

Important file created  : 
- **/beegfs/data/bguinet/these/Clustering3/Candidate_viral_loci_and_viral_protein_clu.tab** (Clustering table)
- **/beegfs/data/bguinet/these/Clustering3/Filtred_Candidate_viral_loci_and_viral_protein_clu.tab** (Filtred clustering table without orphelin clusters or with cluster with > 50 different Hymenoptera species)
--------
