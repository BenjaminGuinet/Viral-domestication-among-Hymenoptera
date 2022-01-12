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

* Snakemake rule: **Snakemake_homologous_analysis**

* Script used : ***Mmseqs2 search & convertalis***

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
- **/beegfs/data/bguinet/these/Results/Candidate_viral_loci_filtred_clu.m8"** (Clustering table)
- **/beegfs/data/bguinet/these/Results/Candidate_viral_loci_filtred_clu_filtred.m8** Filtred clustering table without orphelin clusters or with cluster with > 50 different Hymenoptera species)
- **/beegfs/data/bguinet/these/Results/Candidate_viral_loci_filtred_clu_filtred_aligned.m8** (Filtrer Clustering table with alignment informations)

--------


###########################

## Candidate Gene informations

All the Rules and codes are written in the Snakefile : **Snakefile_add_informations**

### TAX ID information

Retrieving the **taxonomic** and **ontological** informations provided by each of the viral genes candidates to a viral domestication allows us to know more about its origin. We can thus retrieve very valuable information such as the origin of the donor virus, its family, for which product it codes, in which cell compartment it is expressed, etc. For this we query SQlite databases for taxonomy and the UniProt package present in python for gene ontology. 

#### Prerequisites

##### 1) Create sqlite database using taxadb in order to recover the taxa informations

```
pip install taxadb
```

##### 2) Get TAxid of all target protein seqs and add it on the column in the merged silix_blast dataframe

```
/beegfs/data/bguinet/myconda/bin/taxadb download -o taxadb
/beegfs/data/bguinet/myconda/bin/taxadb create -i taxadb --dbname taxadb.sqlite --fast
```

* Snakemake Rule : **Add_taxonomy**

* Script used : **Add_taxid_info.py**

Important file created : **/beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid.m8**

Column added : 

- Taxid
- no_rank 
- subgenus
- genus                         
- species
- subfamily                        
- family
- superfamily 
- infraorder
- suborder
- order
- infraclass
- subclass
- class
- subphylum
- phylum
- kingdom
- superkingdom


The process took **13 minutes** for **26318 taxids.**

The number of cluster is : xxx
The number of viral protein is : xxx

------------

### Gene Ontology information

The function of the candidats genes in wasps can be explored via the presence or absence of key viral genes and their expression. Whereas functional characterization of virus genes is mostly lacking, a wealth of information is available about the functions of homologous genes within the other closly related virus genes. By extension, the functions of the knwon virus genes can be used to predict the functions of the candidat genes.


* Snakemake rule : **Add_ontology**

* Script used : **Add_genomic_ontology.py**

Important file created : **/beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology.m8**

Column added : 

- genomic_structure    
- Taxonomic.lineage.IDs
- Entry                             
- Organism
- Gene.names
- Protein.names
- Gene.ontology..GO
- Gene.ontology..biological.process
- Gene.ontology..molecular.function
- Gene.ontology..cellular.component
- Gene.ontology.IDs
- Status

The process took **54 minutes** for **26319 sequences.**

The number of cluster is : xxx
The number of viral protein is : xxx

------------

### HMMER analysis on candidat loci

Each candidate locus should, if of viral origin, return particular domains known to be present in free viruses when tested via HMMER. This analysis will thus allow us to (i) confirm the viral character of a candidate locus and (ii) study the presence or absence of known domains within the loci, telling us more about the putative function of this locus beyond a simple ontological assignment by taking over the one assigned to a homologous viral protein. 

hmmscan scores a database of HMMs vs. a database of sequences. An HMM database is what needs to be processed by hmmpress.

#### Prerequisites
First we need to dowload the pfam db :

```
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.hmm
```

Then we create the hmmer profile db 

```
/beegfs/data/bguinet/TOOLS/hmmer-3.2.1/src/hmmpress Pfam-A.hmm
```

Generate a fasta file with all Candidate loci filtred and Split the loci into multiple file to accelerate the HMMprocess :
```
mkdir /beegfs/data/bguinet/these/Domain_analysis2/Splited_fasta_loci/
cd /beegfs/data/bguinet/these/Domain_analysis2/Splited_fasta_loci/
awk -v size=500 -v pre=prefix -v pad=5 '/^>/{n++;if(n%size==1){close(f);f=sprintf("%s.%0"pad"d",pre,n)}}{print>>f}' /beegfs/data/bguinet/these/Clustering3/Candidate_loci_filtred.aa 
```

* Snakemake Rules : **HMMER_analysis1, HMMER_analysis2**
(we keep only domain best hits based on bic score and remove all hits with evalue < 0.0005 and cov > 0.35).

* Script used : **Hmmscan,HmmPy.py*

Important file created : **/beegfs/data/bguinet/these/Domain_analysis2/Splited_fasta_loci/Hmmscan_loci_analysis_parsed_{number}.tab**

We will then use the taxonomic information for each candidate loci (when available) in the **Add_genomic_env.py** script. 

------------------

### Merge HMMER and 2LCA analysis

* Snakemake rule : **Merge_2LCA_HMMER_analysis**

* Script used : **Add_2LCA_HMMER.py**

Important file created : **/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER.m8**

Column added : 

- Domain_name
- Pfam_acc
- Hmmscan_evalue
- Domain_description
- genus_2LCA
- family_2LCA
- order_2LCA
- kingdom_2LCA
- superkingdom_2LCA
- kingdom_2LCA

___________________

## Evaluation of viral integration through the study of the Genomic environment  


![Image description](Contaminants.jpeg)

Horizontal gene transferst (HGT) is defined as the accidental acquisition of DNA from another species independently of reproduction and evolutionary distance. To date, the majority of cases of horizontal gene transferst have been found in bacteria and unicellular eukaryotes (Koonin, 2016), as these organisms feed on the phagocytosis of bacteria and other microorganisms, resulting in frequent contact with prokaryotic DNA, which may predispose them to incorporate genetic material into their genomes. Conversely, relatively few cases of functional HGT in eukaryotes have been detected or studied in detail. Indeed, HGT events appear to be rarer in eukaryotes since, (i) DNA is protected by a nucleus, which makes contact between exogenous genetic elements and the genome more difficult, and (ii) in metazoan species at least, the germ line is usually separated from the rest of the cells, so a gene acquired in a metazoan has no chance of being fixed in a population if it does not reach the germ line (Danchin, 2016). Despite these constraints, genes conferring a function have been acquired by HGTs in both vertebrates and invertebrates (e.g. in insects (Moran & Jarvik, 2010), nematodes (Scholl et al., 2003), fish (Sun et al., 2015), terrestrial plants (Quispe-Huamanquispe, Gheysen & Kreuze, 2017) or humans (Crisp et al., 2015)). However, several initial reports that drew much attention to HGTs in humans could not be confirmed by further examination (Crisp et al., 2015), which casts doubt on the true extent of HGTs. Indeed, contamination of DNA samples or low quality genome assembly can easily lead to apparent HGA.

Thus, among potential candidates for horizontal transfer followed by domestication, there may be false positives. These false positives can then correspond to two things. (i) It may be a natural contamination, i.e. we find a viral sequence not integrated in the eukaryotic genome but present in the eukaryotic cell. It will then probably be a free-living virus present in the cell for which its genome has been sequenced and then assembled during sequencing. This type of contaminant is nevertheless of interest since it potentially allows us to know that a virus infects the eukaryote under study, then it is also possible to reconstruct its viral genome from the viral scaffolds. (ii) It may be a laboratory contaminant, i.e. during the laboratory handling steps, external viral particles may have infected the samples. 

One way to make the hypothesis of viral endogenization of these candidate genes more likely is to study the genomic environment around these genes.

The sequencing data available allows us to obtain coverage along the genomes. For this we perform a mapping using the Hisat2 software. The idea being to discriminate the genomic environments of scaffolds, starting from the principle that each genomic entity has a tendency to have a cover and a composition in its own base, corresponding for the first to the abundance of this compartment in the extract of DNA, and for the second the average number of G + C bases. Thus, it is possible to sort the scaffolds according to their coverage and their base composition. Here we distinguish the scaffolds comprising a depth and a base composition of the “nuclear gene” type as being those which show no significant difference compared to the coverage or a base composition observed by the busco genes.

A non-parametric pvalue is calculated, for this we sampled 100 loci of the size of the scaffold of interest within the set of typical BUSCO eukaryotic scaffolds (Chimera). These 100 samples represent a distribution of the typical coverage expected under the hypothesis H0 that it is a eukaryotic scaffold. The pvalue then corresponds to the proportion of BUSCO coverage values ​​more extreme than that observed within the candidate scaffold (bilateral test). We reject the H0 hypothesis if it is a eukaryotic scaffold with a 5% probability of error.

Another way to asses wether one viral candidat gene is likely endogenised or not into the eucaryote genome is to look among the known Hymenoptera genes such as the BUSCO genes, which are found near candidate genes, that is to say on the same scaffold. Such an observation allows us to think that it is indeed a locus integrated in the genome of the Hymenoptera.

![Image_description](Scaffold_origin_assessment.png)

_______________________

## Eucaryotic genes in scaffolds

If a scaffold containing genes of viral origin that are candidates for domestication while also containing genes that are typically eukaryotic, then this gives an additional clue to the endogenous aspect of the gene in question since it will not support the hypothesis that the scaffold is indeed of eukaryotic origin. 

For this analysis we will therefore first need to predict the genes in the genomes of our species. 

For this we will use an *ab initio* approach of the **Augustus** type which will work using the training parameters of the previous BUSCO analyses for these genomes. 

### 1- Gene prediction allong genomes**

#### 1-1- Extract the scaffolds containing candidates** 

* Snakemake file : **Snakefile_add_Augustus_repeat_analysis**
* Snakemake rule : **Extract_scaffolds*

* Script used : **Extract_candidat_scaffolds.py**

Important file created : **/beegfs/data/bguinet/these/Genomes/{species_names}/All_scaffold_containing_{species_names}_candidates.fna** (All scaffolds fasta sequences containing)

#### 1-2- Run augustus on each fasta files** 

* Snakemake file : **Snakefile_add_Augustus_repeat_analysis**
* Snakemake rule : **Augustus_analysis*

* Script used : **Augustus3.3*

Important file created : **/beegfs/data/bguinet/these/Genomes/{species_names}/run_augustus/run_augustus_new.out** (All augustus predictions)


### 2- Taxonomic assignation for predicted augustus genes**

* Snakemake file : **Snakefile_add_Augustus_repeat_analysis**
* Snakemake rule : **MMseqs_augustus*

* Script used : **Mmseqs2 search*


Important file created : **/beegfs/data/bguinet/these/Augustus_search_analysis/Augustus_Uniprot_tax_result.tsv** (All Nr hits for each augustus predictions)

______________________

## Repeat elements in scaffolds

Since few sequenced viral genomes have ETs, or rarely more than one, then the probability that a viral gene flanked by an ET (which would have jumped into the viral genome beforehand) will integrate into the host genome is low. (Gilbert et al 2017). It is thought that this is because the majority of ET integrations in viral genomes are deleterious and are therefore eliminated by negative selection.  This is even lower if 2 ETs are found in the flanks of the EVEs and/or if it can be shown that the ET(s) present in the flanks are present elsewhere in the host genome.

* Snakemake file : **Snakefile_add_Augustus_repeat_analysis**
* Snakemake rule : **Repeat_analysis*

* Script used : **Mmseqs2 search*

Important file created : **/beegfs/data/bguinet/these/Repeat_env_analysis/All_scaffold_repeat_search_result_dna_strand_V_reduced.m8** (All repeat hits within scaffols containing candidates)

______________________

## Merge Repeat and Augustus results

We then merge Augustus and repeat results by counting the number of **Arthopoda/insecta** genes and **repeat elements** within scaffolds containing candidates: 

* Snakemake file : **Snakefile_add_Augustus_repeat_analysis**
* Snakemake rule : **Merge_augustus_repeat_file*

* Script used : **Add_augustus_repeat_to_blast.py*

Important file created : **/beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat.m8** (All eucaryotic and repeat count added into the file)
columnnames : **count_eucaryote; count_repeat**

Column added : 

- count_eucaryote
- count_repeat

_________________________

### Scaffolds coverage and GC content 

#First we will have to generate for each genome a chimeric sequence wich simply corresponds to one chimeric sequence composed of all the scaffold where have been found at least once a Busco sequence.

#The idea is to create a chimeric huge fasta file for each genome species. Each chimeric sequence will be the sum of all scaffolds that had at least one BUSCO sequence into it (which means that this scaffold is likely from eucaryot origin)
#To do so we will get the BUSCO summary table generated by the tblasn analysis, this will give us the scaffold that have at least on BUSCO.
#Because we can aldo find virus loci into the same scaffold as BUSCO, will also have to remove from these scaffolds the putative viral loci.
#To do so we will get the viral summary table generated by the mmseqs2 analysis, this will five us the coordinates of the viral loci into each scaffold.
#We will use these coordinates in order to remove them from the scaffolds of interest.

* Snakemake file : **Snakefile_Add_genomic_env_analysis**
* Snakemake rule : **Add_genomic_env**

* Script used : **Get_busco_scaff.py ,Add_genomic_env_last.py**

Important file created : **/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_env.m8**

Column added : 

- GC_content_BUSCO
- GC_content_scaffold
- pvalue_gc    
- FDR_pvalue_gc (if True, candidate GC significantly != busco GC)
- Number_busco_loci
- Number_viral_loci
- cov_depth_BUSCO
- cov_depth_candidat
- pvalue_cov
- FDR_pvalue_cov (if True, candidate coverage significantly != busco coverage)


The process took **54 minutes** for **26319 sequences.**

The number of cluster is : 3972
__________________________________
