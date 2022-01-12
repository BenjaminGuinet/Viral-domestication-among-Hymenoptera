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
# This BUSCO step is important in mulitple following processes (check assembly quality, build a species phylogeny and build mutliple normal distributions such as 
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



