#!/bin/bash
if [ $# -ne 3 ]
then
   echo "Please provide in the first argument the directory where are found all the Species directories"
   echo "Please provide in the second argument the directory where to found the file Create_bed_file.py"
   echo "Please provide in the third argument the the directory where are found all the Species directories"
   echo 'Exemple usage: bash Recover_loci_sequences.sh /beegfs/data/bguinet/these/Genomes/ /beegfs/home/bguinet/these_scripts/ /beegfs/data/bguinet/these/Species_genome_names.txt'
   exit 1
fi


cat $3 | while read line; do python3 $2Create_bed_file.py --table $1${line}/run_mmseqs2_V/result_mmseqs2_summary_V.m8 --out $1${line}/Recover_fasta_loci.bed; echo $'\n';
done


cat $3 | while read line; do bedtools getfasta -fi $1${line}/${line}.fa -bed $1${line}/Recover_fasta_loci.bed -s -fo $1${line}/Fasta_viral_loci_seq_${line}.fa ; sed -i 's@)@):'${line}'@g' $1${line}/Fasta_viral_loci_seq_${line}.fa; echo "Loci sequences recovered and the output is written at : $1${line}/Fasta_viral_loci_seq_${line}.fa";
done
echo $'\n'

#Allows to get fasta sequences of viral loci with seq name such as :
#>scaffold_number:start-end(strand):Species_name

# -s Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.


#Now we will concatenate all viral loci sequences together

echo "${1%/*/}/"Viral_sequences_loci
echo "${1%/*/}/"
mkdir "${1%/*/}/"Viral_sequences_loci
echo $'Copying of files into a new file in order to contatenate them: \n'
cat $3 | while read line; do cp -v $1${line}/Fasta_viral_loci_seq_${line}.fa "${1%/*/}/"Viral_sequences_loci/Fasta_viral_loci_seq_${line}.fa; done
echo $'\n'



cat "${1%/*/}/"Viral_sequences_loci/Fasta_viral_loci_seq_* > "${1%/*/}/"Viral_sequences_loci/All_fasta_viral_loci.fna
echo "Concatenate succeed, ouput file written at : "${1%/*/}/"Viral_sequences_loci/All_fasta_viral_loci.fna"


