import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Print out a message when the program is initiated.
#print('-----------------------------------------------------------------------------------\n')
#print('                        ###Translate DNA to Amino Acide#######.\n')
#print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-f", "--fasta_file", help="the fasta nucleotide file")
parser.add_argument("-o", "--out_file", help="the output directory")
args = parser.parse_args()

dnaseq= args.fasta_file
Out_file=args.out_file


# List of functions:
record_makers = []
with open(Out_file, 'w') as aa_fa:
        for dna_record in SeqIO.parse(dnaseq, 'fasta'):
                if len(dna_record.seq) > 130:
                        # Translate the sequence
                        aa_seqs = dna_record.translate()
                        # write new record
                        print('>',dna_record.id, file= aa_fa,sep="")
                        print(aa_seqs.seq,file=aa_fa)

print("The translated file has been saved to :", Out_file)
