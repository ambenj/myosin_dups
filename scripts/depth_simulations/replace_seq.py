# Run with:
# python3 replace_seq.py <ref_chr> <ref_start> <ref_end> <ref.fa> <replacement_seq.fa> <out.fa>

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

# user input
CHR = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])
REF_FA = sys.argv[4]
NEW_SEQ_FA = sys.argv[5]
OUT = sys.argv[6]

# read in new sequence (must be a fasta file with only one record)
new_seq_record = SeqIO.read(open(NEW_SEQ_FA, mode='r'), "fasta")

# Parse reference genome
with open(OUT, 'w') as f_out:
    for seq_record in SeqIO.parse(open(REF_FA, mode='r'), 'fasta'):
        # remove .id from .description record (remove all before first space)
        seq_record.description=' '.join(seq_record.description.split()[1:])
        print('SequenceID = '  + seq_record.id)
        print('Description = ' + seq_record.description + '\n')

        # Look for target sequence
        if seq_record.id == CHR:
            print("Found", CHR)

            # Change sequence at specified position at sequence
            seq_record.seq = seq_record.seq[:start] + new_seq_record.seq + seq_record.seq[end:]
            print("New seq at", CHR, start, end,":", seq_record.seq[start:start+20], "...")
            print("Replacement seq starts with", new_seq_record.seq[0:20], "...")

        # write new fasta file
        r=SeqIO.write(seq_record, f_out, 'fasta')
        if r!=1: print('Error while writing sequence:  ' + seq_record.id)
