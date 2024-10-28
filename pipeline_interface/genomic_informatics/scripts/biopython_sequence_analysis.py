# Script for Biopython sequence analysis
from Bio import SeqIO

def parse_sequences(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"ID: {record.id}")
        print(f"Sequence: {record.seq}")

if __name__ == "__main__":
    parse_sequences('../data/cancer_sequences.fasta')
