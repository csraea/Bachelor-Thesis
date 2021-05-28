from Bio import pairwise2
import os
import sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def compare(filename1, filename2, mode):
    name1, ext = os.path.splitext(filename1)
    if ext != ".fasta" and ext != ".fa" and ext != ".fna":
        print("Error: specified file is not FASTA format")
        sys.exit()

    sequence = SeqIO.read(filename1, "fasta")
    genome1 = sequence.seq


    name2, ext = os.path.splitext(filename2)
    if ext != ".fasta" and ext != ".fa" and ext != ".fna":
        print("Error: specified file is not FASTA format")
        sys.exit()

    sequence = SeqIO.read(filename1, "fasta")
    genome2 = sequence.seq

    if mode == 'v':
        print('Sequence lengths:')
        print(os.path.basename(name1) + ':', len(genome1.seq))
        print(os.path.basename(name2) + ':', len(genome2.seq))

    # Alignments using pairwise2 alghoritm
    res = pairwise2.align.globalxx(genome1, genome2, one_alignment_only=True, score_only=True)
    print(os.path.basename(name1) +  "/" + os.path.basename(name2) + 'similarity (%):', res / len(genome1) * 100)