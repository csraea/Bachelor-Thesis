from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import sys

def getStats(filename, mode, start, end):

    flag = False

    name, ext = os.path.splitext(filename)
    if ext != ".fasta" and ext != ".fa" and ext != ".fna":
        print("Error: specified file is not FASTA format")
        sys.exit()

    covid19 = SeqIO.read(filename, "fasta")
    DNA = covid19.seq

    if start >= end and start != 0:
        print("Error: START is greater then END")

        if mode == 'v':
            return
        else:
            sys.exit()

    if start > len(DNA):
        flag = True
        start = 0
    if end > len(DNA) or end == 0:
        flag = True
        end = len(DNA)

    truncated = DNA[start:end]

    A = 0
    T = 0
    G = 0
    C = 0

    for n in truncated:
        if n == 'A' or n == 'a':
            A += 1
        elif n == 'T' or n == 't':
            T += 1
        elif n == 'G' or n == 'g':
            G += 1
        elif n == 'C' or n == 'c':
            C += 1
        else:
            # skip unrecognized symbol
            pass

    if mode == 'v' and flag:
        print("Frequencies of nucleotides on the interval [" + str(start) + ";" + str(end) + "]:\n")
    else:
        print("Frequencies of nucleotides:")

    print("A:\t" + str(A))
    print("T:\t" + str(T))
    print("G:\t" + str(G))
    print("C:\t" + str(C))
    print("Total:\t" + str(A+T+G+C))
    
    if mode == 'v' and flag:
        print("GC-content on interval [" + str(start) + ";" + str(end) + "]:\t" + str(round((G+C)/(A + T + G + C), 4)) + "%")
    else:
        print("GC-content:\t" + str(round((G+C)/(A + T + G + C), 4)) + "%")

    return
