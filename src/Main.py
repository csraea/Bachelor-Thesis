#!/usr/bin/python

import sys
import argparse
import os
from os import listdir
from os.path import isfile, join

import SeqCollector
import StatGenerator
import GatesVisualization
import MatrixVisualization
import HMatrixVisualization
import ORFPlotter
import Comparison

DATA_PATH = "./data/"
FASTA = "SARS-CoV-2.fasta"
GENBANK = "SARS-CoV-2.gbk"

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--mode", help="Execution mode: quiet / verbose", default="v", choices=['q', 'v'], dest="mode", required=False)
parser.add_argument("-d", "--download", help="Download SARS-CoV-2 genome associated files", action="store_true", dest="download", required=False)
parser.add_argument("-g", "--gates", help="Perform Gates' visualization. Parameter is an input sequence filename", dest="gates", required=False)
parser.add_argument("-o", "--orf", help="Plot ORFs of the genome. Parameter is an input sequence filename", dest="orf", required=False)
parser.add_argument("-s", "--stat", help="Obtain genome statistical data, including the distribution of nucleotides and a GC-content", action="store_true", dest="stat", required=False)
parser.add_argument("-x", "--matrix", help="Plot the nucleotide sequence into a 2D matrix. Parameter is the input sequence filename", dest="matrix", required=False)
parser.add_argument("-i", "--hmatrix", help="Plot the nucleotide sequense into a Hashed 2D matrix. Parameter is the output figure filename", dest="hash", required=False)
parser.add_argument("-c", "--compare", help="Compare specified genome sequences using the pairwise2 algorithm", dest="comp", nargs=2, required=False)
parser.add_argument("-S", "--size", help="Size of the picture side in pixels", type=int, dest="size", required=False)
parser.add_argument("-p", "--pos", help="Start and end positions of the nucleotide sequence to perform an action (0 for defaults)", type=int, nargs=2, dest="pos", required=False)
parser.add_argument("-a", "--all", help="Perform all possible actions but comprasion in the default mode", action="store_true", dest="all", required=False)
parser.add_argument("-n", "--name", help="Input sequnce filename", dest="name", required=False)

args = parser.parse_args()

def verifyArgs():
    # if args.mode == 'v' and (args.all is True or args.download is True or args.gates is not None or args.orf is not None or args.prot is not None or args.stat is True or args.matrix is not None or args.hash is not None or args.comp is not None or args.size is not None or args.output is not None or args.pos is not None):
    #     parser.error("verbose mode does not support other command line arguments")
    #     sys.exit()
    # elif args.download and (args.all):
    #     parser.error("--all already includes --download")
    #     sys.exit()

    if args.stat and args.pos is None:
        parser.error("--stat requires -P (--pos) argument")
        sys.exit()
    if args.gates and args.pos is None:
        parser.error("--gates requires -P (--pos) argument")
        sys.exit()
    if args.hash and args.size is None:
        parser.error("--hmatrix requires -S (--size) argument")
        sys.exit()
    pass 

def welcomeBanner():
    print("+----------------------------------+")
    print("|--- Welcome to the Visualizer! ---|")
    print("+----------------------------------+")

def mainMenu():
    print(color.BOLD + "Choose the option: " + color.END)
    print("1. Download SARS-CoV-2 genome sequence & associated files")
    print("2. Plot sequence statistics")
    print("3. Gates' visualization")
    print("4. 2D Matrix visualization")
    print("5. Improved 2D Matrix visualization")
    print("6. Plot ORFs")
    print("7. Compare genomes")
    print("8. Exit")

    while 1 :
        opt = input(color.BOLD + "Choice: " + color.END)
        if opt.isdigit() and int(opt) > 0 and int(opt) < 9:
            return int(opt)
        elif opt == 'q' or opt == 'Q' or opt == 'quit' or opt == 'QUIT':
            return 8
        else:
            print("Incorrect option!")

def vObtainFiles(msg):
    seqExist = False

    print(color.BOLD + msg + color.END)
    onlyfiles = [f for f in listdir(DATA_PATH) if isfile(join(DATA_PATH, f))]
    idx = 1
    for f in onlyfiles:
        name, ext = os.path.splitext(DATA_PATH + f)
        if ext == ".fasta" or ext == ".fa" or ext == ".fna":
            print(str(idx) + ". " + f)
            seqExist = True
            idx += 1

        
    if not seqExist:
        print("No genome sequences are available!")
        return None

    while 1 :
        opt = input(color.BOLD + "Choice: " + color.END)
        if opt.isdigit() and int(opt) > 0 and int(opt) < idx:
            break
        elif opt == 'q' or opt == 'Q' or opt == 'quit' or opt == 'QUIT':
            sys.exit()
        else:
            print("Incorrect option!")
    idx = 1
    for f in onlyfiles:
        name, ext = os.path.splitext(DATA_PATH + f)
        if idx == int(opt) and (ext == ".fasta" or ext == ".fa" or ext == ".fna"):
            return f
        elif ext == ".fasta" or ext == ".fa" or ext == ".fna":
            idx += 1
        else: pass

    return None


def vObtainFiles2(msg):
    seqExist = False

    print(color.BOLD + msg + color.END)
    onlyfiles = [f for f in listdir(DATA_PATH) if isfile(join(DATA_PATH, f))]
    idx = 1
    for f in onlyfiles:
        name, ext = os.path.splitext(DATA_PATH + f)
        if ext == ".gbk" or ext == ".genbank" or ext == ".gb":
            print(str(idx) + ". " + f)
            seqExist = True
            idx += 1

        
    if not seqExist:
        print("No genome sequences are available!")
        return None

    while 1 :
        opt = input(color.BOLD + "Choice: " + color.END)
        if opt.isdigit() and int(opt) > 0 and int(opt) < idx:
            break
        elif opt == 'q' or opt == 'Q' or opt == 'quit' or opt == 'QUIT':
            sys.exit()
        else:
            print("Incorrect option!")
    idx = 1
    for f in onlyfiles:
        name, ext = os.path.splitext(DATA_PATH + f)
        if idx == int(opt) and (ext == ".gbk" or ext == ".genbank" or ext == ".gb"):
            return f
        elif ext == ".gbk" or ext == ".genbank" or ext == ".gb":
            idx += 1
        else: pass

    return None

def vObtainInterval():
    print("Specify the interval (0 for the entire genome)")
    
    while 1:
        start = input("Start:\t")
        if start.isdigit() and int(start)  >= 0:
            end = input("End:\t")
            if end.isdigit() and int(end) >= 0:
                break
            else:
                print("Error: position cannot be less then 0.")
        else:
            print("Error: position cannot be less then 0.")
    return int(start), int(end) 
    
def vObtainSize():
    print("Specify the size of an image side (preferably a power of 2, >= 512):")
    
    while 1:
        size = input("Size (px):\t")
        if size.isdigit() and int(size)  >= 512:
            return int(size)
        else:
            print("Error: invalid size.")

def main():
    verifyArgs()
    
    opt = None
    if args.mode == 'v':
        welcomeBanner()
        opt = mainMenu()


    while 1:
        # 1 Download contents
        if opt == 1 or args.download or args.all:
            if args.mode == 'v':
                print("Necessary files are being downloaded...")
        
            SeqCollector.downloadFiles()

            if args.mode == 'v':
               print("Done!")

        # 2 Plot sequence statistics
        if opt == 2 or args.stat or args.all:
            filename = None

            start = -1
            end = -1

            if args.mode == 'v':
                filename = vObtainFiles("Choose the sequence to plot the statistics of:")
                start, end = vObtainInterval()
            elif args.all:
                filename = FASTA
                start = 0
                end = 0
            else:
                if args.pos is not None:
                   start = args.pos[0]
                   end = args.pos[1]
                filename = args.name

            if filename is not None:
                StatGenerator.getStats(DATA_PATH + filename, args.mode, start, end)
                pass
        if opt == 3 or args.gates or args.all:
            filename = None

            start = -1
            end = -1

            if args.mode == 'v':
                filename = vObtainFiles("Choose the sequence to visualize using Gates' method:")
                start, end = vObtainInterval()
            elif args.all:
                filename = FASTA
                start = 0
                end = 0
            else:
                if args.pos is not None:
                   start = args.pos[0]
                   end = args.pos[1]
                filename = args.gates

            if filename is not None:
                GatesVisualization.visualize(DATA_PATH + filename, args.mode, start, end)
                if args.mode == 'v':
                    print("Done")
            pass
        if opt == 4 or args.matrix or args.all:
            filename = None

            if args.mode == 'v':
                filename = vObtainFiles("Choose the sequence to visualize using 2D Matrix method:")
            elif args.all:
                filename = FASTA
            else:
                filename = args.matrix

            if filename is not None:
                MatrixVisualization.visualize(DATA_PATH + filename)
                if args.mode == 'v':
                    print("Done")
            pass
        if opt == 5 or args.hash or args.all:
            filename = None
            size = 512
            if args.mode == 'v':
                filename = vObtainFiles("Choose the sequence to visualize using 2D HMatrix method:")
                size = vObtainSize()
            elif args.all:
                filename = FASTA
            else:
                filename = args.hash
                size = args.size
            if filename is not None:
                HMatrixVisualization.visualize(DATA_PATH + filename, args.mode, size)
                if args.mode == 'v':
                    print("Done")
            pass
        if opt == 6 or args.orf or args.all:
            filename = None
            if args.mode == 'v':
                filename = vObtainFiles2("Choose the annotation file to visualize ORFs:")
            elif args.all:
                filename = GENBANK
            else:
                filename = args.orf
            if filename is not None:
                ORFPlotter.visualize(DATA_PATH + filename, args.mode)
                if args.mode == 'v':
                    print("Done")
            pass
        if opt == 7 or args.comp:
            filename1 = None
            filename2 = None

            if args.mode == 'v':
                filename1 = vObtainFiles("Choose the first sequence to compare:")
                filename2 = vObtainFiles("Choose the second sequence to compare:")
            else:
                filename1 = args.comp[0]
                filename2 = args.comp[1]

            if filename is not None:
                Comparison.compare(DATA_PATH + filename1, DATA_PATH + filename2, args.mode)
                if args.mode == 'v':
                    print("Done")
            pass
        if opt == 8:
            sys.exit()

        
        if args.mode == 'v':
            input()
            opt = mainMenu()
        else:
            sys.exit()

if __name__ == "__main__":
    main()