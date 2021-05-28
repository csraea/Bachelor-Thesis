from PIL import ImageDraw
import PIL
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import os
import sys

w = 0
h = 0
x1 = 0
y1 = 0

def save(outFileName, image):
    
    if not os.path.exists("out"):
        os.mkdir("out")
    os.chdir("out")
    image.save(os.path.basename(outFileName) + "-Gates.png")
    os.chdir("..")

def visualize(filename, mode, start, end):
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
        start = 0
    if end > len(DNA) or end == 0:
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

    image1 = PIL.Image.new("RGB", (int((A+T)/4+1000), int((G+C)/4+1000)), (255, 255, 255))
    draw = ImageDraw.Draw(image1)

    y1 = int((G+C)/8+500)
    x1 = int((A+T)/8+500)
    
    progressY = 0
    progressX = 0
    color = None

    for n in truncated:
        if n == 'A' or n == 'a':
            progressY = -1
            progressX = 0
            color = "red"
        if n == 'T' or n == 't':
            progressY = 1
            progressX = 0
            color = "blue"
        if n == 'G' or n == 'g':
                progressY = 0
                progressX = -1
                color = "green"
        if n == 'C' or n == 'c':
                progressY = 0
                progressX = 1
                color = "yellow"
        if n != '\n':
                draw.line([x1, y1,  x1 + progressX, y1 + progressY], fill = color, width = 1)
                x1 += progressX
                y1 += progressY

    save(name, image1)
