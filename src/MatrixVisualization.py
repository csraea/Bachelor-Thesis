import numpy as np
from PIL import Image
import os
import sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import math

def visualize(filename):

    name, ext = os.path.splitext(filename)
    if ext != ".fasta" and ext != ".fa" and ext != ".fna":
        print("Error: specified file is not FASTA format")
        sys.exit()

    covid19 = SeqIO.read(filename, "fasta")
    DNA = covid19.seq

    length_side = math.ceil(math.sqrt(len(DNA)))

    data = np.zeros((length_side, length_side, 3), dtype=np.uint8)

    COLORS = {"A": [255, 0, 0], "T": [0, 0, 255], "G": [0,255,127], "C": [255, 255, 0], "a": [255, 0, 0], "t": [0, 0, 255], "g": [0,255,127], "c": [255, 255, 0]}

    current_base = 0
    for x in range(length_side):
        for y in range(length_side):
            try:
                data[x:x+1, y:y+1] = COLORS[DNA[current_base]]
            except IndexError:
                pass
            current_base += 1
    img = Image.fromarray(data)

    scaled_img = img.resize((length_side*8,length_side*8), 4)
    if not os.path.exists("out"):
        os.mkdir("out")
    os.chdir("out")
    scaled_img.save(os.path.basename(name) + "-Matrix.png")
    os.chdir("..")

