#!/usr/bin/python

import os
import sys
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
import numpy as np

def save(outFileName, plt):
    if not os.path.exists("out"):
        os.mkdir("out")
    os.chdir("out")
    plt.savefig(os.path.basename(outFileName) + "-ORFs.png")
    os.chdir("..")

def visualize(filename, mode):
    name, ext = os.path.splitext(filename)
    if ext != ".gbk" and ext != ".genbank" and ext != ".gb":
        print("Error: specified file is not of GenBank format")
        sys.exit()

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(17, 9), sharex=True, gridspec_kw={"height_ratios": [4, 1]}
    )
    ax1.set_title(os.path.basename(name) + ': Open reading frames', size=22, weight='bold')

    # PLOT THE RECORD MAP
    record = SeqIO.read(filename, "genbank")
    gRecord = BiopythonTranslator().translate_record(record)
    gRecord.plot(ax=ax1, with_ruler=False, strand_in_label_threshold=4)

    # PLOT THE LOCAL GC CONTENT (using 50bp windows)
    gcContent = lambda s: 100.0 * len([c for c in s if c in "GC"]) / 50
    xxx = np.arange(len(record.seq) - 50)
    yyy = [gcContent(record.seq[x : x + 50]) for x in xxx]
    ax2.fill_between(xxx + 25, yyy, alpha=0.3)
    ax2.set_ylim(bottom=0)
    ax2.set_ylabel("GC(%)")
    
    save(name, plt)

