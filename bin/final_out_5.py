#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import re as re
import numpy as np

parser = argparse.ArgumentParser(description='tranlate to proteins')
parser.add_argument('--BED', type=str, help='')
parser.add_argument('--BED_old', type=str, help='')
parser.add_argument('--Out', type=str, help='output file')
args = parser.parse_args()

Final_anno = pd.read_csv(args.BED, sep='\t', header = None)
Final_5 = pd.read_csv(args.BED_old, sep='\t')
Final_anno = Final_anno.drop([len(Final_anno.columns)-2,len(Final_anno.columns)-3,len(Final_anno.columns)-5,len(Final_anno.columns)-6,len(Final_anno.columns)-7], axis = 1)
names = Final_5.columns.tolist()
names.append("Annotation_2")
names.append("NT_Overlap_2")
Final_anno.columns = names
#Final_anno = Final_anno.drop_duplicates(["Chr", "Peptide_START", "Peptide_STOP", "Transcript"])
Final_anno.to_csv(args.Out, sep="\t", index=False)