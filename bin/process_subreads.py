#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re as re

parser = argparse.ArgumentParser(description='tranlate to proteins')
parser.add_argument('--tsv_file', type=str, help='')
parser.add_argument('--tsv_out', type=str, help='output file')
parser.add_argument('--BAM_cov', type=int, help='output file')
args = parser.parse_args()

BED = pd.read_csv(args.tsv_file, sep = "\t", engine='python', header = None, error_bad_lines=False)

BED["SAMPLE"] = BED[3].str.split("STRG").str[1].fillna(method='ffill')
BED["SAMPLE"] = "STRG" + BED["SAMPLE"]
BED["Sequence"] = BED[0].str.split(" ").str[0]
BED["Annotation"] = BED[5].fillna(method='ffill')
BED["Chromosome"] = BED[0].str.split("chr").str[1].fillna(method='ffill')
BED["Strand"] = BED[6].fillna(method='ffill')
BED["Chromosome"] = "chr" + BED["Chromosome"]
BED["START"] = BED[1].fillna(method='ffill').astype(int)
BED["STOP"] = BED[2].fillna(method='ffill').astype(int)
BED = BED[BED["Sequence"].str.contains("chr")==False]
BED["Count"] = BED.groupby("Sequence")["Sequence"].transform('count')
BED =  BED.iloc[:,8:]
idx = BED.groupby(['Sequence'])['Count'].transform(max) == BED['Count']
BED = BED[idx]
BED_rmdup = BED.drop_duplicates("Sequence")
BED_rmdup = BED_rmdup[BED_rmdup["Sequence"].str.contains(">")==False]
BED_rmdup = BED_rmdup[BED_rmdup["Sequence"].str.contains("-")==False]
BED_rmdup["Count"].value_counts()
Test_1 = BED_rmdup.iloc[:,[3,5,6]]
Test_2 = BED_rmdup.iloc[:,[0,1,4,2,7]]
C = Test_1.join(Test_2)
C["Total_count"] = C.groupby("SAMPLE")["Count"].transform('sum')
idx = C.groupby(['SAMPLE'])['Count'].transform(max) == C['Count']
C = C[idx]
C = C.drop_duplicates("SAMPLE")
C = C[~C.Sequence.str.contains(r'[0-9]')]
C = C[C.Sequence.str.isupper()]
C = C[C["Total_count"] > args.BAM_cov]

def reverse_strand(Input):
    if Input[5] == "-":
        Input[4] = str(Seq(Input[4]).reverse_complement())
    return(Input)

D = C.apply(reverse_strand, 1)


D.to_csv(args.tsv_out, sep="\t", index=False, header=False)