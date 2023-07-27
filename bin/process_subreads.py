#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio.Seq import Seq
import re as re
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--tsv_file', type=str)
parser.add_argument('--tsv_out', type=str)
parser.add_argument('--BAM_cov', type=int)
args = parser.parse_args()

data = []

# Read the file line by line and process each line
with open(args.tsv_file, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        # Check the number of values in the row
        num_values = len(row)
        
        # If the row has only 1 value, add NaN to the data list
        if num_values == 1:
            data.append([row[0], pd.NA])
        # If the row has 8 values, add all of them to the data list
        elif num_values == 8:
            data.append(row)
        # For any other number of values, raise an error or handle accordingly
        else:
            continue

# Create a DataFrame from the processed data list
BED = pd.DataFrame(data)

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