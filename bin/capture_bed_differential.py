#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='tranlate to proteins')
parser.add_argument('--GTF', type=str, help='')
parser.add_argument('--Out', type=str, help='')
args = parser.parse_args()

mTEC_all = pd.read_csv(args.GTF, sep='\t', header=None)

mTEC_all["ID"] = mTEC_all[8].str.split("transcript_id \"").str[1].str.split("\"").str[0]
mTEC_all['TPM'] = mTEC_all[8].str.split("TPM \"").str[1].str.split("\"").str[0].astype(float)
mTEC_all['TPM'] = mTEC_all['TPM'].fillna(method='ffill')
mTEC_all_exon = mTEC_all[mTEC_all[2] == "exon"]
mTEC_all_exon["Coverage"] = mTEC_all_exon[8].str.split("cov \"").str[1].str.split("\"").str[0].astype(float)
mTEC_all_bed = mTEC_all_exon[[0,3,4,"TPM","Coverage",6,"ID"]]
mTEC_all_bed["Count"] = 1
mTEC_all_bed.to_csv(args.Out, sep="\t", index=False, header=None)