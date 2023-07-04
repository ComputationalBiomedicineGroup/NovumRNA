#!/usr/bin/env python3

import pandas as pd
import csv
import argparse

parser = argparse.ArgumentParser(description='tranlate to proteins')
parser.add_argument('--GTF', type=str, help='peptides from whippet downstream')
parser.add_argument('--ID', type=str, help='peptides from whippet downstream')
parser.add_argument('--Out', type=str, help='output file')
args = parser.parse_args()

gtf_file = args.GTF
IDs = pd.read_csv(args.ID, sep = "\t", engine='python', header = None, error_bad_lines=False)

test_list = IDs[0].tolist()
df = pd.read_csv(gtf_file, sep = "\t", engine='python', header = None)
df["transcript"] = df[8].str.split("transcript_id \"").str[1].str.split("\"").str[0]
subset = df[df["transcript"].str.contains('|'.join(test_list))]
subset = subset.drop("transcript", axis=1)
subset.to_csv(args.Out, sep='\t', header= None, index = False, quoting=csv.QUOTE_NONE)