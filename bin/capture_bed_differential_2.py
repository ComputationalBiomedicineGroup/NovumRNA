#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--BED', type=str)
parser.add_argument('--Out', type=str)
parser.add_argument('--tpm_max_diff', type=int)
parser.add_argument('--cov_max_diff', type=int)
args = parser.parse_args()


Merged = pd.read_csv(args.BED, sep='\t', header=None)

Merged_low = Merged[(Merged[4] <= args.tpm_max_diff) & (Merged[7] <= args.cov_max_diff)]
Merged_low_export = Merged_low[[0,1,2,3,4,5]]

Merged_low_export[3] = "Differential"

Merged_low_export.to_csv(args.Out, sep='\t', header=None, index=False)