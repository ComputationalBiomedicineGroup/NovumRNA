#!/usr/bin/env python3

import pandas as pd
import argparse
import io
from itertools import permutations

parser = argparse.ArgumentParser()
parser.add_argument('--HLA_HD_out', type=str)
parser.add_argument('--file_out', type=str)
args = parser.parse_args()

data_file = args.HLA_HD_out

# Delimiter
data_file_delimiter = '\t'

# The max column count a line in the file could have
largest_column_count = 0

# Loop the data lines
with open(data_file, 'r') as temp_f:
    # Read the lines
    lines = temp_f.readlines()

    for l in lines:
        # Count the column count for the current line
        column_count = len(l.split(data_file_delimiter)) + 1

        # Set the new most column count
        largest_column_count = column_count if largest_column_count < column_count else largest_column_count

# Close file
temp_f.close()

# Generate column names (will be 0, 1, 2, ..., largest_column_count - 1)
column_names = [i for i in range(0, largest_column_count)]

# Read csv
DRB = pd.read_csv(data_file, header=None, delimiter=data_file_delimiter, names=column_names, skipfooter = 8, skiprows = 3, engine='python')
# This list is from the command line tool option -list, it gives the avilable alleles to not run in errors later
HLA = pd.read_csv("/scripts/Valid_HLAII_alleles.txt", sep = ",", header=None)

DRB[1] = DRB[1].str.replace('HLA-', '').str.rsplit(':', n=1).str[0]
DRB[2] = DRB[2].str.replace('HLA-', '').str.rsplit(':', n=1).str[0]

DRB_1 = DRB[1][(DRB[1] != 'Not typed') & (DRB[1] != '-')]
DRB_2 = DRB[2][(DRB[2] != 'Not typed') & (DRB[2] != '-')]

DRB_both = DRB_1.tolist() + DRB_2.tolist()

allowed_prefixes = ['DQB1', 'DQA1', 'DPA1', 'DPB1']
allowed_data = [item for item in DRB_both if any(prefix in item for prefix in allowed_prefixes)]
rest_data = [item for item in DRB_both if not any(prefix in item for prefix in allowed_prefixes)]

combinations_list = ["-".join(perm) for perm in permutations(allowed_data, 2)]
DRB_both = combinations_list + rest_data

if DRB_both:
    Final_list = HLA[HLA[0].isin(DRB_both)]
else:
    Final_list = pd.DataFrame(data=[])

Final_list.to_csv(args.file_out, header= None, index = False)