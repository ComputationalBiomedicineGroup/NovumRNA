#!/usr/bin/env python3

import pandas as pd
import argparse
import io

parser = argparse.ArgumentParser(description='Refine HLA-HD output')
parser.add_argument('--HLA_HD_out', type=str, help='HLA-HD output')
parser.add_argument('--file_out', type=str, help='output file')
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
HLA = pd.read_csv("/home/ausserh/projects/2021/CRCnoncanonical/nonCanonicalNeoAG/possible_MHCII_alleles.txt", sep = ",", header=None)

DRB[1] = DRB[1].str.replace('HLA-', '')
DRB[1] = DRB[1].str.replace('*', '')
DRB[1] = DRB[1].str.replace(':', '')
DRB[2] = DRB[2].str.replace('HLA-', '')
DRB[2] = DRB[2].str.replace('*', '')
DRB[2] = DRB[2].str.replace(':', '')

DRB_1 = DRB[1][(DRB[1] != 'Not typed') & (DRB[1] != '-')]
DRB_2 = DRB[2][(DRB[2] != 'Not typed') & (DRB[2] != '-')]

DRB_both = DRB_1.tolist() + DRB_2.tolist()

if DRB_both:
    y = []
    for i in DRB_both:
        if "DRB" in i :
            a = i[:4] + '_' + i[4:8]
        else:
            a = i[:8]
        y.append(a)

    df = pd.DataFrame(data=y)
    
    y = []
    for i in df[0]:
        if HLA[0].str.contains(i).any():
            Value = HLA[HLA[0].str.contains(i)]
            y.append(Value[0].iloc[0])

    df = pd.DataFrame(data=y)
else:
    df = pd.DataFrame(data=[])

df.to_csv(args.file_out, header= None, index = False)