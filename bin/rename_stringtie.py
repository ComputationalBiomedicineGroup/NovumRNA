#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import csv
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import re as re
import numpy as np

parser = argparse.ArgumentParser(description='rename stringtie trnscript IDs')
parser.add_argument('--stringtie_out', type=str, help='stringtie GTF')
parser.add_argument('--gtf_out', type=str, help='output file')
parser.add_argument('--vaf_out', type=str, help='output file')
args = parser.parse_args()

string_1 = pd.read_csv(args.stringtie_out, sep='\t', skiprows= 2, header=None)
Chr = args.stringtie_out

def add_chr(Input_1):
    Chr_new = Chr.split("_")[0]
    Name = Chr.split("_")[1]
    my_string = Input_1[8]
    index = my_string.find('STRG')
    final_string = my_string[:index] + Chr_new + "_" + Name + "_" + my_string[index:]
    
    my_string = final_string
    index = my_string.find('"STRG')
    final_string_1 = my_string[:index+1] + Chr_new + "_" + Name + "_" + my_string[index+1:]
    return(final_string_1)

Test = string_1.apply(add_chr, 1)
string_1[8] = Test

string_1.to_csv(args.gtf_out, sep='\t', quoting=csv.QUOTE_NONE, header= None, index = False)

# VAF features

string_1["ID"] = string_1[8].str.split("transcript_id \"").str[1].str.split("\"").str[0]
string_1["ID_short"] = string_1["ID"].str.rsplit('.', 1).str.get(0)
string_1["TPM"] = string_1[8].str.split("TPM \"").str[1].str.split("\"").str[0].astype(float)
string_1["TPM"] = string_1["TPM"].fillna(method='ffill')
string_1["Cov"] = string_1[8].str.split("cov \"").str[1].str.split("\"").str[0].astype(float)
string_1['Cov_transcript'] = np.nan
string_1.loc[string_1[2] == 'transcript', 'Cov_transcript'] = string_1.loc[string_1[2] == 'transcript', 'Cov']
string_1['Cov_transcript'].fillna(method='ffill', inplace=True)
string_1["TPM_iso_max"] = string_1.groupby("ID_short")["TPM"].transform('max')
string_1['isoform_count'] = np.nan
string_1.loc[string_1[2] == 'transcript', 'isoform_count'] = string_1.loc[string_1[2] == 'transcript', 'ID_short'].groupby(string_1['ID_short']).transform('count')
string_1['isoform_count'].fillna(method='ffill', inplace=True)
string_1['isoform_count'] = string_1['isoform_count'].astype(int)
string_1["TPM_iso_vaf"] = round(string_1["TPM"]/string_1["TPM_iso_max"], 1)
string_1["Cov_within_vaf"] = round(string_1["Cov"]/string_1["Cov_transcript"], 1)
string_1["ID"] = 'STRG' + string_1["ID"].str.split("STRG", n=1).str[1]
string_1 = string_1[string_1[2] == "exon"]

string_1.to_csv(args.vaf_out, sep='\t', quoting=csv.QUOTE_NONE, header= None, index = False)



