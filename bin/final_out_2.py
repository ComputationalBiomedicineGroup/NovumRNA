#!/usr/bin/env python3

import argparse
import pandas as pd
import re as re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--Filtering', type=str)
parser.add_argument('--VAF', type=str)
parser.add_argument('--BED', type=str)
parser.add_argument('--Translation', type=str)
parser.add_argument('--BED_2', type=str, help='')
parser.add_argument('--BIND', type=str)
parser.add_argument('--Specific', type=str)
parser.add_argument('--Out', type=str)
parser.add_argument('--Out_header', type=str)
args = parser.parse_args()

Filtering = pd.read_csv(args.Filtering, sep='\t', header=None)
Filtering.rename(columns={1:'Peptide'}, inplace=True)
Filtering["Annotation"] = Filtering[0].str.split("_").str[3].str.split(" ").str[0]
Filtering["Transcript"] = Filtering[0].str.split("_").str[0]
Filtering["Peptide_START"] = Filtering[0].str.split("_").str[1].astype(int)
Filtering["Peptide_STOP"] = Filtering[0].str.split("_").str[2].astype(int)
Filtering = Filtering.drop([0, 2], axis = 1)

VAF = pd.read_csv(args.VAF, sep='\t', header=None)
VAF = VAF[[3,4,6,9,15,16,17,22]]
VAF_Filtering = pd.merge(how = "inner", left=Filtering, right=VAF, left_on=["Transcript"], right_on=[9])
VAF_Filtering["Remove"] = np.where((VAF_Filtering["Peptide_START"] >= VAF_Filtering[3]) & (VAF_Filtering["Peptide_STOP"] <= VAF_Filtering[4]), "Fine", "Remove")
VAF_Filtering = VAF_Filtering[VAF_Filtering["Remove"] != "Remove"]
VAF_Filtering = VAF_Filtering.drop(["Remove", 3, 4, 6, 9], axis = 1)
VAF_Filtering.columns = ["Peptide", "Annotation", "Transcript", "Peptide_START", "Peptide_STOP", "isoform_count", "TPM_iso_perc", "Cov_within_VAF", "Gene"]

BED = pd.read_csv(args.BED, sep='\t')
BED = BED[["Chr", "E_START", "E_STOP", "STRAND", "Annotation", "NT_Overlap", "Overlap_perc", "ID", "E_Coverage", "TPM"]]
BED["Transcript"] = BED["ID"].str.split("_").str[2]
BED = BED.drop("ID", axis = 1)

BED_Filtering = pd.merge(how = "inner", left=VAF_Filtering, right=BED, left_on=["Transcript", "Annotation"], right_on=["Transcript", "Annotation"])
BED_Filtering["Remove"] = np.where((BED_Filtering["Peptide_START"] >= BED_Filtering["E_START"]) & (BED_Filtering["Peptide_STOP"] <= BED_Filtering["E_STOP"]) & (BED_Filtering["STRAND"] == "+")
         | (BED_Filtering["Peptide_START"] <= BED_Filtering["E_STOP"]) & (BED_Filtering["Peptide_STOP"] >= BED_Filtering["E_START"]) & (BED_Filtering["STRAND"] == "-"), "Fine", "Remove")
BED_Filtering = BED_Filtering[BED_Filtering["Remove"] != "Remove"]
BED_Filtering = BED_Filtering.drop("Remove", axis = 1)

Translation = pd.read_csv(args.Translation, sep='\t')
Translation = Translation[["ID", "Translation_ref"]]
Translation = Translation.drop_duplicates()
Translation["Transcript"] = Translation["ID"].str.split("_").str[2]
Translation["ID"] = Translation["ID"].str.split('_STRG').str[0]

BED_Filtering_Translation = pd.merge(how = "inner", left=BED_Filtering, right=Translation, left_on="Transcript", right_on="Transcript")

BED_2 = pd.read_csv(args.BED_2, sep='\t', header=None)

def switch(Input):
    if Input[5] == "-":
        a = Input[2]
        Input[2] = Input[1]
        Input[1] = a
    return(Input)

BED_2 = BED_2.apply(switch, axis=1)

BED_2[3] = "STRG" + BED_2[3].str.split("STRG").str[1].str.split("_").str[0]

BED_2.columns = ["Chr", "Peptide_START", "Peptide_STOP", "Transcript","Peptide_nt", "Strand", "Annotation", "BAM_reads", "BAM_reads_all"]
BED_2 = BED_2[["Peptide_START", "Peptide_STOP", "Transcript", "Annotation", "BAM_reads", "BAM_reads_all"]]

BED_Filtering_Translation_2 = pd.merge(left=BED_Filtering_Translation, right=BED_2, left_on=["Transcript", "Peptide_START", "Peptide_STOP", "Annotation"], 
                         right_on=["Transcript", "Peptide_START", "Peptide_STOP", "Annotation"], how='left')

BIND = pd.read_csv(args.BIND, sep='\t')
BIND = BIND[["HLA Allele", "Best Percentile", "Epitope Seq"]]
# Use pivot to make 'HLA Allele' values as columns
pivoted_df = BIND.groupby(['Epitope Seq', 'HLA Allele'])['Best Percentile'].min().unstack(fill_value=0).reset_index()
# If you want to fill NaNs with 0
pivoted_df = pivoted_df.fillna(0)
pivoted_df.columns = ['Rank_' + str(col) if col != 'Epitope Seq' else col for col in pivoted_df.columns]
pivoted_df = pivoted_df.rename(columns={'Epitope Seq': 'Peptide'})
BIND_2 = pivoted_df.drop_duplicates("Peptide")

BED_Filtering_Translation_2_BIND = pd.merge(left=BED_Filtering_Translation_2, right=BIND_2, left_on=["Peptide"], 
                         right_on=["Peptide"], how='left')

Specific = pd.read_csv(args.Specific, sep='\t')
Specific = Specific[["Full", "REF", "REF_NT", "NEW_NT", "Mismatch"]]
Specific = Specific.drop_duplicates("Full")

Final = pd.merge(left=BED_Filtering_Translation_2_BIND, right=Specific, left_on=["Peptide"], 
                         right_on=["Full"], how='left')

Final["Mismatch"] = Final["Mismatch"].str.replace(r"[\[\]()' ]", "", regex=True)

def process_mismatch(row):
    if not row or row.strip() == '':
        return '', '', ''
    parts = row.split(',')
    # Splitting parts into three groups
    SNP, Ref_NT, SNP_pos = [], [], []
    for i in range(0, len(parts), 3):
        # Add checks to prevent index errors
        SNP.append(parts[i] if i < len(parts) else '')
        Ref_NT.append(parts[i+1] if i+1 < len(parts) else '')
        SNP_pos.append(parts[i+2] if i+2 < len(parts) else '')
    return ','.join(SNP), ','.join(Ref_NT), ','.join(SNP_pos)

Final[['SNP', 'Ref_NT', 'SNP_pos']] = Final['Mismatch'].apply(process_mismatch).tolist()

Final = Final.dropna()
Final_2 = Final[["Chr", "Peptide_START", "Peptide_STOP", "Peptide", "Transcript", "STRAND", "ID", "SNP", "Ref_NT", "SNP_pos"]]
Final_3 = Final.drop(["Chr", "Peptide_START", "Peptide_STOP", "Transcript", "ID", "Peptide", "STRAND", "Full", "Mismatch", "SNP", "Ref_NT", "SNP_pos"], axis = 1)
Final_4 = pd.concat([Final_2, Final_3], axis=1)

def switch_2(Input):
    if Input["STRAND"] == "-":
        a = Input[2]
        Input[2] = Input[1]
        Input[1] = a
    return(Input)

Final_5 = Final_4.apply(switch_2, axis=1)

rank_columns = Final_5.filter(like='Rank')
Final_5 = Final_5[(rank_columns <= 2).any(axis=1)]

Final_5.to_csv(args.Out, sep="\t", index=False, header = None)
Final_5.to_csv(args.Out_header, sep="\t", index=False)