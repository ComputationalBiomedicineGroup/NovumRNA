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
parser.add_argument('--Final', type=str, help='')
parser.add_argument('--Out_novel_sb', type=str, help='Final_out')
parser.add_argument('--Out_novel_wb', type=str, help='output file')
parser.add_argument('--Out_diff_sb', type=str, help='output file')
parser.add_argument('--Out_diff_wb', type=str, help='output file')
args = parser.parse_args()


Final = pd.read_csv(args.Final, sep='\t')
filter_col = [col for col in Final if col.startswith('Rank')]
Final["Min_Rank"]= Final[filter_col].min(axis = 1)
Final_diff = Final[(Final["Annotation"] == "Differential") & (Final["Transcript_ref"] == "Protein")]
Final_novel = Final[(Final["Annotation"] != "Differential") & (Final["Transcript_ref"] == "Protein")]
Final_novel_wb = Final_novel[(Final_novel["Min_Rank"] <= 2) & (Final_novel["Min_Rank"] > 0.5)]
Final_novel_sb = Final_novel[Final_novel["Min_Rank"] <= 0.5]
Final_diff_wb = Final_diff[(Final_diff["Min_Rank"] <= 2) & (Final_diff["Min_Rank"] > 0.5)]
Final_diff_sb = Final_diff[Final_diff["Min_Rank"] <= 0.5]

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def Scoring(data):
    data["TPM_score"] = NormalizeData(data.TPM)
    data["E_Coverage_score"] = NormalizeData(data.E_Coverage)
    data["BAM_reads_score"] = NormalizeData(data.BAM_reads)
    data["Min_Rank_score"] = NormalizeData(data.Min_Rank)
    data["Min_Rank_score"] = 1 - data["Min_Rank_score"]
    data["Score"] = (data["TPM_score"] + data["E_Coverage_score"] + 
                         data["BAM_reads_score"] + data["Min_Rank_score"])/4
    data["Score_2"] = (data["TPM_score"] + data["E_Coverage_score"] + 
                         data["BAM_reads_score"])/3
    data["Score"] = round(data["Score"], 2)
    #data = data.drop(["TPM_score", "E_Coverage_score", "BAM_reads_score", "Min_Rank_score"], axis = 1)
    data = data.sort_values(by="Score", ascending=False)
    return data

Final_novel_wb_score = Scoring(Final_novel_wb)
Final_novel_sb_score = Scoring(Final_novel_sb)
Final_diff_wb_score = Scoring(Final_diff_wb)
Final_diff_sb_score = Scoring(Final_diff_sb)

Final_novel_wb_score.to_csv(args.Out_novel_wb, sep="\t", index=False)
Final_novel_sb_score.to_csv(args.Out_novel_sb, sep="\t", index=False)
Final_diff_wb_score.to_csv(args.Out_diff_wb, sep="\t", index=False)
Final_diff_sb_score.to_csv(args.Out_diff_sb, sep="\t", index=False)