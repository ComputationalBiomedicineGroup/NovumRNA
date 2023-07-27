#!/usr/bin/env python3

import pandas as pd
import argparse
import mapply
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--BED', type=str)
parser.add_argument('--GTF', type=str)
parser.add_argument('--Out_anno', type=str)
parser.add_argument('--Out_bed', type=str)
parser.add_argument('--TPM_min_diff', type=int)
parser.add_argument('--Cov_min_diff', type=int)
parser.add_argument('--TPM_min_novel', type=int)
parser.add_argument('--Cov_min_novel', type=int)
args = parser.parse_args()

GTEx_bed = pd.read_csv(args.BED, sep='\t', header=None)
GTEx_bed_transcripts = GTEx_bed[GTEx_bed[2] == "transcript"]
GTEx_bed_exon = GTEx_bed[(GTEx_bed[2] == "exon") | (GTEx_bed[2] == "CDS")]
GTEx_bed_exon["Length"] = GTEx_bed_exon[11] - GTEx_bed_exon[10]
GTEx_bed_exon["Overlap_perc"] = GTEx_bed_exon[15]/(GTEx_bed_exon["Length"]/100)
GTEx_bed_exon = GTEx_bed_exon[GTEx_bed_exon[15] > 23]
GTEx_bed_exon["Cov"] = GTEx_bed_exon[8].str.split("cov \"").str[1].str.split("\"").str[0].astype(float)
GTEx_bed_exon["ID"] = GTEx_bed_exon[8].str.split("transcript_id \"").str[1].str.split("\"").str[0]
GTEx_bed_transcripts["ID"] = GTEx_bed_transcripts[8].str.split("transcript_id \"").str[1].str.split("\"").str[0]
GTEx_bed_transcripts['Counts'] = GTEx_bed_transcripts.groupby('ID')['ID'].transform('count')

GTEx_bed_transcripts["Lenght"] = GTEx_bed_transcripts[4] - GTEx_bed_transcripts[3]
GTEx_bed_transcripts['Sum'] = GTEx_bed_transcripts.groupby('ID')['Lenght'].transform('sum')
GTEx_bed_transcripts["Single"] = ((GTEx_bed_transcripts["Sum"] - GTEx_bed_transcripts["Lenght"]) == GTEx_bed_transcripts["Lenght"])
GTEx_bed_transcripts = GTEx_bed_transcripts[GTEx_bed_transcripts["Single"] == False]
GTEx_bed_transcripts["TPM"] = GTEx_bed_transcripts[8].str.split("TPM \"").str[1].str.split("\"").str[0].astype(float)

GTEx_bed_transcripts = GTEx_bed_transcripts[[2,3,4,"ID","TPM"]]
GTEx_bed_transcripts.columns = ["Transcript", "T_START", "T_STOP", "ID", "TPM"]
GTEx_bed_exon = GTEx_bed_exon.drop(columns=[5,7,8,13], axis = 1)
GTEx_bed_exon.columns = ["Chr", "Assembly", "Exon", "E_START", "E_STOP", "STRAND", "Ref_Chr", "Ref_START", "Ref_STOP",
                        "Annotation", "Ref_STRAND", "NT_Overlap", "Ref_length", "Overlap_perc", "E_Coverage", "ID"]
GTEx_both_GDC = pd.merge(how="inner", left=GTEx_bed_exon, right=GTEx_bed_transcripts, left_on="ID", right_on="ID")
GTEx_both_GDC = GTEx_both_GDC[GTEx_both_GDC["TPM"] >= args.TPM_min_novel]
GTEx_both_GDC = GTEx_both_GDC[~((GTEx_both_GDC['Annotation']=='Differential') & (GTEx_both_GDC['E_Coverage'] <= args.Cov_min_diff))]
GTEx_both_GDC = GTEx_both_GDC[~((GTEx_both_GDC['Annotation']=='Differential') & (GTEx_both_GDC['TPM'] <= args.TPM_min_diff))]
GTEx_both_GDC = GTEx_both_GDC[GTEx_both_GDC["E_Coverage"] >= args.Cov_min_novel]
GTEx_both_GDC.to_csv(args.Out_anno, sep="\t", index=False)
GTEx_both_GDC_2 = GTEx_both_GDC.drop_duplicates()
GTEx_both_GDC_2 = GTEx_both_GDC_2.groupby(["E_START", "E_STOP"], group_keys=False).apply(lambda x: x.loc[x.Overlap_perc.idxmax()]).reset_index(drop=True)
Annotation = GTEx_both_GDC_2
print(Annotation)
Annotation_2 = Annotation[(Annotation["Overlap_perc"] < 99) | (Annotation["Annotation"] == "Differential")]
GTF = pd.read_csv(args.GTF, sep='\t', header=None)
GTF["Transcript"] = GTF[8].str.split("transcript_id \"").str[1].str.split("\"").str[0]
GTF["Coverage"] = GTF[8].str.split("cov \"").str[1].str.split("\"").str[0]
Annotation_retained = Annotation[(Annotation["Overlap_perc"] >= 99) & (Annotation["Annotation"] != "Differential")]

mapply.init(
    n_workers=8,
    chunk_size=100,
    max_chunks_per_worker=10,
    progressbar=True
)

def pre_mRNA(Input):
    a = Input["ID"]
    Anno_sub = Annotation[Annotation["ID"] == a]
    GTF_sub = GTF[GTF["Transcript"] == a]
    if len(Anno_sub["E_Coverage"]) == 1:
        GTF_mean = GTF_sub["Coverage"].astype(float).mean()
        for i,j,k in zip(Anno_sub["E_Coverage"],Anno_sub["ID"],Anno_sub["TPM"]):
            if (k > 2):
                return(j)
    else:
        return("X")            

Result = Annotation_retained.mapply(pre_mRNA, axis=1)

try:
    len(Result.columns)
    Result_2 = Result.replace('None', np.nan).bfill(axis=1).iloc[:, 0]
    Final = Result_2.dropna().values.tolist()
    Final = list(set(Final))
    Kept = Annotation_retained[Annotation_retained["ID"].isin(Final)]
    Final_out = pd.concat([Annotation_2, Kept])
except Exception as exception:
    Final = Result.dropna().values.tolist()
    Final = list(set(Final))
    Kept = Annotation_retained[Annotation_retained["ID"].isin(Final)]
    Final_out = pd.concat([Annotation_2, Kept])


Final_out.to_csv(args.Out_anno, sep="\t", index=False)

def get_overlap(Input):
    x = range(Input["E_START"],Input["E_STOP"])
    y = range(Input["Ref_START"],Input["Ref_STOP"])
    xs = set(x)
    Test = xs.intersection(y)
    if len(Test) > 0:
        Input["Overlap_START"] = min(Test)
        Input["Overlap_STOP"] = max(Test)
        Input["ID"] = Input["ID"] + "_" + str(Input["Overlap_START"])
        return(Input)

Overlap = Final_out.apply(get_overlap, axis=1)
Overlap_short = Overlap[["Chr", "Overlap_START", "Overlap_STOP", "ID", "Annotation", "STRAND"]]
#Overlap_short = Overlap_short.drop_duplicates(["Chr", "Overlap_START", "Overlap_STOP", "STRAND"])

Overlap_short.to_csv(args.Out_bed, sep="\t", index=False, header=None)
