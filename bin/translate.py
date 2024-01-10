#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import mapply
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import re as re
import difflib
import itertools
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--Short_GTF', type=str)
parser.add_argument('--Stringtie_tsv', type=str)
parser.add_argument('--Gencode_GTF', type=str)
parser.add_argument('--Reference_tsv', type=str)
parser.add_argument('--Regions_bed', type=str)
parser.add_argument('--peptide_length', nargs="+", default=[9, 15], type=int)
parser.add_argument('--closest_out', type=str)
parser.add_argument('--matching_out', type=str)
parser.add_argument('--Bed_out', type=str)
parser.add_argument('--Overlaps_out', type=str)
parser.add_argument('--Peptides_out', type=str)
parser.add_argument('--Regions_fasta', type=str)
parser.add_argument('--Ref_range', type=int)
args = parser.parse_args()

GTF_07H103 = pd.read_csv(args.Short_GTF, sep='\t', header=None)
GTF_07H103 = GTF_07H103[GTF_07H103[2] == "transcript"]
FASTA_07H103 = pd.read_csv(args.Stringtie_tsv, sep='\t', header=None)
FASTA_07H103.columns = ["ID", "Sequence", "None"]
FASTA_07H103 = FASTA_07H103.drop("None", axis = 1)
GTF_07H103["ID"] = GTF_07H103[8].str.split("transcript_id \"").str[1].str.split("\"").str[0]
GTF_07H103 = pd.merge(how="inner", left=GTF_07H103, right=FASTA_07H103, left_on="ID", right_on="ID")
GTF_gencode = pd.read_csv(args.Gencode_GTF, sep='\t', header=None, skiprows= 7)
GTF_gencode["ID"] = GTF_gencode[8].str.extract(r'transcript_id "([^"]+)"')
GTF_gencode_names = pd.read_csv(args.Reference_tsv, sep='\t', header=None)
GTF_gencode_names.columns = ["ID", "Protein", "None"]
GTF_gencode_names = GTF_gencode_names.drop("None", axis = 1)
print(GTF_gencode_names["ID"])
GTF_gencode_names["ID"] = GTF_gencode_names["ID"].str.split("_").str[0]
print(GTF_gencode_names["ID"])
GTF_gencode_names["ID"] = GTF_gencode_names["ID"].str.split("|").str[1]
GTF_gencode = GTF_gencode[GTF_gencode[2] == "transcript"]
GTF_gencode = pd.merge(how="inner", left=GTF_gencode, right=GTF_gencode_names, left_on="ID", right_on="ID")
GTF_gencode = GTF_gencode[GTF_gencode[0].str.contains("chr")]
GTF_gencode = GTF_gencode.fillna("NONE")
GTF_gencode = GTF_gencode[GTF_gencode["Protein"] != "NONE"]
GTF_07H103 = GTF_07H103[GTF_07H103[0].str.contains("chr")]
GTF_07H103 = GTF_07H103[GTF_07H103[0].str.contains("_") == False]
GTF_07H103 = GTF_07H103[GTF_07H103[6].str.contains("\.") == False]

mapply.init(
    n_workers=16,
    chunk_size=100,
    max_chunks_per_worker=10,
    progressbar=True
)

def get_overlap(Merged):
    s1 = Merged[0]
    s2 = Merged[1].split("_")[0]
    s = difflib.SequenceMatcher(None, s1, s2, autojunk = False)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
    Final = s1[pos_a:]
    
    try:
        i = Final.index('4')
    except ValueError as e:
        i = len(Final)
        
    j = s1.rfind('4', 0, pos_a)    
    if j != -1:
        Final = s1[j+1:pos_a]+Final[:i]
    else:
        Final = s1[:pos_a]+Final[:i]
    return s1, s2, len(s1[pos_a:pos_a+size]), Final, j, i, pos_a

def get_matching(Input):
    rna_seqs = [Seq(Input.Sequence)]
    gene_id = Input.Transcript
    gene_id = Input.Transcript.split('STRG')[1]
    gene_id = "STRG" + gene_id
    aa_seqs = (s[i:].translate(to_stop=False, cds=False, stop_symbol='4') for i in range(3) for s in rna_seqs)
    aa_seqs2 = ','.join(str(v) for v in aa_seqs)
    aa_seqs2 = aa_seqs2.split(',')
    aa_seqs2 = list(filter(None, aa_seqs2))
    if Input.Translation_ref != "No_ref":
        Merged = list(itertools.product(aa_seqs2, [Input.Ref_protein]))
        Test_4 = pd.DataFrame(map(get_overlap, Merged))
        Test_max = Test_4[Test_4[2]==Test_4[2].max()]
        Input["Translated"] = str(Test_max[3].iloc[0])
        Input["Shift"] = Test_max[3].index.values[0]
        Input["Overlap"] = Test_max[2].iloc[0]
        Input = Input.values.tolist()
    else:
        try:
            Kozak_1 = str(rna_seqs[0]).index('GCCACCATGG')
            Kozak_1 = int(round(Kozak_1/3))
            print(Kozak_1, "kozak")
        except ValueError as e:
            Kozak_1 = 10000
            
        try:
            Kozak_2 = str(rna_seqs[0]).index('GCCGCCATGG')
            Kozak_2 = int(round(Kozak_2/3))
            print(Kozak_2, "kozak")
        except ValueError as e:
            Kozak_2 = 10000
            
        x = -1
        Positions = []
        for s in aa_seqs2:
            x = x +1
            try:
                i = s.index('M')
                a = s[i:]
                all_positions = [pos for pos, char in enumerate(s) if char == 'M']
                arr = np.asarray(all_positions)
                Kozak_3 = (np.abs(arr - Kozak_1)).argmin()
                Kozak_4 = (np.abs(arr - Kozak_2)).argmin()
                Koz_1 = Kozak_1 - arr[Kozak_3]
                Koz_2 = Kozak_2 - arr[Kozak_4]
                if abs(Koz_1) < 50:
                     b = s[arr[Kozak_3]:]  
                else:
                     b = "X"
                if abs(Koz_2) < 50:
                     b = s[arr[Kozak_4]:]         
            except ValueError as e:
                i = 10000
                a = s
                b = "X"
            Positions.append([x, i, a, b, len(b)])
        Positions = pd.DataFrame(Positions)
        Final = Positions[Positions[4]==Positions[4].max()]
        Final_str = Final[3].iloc[0]
        if len(Final_str) < 8:
            Final = Positions[Positions[1]==Positions[1].min()]
            Final_str = Final[2].iloc[0]
        else:
            print("MALAKA")
        try:
            j = Final_str.index('4')
        except ValueError as e:
            j = len(Final_str)
        
        Final_1 = Final_str[:j]
        Input["Translated"] = Final_1
        Input["Shift"] = Final[0].iloc[0]
        Input["Overlap"] = 0
        Input = Input.values.tolist()
    return(Input)

def closest_value(Input):
  input_value = Input[3]
  Name = Input["ID"]
  GTF_strand = GTF_gencode[(GTF_gencode[6] == Input[6]) & (GTF_gencode[0] == Input[0])]
  input_list = GTF_strand[3].to_list()  
  arr = np.asarray(input_list)
  i = (np.abs(arr - input_value)).argmin()
  Test_1 = GTF_strand[GTF_strand[3] == arr[i]]
  Test_2 = Test_1["ID"].iloc[0]
  Diff = Input[3] - abs(arr[i])
  if (abs(Diff) <= args.Ref_range) & (Test_1["Protein"].iloc[0] != "NONE"):
            Confidence = Test_2
  elif (abs(Diff) > args.Ref_range) & (Test_1["Protein"].iloc[0] != "NONE"):
            Confidence = "No_ref"        
  return [Input[3], Input["Sequence"], Test_1[3].iloc[0], Test_1[6].iloc[0], Test_1["Protein"].iloc[0], Diff, Name, Test_2, Confidence]

seq_dict = {}
for rna_record in SeqIO.parse(args.Regions_fasta, 'fasta'):
        # use both fwd and rev sequences
        rna_seqs = [rna_record.seq]

        # generate all translation frames
        aa_seqs = ([s[i:].translate(to_stop=False,cds=False), i] for i in range(3) for s in rna_seqs)
        aa_seqs_test = ([s[i:].translate(to_stop=False,cds=False), i] for i in range(3) for s in rna_seqs)
        aa_seqs2 = ','.join(str(v) for v,u in aa_seqs)
        aa_seqs2 = aa_seqs2.split(',')
        for seq in aa_seqs2:
            if len(seq) > 7:
                gene_id = rna_record.id.split("STRG")[1].split("_")[0]
                gene_id = "STRG" + gene_id
                if not seq_dict.get(gene_id): 
                    seq_dict[gene_id] = [] 
                for v,u in aa_seqs_test:
                    if seq in v:
                        a = u
                        break    
                seq_dict[gene_id].append([rna_record.id, seq, a])


Result = GTF_07H103.mapply(closest_value, axis = 1)

try:
    len(Result.columns)
    Test = Result.replace('None', np.nan).bfill(axis=1).iloc[:, 0]
    df = pd.DataFrame.from_records(Test.reset_index(drop=True))
    df.columns = ["my_START", "Sequence", "Ref_START", "STRAND", "Ref_protein", "Diff", "Transcript", "ref_Transcript", "Translation_ref"]
    df = df.sort_values(by=['Translation_ref'])
    df = df.drop_duplicates(["Ref_protein", "Transcript"], keep="first")
except Exception as exception:
    df = pd.DataFrame.from_records(Result.reset_index(drop=True))
    df.columns = ["my_START", "Sequence", "Ref_START", "STRAND", "Ref_protein", "Diff", "Transcript", "ref_Transcript", "Translation_ref"]
    df = df.sort_values(by=['Translation_ref'])
    df = df.drop_duplicates(["Ref_protein", "Transcript"], keep="first")

df.to_csv(args.closest_out, sep='\t', index = False)

Test_5 = df.mapply(get_matching, axis=1)

try:
    len(Test_5)
    Test = Test_5.replace('None', np.nan).bfill(axis=1).iloc[:, 0]
    df = pd.DataFrame.from_records(Test.reset_index(drop=True))
except Exception as exception:
    df = pd.DataFrame.from_records(Test_5.reset_index(drop=True))

df.to_csv(args.matching_out, sep='\t', index = False)

records2 = df[[9, 6]]
records2.columns = ["Translated", "Transcript"]


BED = pd.read_csv(args.Regions_bed, sep = "\t", engine='python', header = None)
BED.columns = ["Chr", "Overlap_START", "Overlap_STOP", "ID", "Annotation", "STRAND"]
BED["ID_2"] = BED["ID"].str.rsplit('_', n=1).str[0]
Test_8 = df[[6, 8]]
BED_3 = pd.merge(how="inner", left=BED, right=Test_8, left_on="ID_2", right_on=6)
BED_3 = BED_3.drop(["ID_2", 6], axis = 1)
BED_3.rename(columns={8:'Translation_ref'}, inplace=True)

def get_overlap_2(s1, s2, edge_2):
    edge = edge_2 -1
    s = difflib.SequenceMatcher(None, s1, s2, autojunk = False)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
    print(pos_a)
    print(pos_b)
    if pos_b -edge > 0:
        A = -edge
    else:
        A = 0 - pos_b
    if pos_b + size + edge < len(s2):
        B = edge
    else:
        B = len(s2) - (pos_b + size)
    if len(s2[pos_b:pos_b+size]) > edge:    
        return s2[pos_b+A:pos_b+size+B], s2[pos_b:pos_b+size], A, B
    else:
        return s2[pos_b:pos_b+size], s2[pos_b:pos_b+size], A, B

# START and STOP of region within translated full transcript based on genomic location 
def overlap_position_strand(Input):
    if Input["STRAND"] == "+":
        Input["Overlap_START_X"] = Input["Full"].index(Input["Original"]) + Input["minus_shift"]
        Input["Overlap_STOP_X"] = len(Input["Overlap"])# + Input["plus_shift"]
        Input["Overlap_START_1"] = Input["Overlap_START"] + (Input["Overlap_START_X"]*3)
        Input["Overlap_STOP_1"] = Input["Overlap_START_1"] + (Input["Overlap_STOP_X"]*3)
    else:
        Input["Overlap_START_X"] = Input["Full"].index(Input["Original"]) + Input["minus_shift"]
        Input["Overlap_STOP_X"] = len(Input["Overlap"])# + Input["plus_shift"]
        Input["Overlap_STOP_1"] = Input["Overlap_STOP"] - (Input["Overlap_START_X"]*3)
        Input["Overlap_START_1"] = Input["Overlap_STOP_1"] - (Input["Overlap_STOP_X"]*3)
    return Input

def get_sliding(Input, length=1):
    Splitted = []
    if Input["STRAND"] == "+":
        peptides = range(0,len(Input["Overlap"]))
        regions = range(Input["Overlap_START_1"],Input["Overlap_STOP_1"],3)
        for (p,j) in zip(peptides, regions):
            x = j + Input["Shift"]
            y = j + length*3 + Input["Shift"]
            if (p+length <= len(Input["Overlap"])) & (y <= Input["Overlap_STOP_1"]):
                p2 = Input["Overlap"][p:p+length]
                Splitted.append([Input["Chr"], x, y, Input["ID"] + "_" + str(x) +  "_" + str(y), p2, Input["Annotation"], Input["STRAND"], length])            
    if Input["STRAND"] == "-":
        peptides = range(0,len(Input["Overlap"]))
        regions = range(Input["Overlap_STOP_1"],Input["Overlap_START_1"],-3)
        for (p,j) in zip(peptides, regions):
            x = j - Input["Shift"]
            y = j - length*3 - Input["Shift"]
            if (p+length <= len(Input["Overlap"])):
                p2 = Input["Overlap"][p:p+length]
                Splitted.append([Input["Chr"], y, x, Input["ID"] + "_" + str(x) +  "_" + str(y), p2, Input["Annotation"], Input["STRAND"], length])
    return(Splitted)


# Chop to peptides, sequence and genomic location
def all_comes_together(Input, params=[1,2,3,4]):
    appended_data = []
    for i in params:
        seq_2 = []
        for z,j in Input.itertuples(index=False):
            gene_id = j.split('STRG')[1].split("_")[0]
            gene_id = "STRG" + gene_id
            if gene_id in seq_dict:
                for seq in seq_dict[gene_id]:
                    overlap = len(get_overlap_2(str(z).upper(), str(seq[1]).upper(), i)[0])
                    overlap_2 = get_overlap_2(str(z).upper(), str(seq[1]).upper(), i)[0]
                    original = get_overlap_2(str(z).upper(), str(seq[1]).upper(), i)[1]
                    minus_shift = get_overlap_2(str(z).upper(), str(seq[1]).upper(), i)[2]
                    plus_shift = get_overlap_2(str(z).upper(), str(seq[1]).upper(), i)[3]
                    seq_2.append([overlap, str(seq[0]), overlap_2.upper(), str(z).upper(), str(seq[1]).upper(), seq[2], original.upper(), minus_shift, plus_shift])

        Test = pd.DataFrame({'Count': [el[0] for el in seq_2 ], 'Sample': [el[1] for el in seq_2 ], 
                            'Overlap': [el[2] for el in seq_2 ], 'Ref': [el[3] for el in seq_2 ], 'Full': [el[4] for el in seq_2 ], 'Shift': [el[5] for el in seq_2 ],
                            'Original': [el[6] for el in seq_2 ], 'minus_shift': [el[7] for el in seq_2 ], 'plus_shift': [el[8] for el in seq_2 ]})
        print(Test)
        idx = Test.groupby(['Sample'])['Count'].transform(max) == Test['Count']
        Test = Test[idx]
        Test = Test.drop_duplicates(["Sample", "Count"])

        # Match must be at least 7 aas long
        Test = Test[Test["Count"] > i -1]

        #Test = Test.apply(overlap_position, axis=1)

        # Merge with Metadata
        Combined = pd.merge(left=Test, right=BED_3, left_on="Sample", right_on="ID")

        Combined = Combined.apply(overlap_position_strand, 1)

        Combined["Sample"] = Combined["Sample"].str.rsplit('_', n=1).str[0]
        Combined["ID"] = Combined["ID"].str.rsplit('_', n=1).str[0]

        Combined.to_csv(args.Overlaps_out, sep='\t', index = False)

        Test_3 = Combined.apply(get_sliding, length=i, axis=1)
        my_list = [j for i in Test_3 for j in i ]
        Test_4 = pd.DataFrame(my_list)
        Test_4 = Test_4.drop_duplicates()
        appended_data.append(Test_4)
    appended_data = pd.concat(appended_data)
    return appended_data

#Test_4 = records2.apply(all_comes_together, params=(args.peptide_length), axis=1)
Test_4 = all_comes_together(records2, params=(args.peptide_length))

#try:
#    len(Test_4)
#    Test = Test_4.replace('None', np.nan).bfill(axis=1).iloc[:, 0]
#    df = pd.DataFrame.from_records(Test.reset_index(drop=True))
#except Exception as exception:
#    df = pd.DataFrame.from_records(Test_4.reset_index(drop=True))


Test_export = Test_4[[3, 4]]

with open(args.Peptides_out, 'w') as aa_fa:
    for j,i in Test_export.itertuples(index=False):
        if i != "None":
            if len(i) > 12:
                aa_record = SeqRecord(Seq(i), id=j, description="MHCII")
                SeqIO.write(aa_record, aa_fa, 'fasta')
            else:
                aa_record = SeqRecord(Seq(i), id=j, description="MHCI")
                SeqIO.write(aa_record, aa_fa, 'fasta')    
Test_4.to_csv(args.Bed_out, sep='\t', header= None, index = False)