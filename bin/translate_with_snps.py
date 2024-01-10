#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import re as re
import difflib

parser = argparse.ArgumentParser()
parser.add_argument('--fasta_in', type=str)
parser.add_argument('--peptides', type=str)
parser.add_argument('--regions', type=str)
parser.add_argument('--tsv_out', type=str)
parser.add_argument('--tsv_out_meta', type=str)
args = parser.parse_args()

def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2)) 
    return s1[pos_a:pos_a+size]

def get_mismatch(s1, s2):
    x = [(s1[i],s2[i],i+1) for i in range(len(s1)) if s1[i] != s2[i]]
    return x

seq_dict = {}
for rna_record in SeqIO.parse(args.fasta_in, 'fasta'):
        # use both fwd and rev sequences
        rna_seqs = [rna_record.seq]

        # generate all translation frames
        aa_seqs = (s[i:].translate(to_stop=False,cds=False) for i in range(3) for s in rna_seqs)
        aa_seqs2 = ','.join(str(v) for v in aa_seqs)
        aa_seqs2 = aa_seqs2.split(',')
        for seq in aa_seqs2:
                gene_id = rna_record.id.rsplit('_', 1)[0]
                if not seq_dict.get(gene_id): 
                    seq_dict[gene_id] = [] 
                seq_dict[gene_id].append([rna_record.id, seq, str(rna_record.seq)])

records2 = list(SeqIO.parse(args.peptides, "fasta"))

NT_ref = list(SeqIO.parse(args.regions, "fasta"))

seq_2 = []
for rec,nt in zip(records2,NT_ref):
   gene_id = rec.id.split('STRG')[1]
   gene_id = "STRG" + gene_id
   if gene_id in seq_dict:
        for seq in seq_dict[gene_id]:
            overlap = len(get_overlap(str(rec.seq).upper(), str(seq[1]).upper()))
            overlap_2 = get_overlap(str(rec.seq).upper(), str(seq[1]).upper())
            mismatch = get_mismatch(str(seq[2]).upper(), str(nt.seq).upper())
            seq_2.append([overlap, str(seq[0]), overlap_2.upper(), str(seq[1]).upper(), str(rec.seq).upper(), str(nt.seq).upper(), str(seq[2]).upper(), mismatch])

Test = pd.DataFrame({'Count': [el[0] for el in seq_2 ], 'Sample': [el[1] for el in seq_2 ],
                     'Overlap': [el[2] for el in seq_2 ], 'Full': [el[3] for el in seq_2 ], 
                     'REF': [el[4] for el in seq_2 ], 'REF_NT': [el[5] for el in seq_2 ],
                     'NEW_NT': [el[6] for el in seq_2 ], 'Mismatch': [el[7] for el in seq_2 ]})
idx = Test.groupby(['Sample'])['Count'].transform(max) == Test['Count']
Test = Test[idx]
Test = Test.drop_duplicates(["Sample", "Count"])
Test[~Test["Full"].str.contains("\*")]
Test_export = Test[["Sample", "Full"]]
with open(args.tsv_out, 'w') as aa_fa:
    for j,i in Test_export.itertuples(index=False):
        if (i != "None") & (len(i) > 7) & ("*" not in i):
            aa_record = SeqRecord(Seq(i), id=j, description="translated sequence")
            SeqIO.write(aa_record, aa_fa, 'fasta')

Test.to_csv(args.tsv_out_meta, sep='\t', index = False)