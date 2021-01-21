#!/usr/bin/env python3
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import sys
import argparse

def convert_alphabet(aln, alphabet):
    new_aln= aln
    for a in new_aln:
        a.seq.alphabet=alphabet
    return(new_aln)

in_n=snakemake.input["nuc_aln"]
in_g=snakemake.input["gpa_aln"]
out_aln=snakemake.output['raxml_input_aln']
out_prtn=snakemake.output['raxml_input_partitions']
n_aln= AlignIO.read(in_n, 'fasta')
g_aln= AlignIO.read(in_g, 'fasta')
n_aln_l= n_aln.get_alignment_length()
g_aln_l= g_aln.get_alignment_length()
print(n_aln_l)
print(g_aln_l)
# test alignments
try:
    i=set([str(a.id) for a in n_aln]) & set([str(a.id) for a in g_aln])
    u=set([str(a.id) for a in n_aln]) | set([str(a.id) for a in g_aln])
    if len(u-i) > 0:
        raise Exception('The two alignments include different samples')
except Exception as e:
    sys.exit(str(e))
# create the alignment
n_aln.sort()
n_aln=convert_alphabet(n_aln, generic_dna)
g_aln.sort()
g_aln=convert_alphabet(g_aln, generic_dna)
conc_aln= n_aln+g_aln
AlignIO.write(conc_aln, out_aln, 'fasta')
# print the partitions
with open(out_prtn, 'w') as out_fh:
    out_fh.write('DNA, nuc= 1-{}\n'.format(str(n_aln_l)))
    out_fh.write('BIN, gpa= {}-{}\n'.format(
        str(n_aln_l+1), str(n_aln_l+g_aln_l)
    ))
