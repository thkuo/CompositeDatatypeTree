'''
Re-format the alignment files for BEAST2

'''
from Bio import AlignIO
from Bio.Alphabet import generic_dna

nuc_aln_fa='../../results.5/alignment/nuc.var.aln'
nuc_aln_nex='../../results.5/transmission_ana/nuc.var.nex'
AlignIO.convert(nuc_aln_fa, 'fasta', nuc_aln_nex, 'nexus', generic_dna)
gpa_aln_fa='../../results.5/gpa/gpa.var.aln'
gpa_aln_nex='../../results.5/transmission_ana/gpa.var.nex'
AlignIO.convert(gpa_aln_fa, 'fasta', gpa_aln_nex, 'nexus', generic_dna)
# there seems no alphabet for binary states, so the gpa file needs manual
# edits for its type
