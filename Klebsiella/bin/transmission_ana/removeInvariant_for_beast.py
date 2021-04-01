#!/usr/bin/env python3
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
#import string
import argparse
from tqdm import tqdm
import os

def parse_tree_tiplabels(tr_f):
    #' ensure exisiting file
    assert os.path.isfile(tr_f)
    #' consider only the subset of samples or not
    tr= Phylo.read(tr_f, 'newick')
    target_samples= [node.name for node in tr.get_terminals()]
    return(target_samples)

def determine_characters_set(t):
    #' the order should be in line with that used in BEAST
    #' refer to https://github.com/CompEvol/beast2/tree/v2.1.3/src/beast/evolution/datatype
    if t == 'nucleotide':
        return(['a' ,'c', 'g', 't'])
    elif t== 'binary':
        return(['0', '1'])


def main(args):
    aln_f= args.i
    out_f= args.o
    char_type= args.t
    target_samples= [] if args.tr is None else parse_tree_tiplabels(args.tr)

    #' read the alignment
    aln= AlignIO.read(aln_f, 'fasta')
#    print(len(aln))
#    print(aln.get_alignment_length())
    #' filter the samples if specified
    if len(target_samples)>0:
        aln_subset_d= []
        for ix in range(len(aln)): 
            if aln[ix].id in target_samples:
                aln_subset_d.append(aln[ix])
        aln= MultipleSeqAlignment(aln_subset_d)
#    print(len(aln))
#    print(aln.get_alignment_length())
#    print([rec.id for rec in aln ])
    aln.sort()
    sample_num= len(aln)
    col_num= aln.get_alignment_length()

    new_aln= aln[:, 0:1] ## removed afterward
    print([s.id for s in new_aln])
    #' The order is important for BEAST to parse
    state_chars_set=determine_characters_set(char_type)
    inv_counts= {s:0 for s in state_chars_set}
    for n in tqdm(range(col_num)):
        col_str= aln[:, n]
        col_str= col_str.lower()
        #' determine all unique states
        u_states=set(col_str) 
        if (len(u_states) == 1) and (all([s in state_chars_set for s in u_states])):
            state= list(u_states)[0]
            inv_counts[state]=inv_counts[state]+1
        else:
            col_aln= aln[:, n:n+1]
            new_aln=new_aln+col_aln
    new_aln= new_aln[:, 1:]

    #' ensure the sum
    assert (sum(inv_counts.values())+new_aln.get_alignment_length() 
            == aln.get_alignment_length())

    with open(out_f, 'w') as out_fh:
        AlignIO.write(new_aln, out_fh, "fasta")
    out_nex= out_f+'.nex'
    AlignIO.convert(out_f, 'fasta', out_nex, 'nexus', generic_dna)
    weights_out_f= out_f+'.inv_weights'
    with open(weights_out_f, 'w') as weights_out_fh:
        weights_out_fh.write('#'+
        ' '.join([s for s in state_chars_set])+'\n')
        weights_out_fh.write(
        ' '.join([str(inv_counts[s]) for s in state_chars_set]))

if __name__== '__main__':
    # read an alignment
    parser= argparse.ArgumentParser()
    parser.add_argument('--in', required= True, dest= 'i', 
        help='input alignment')
    parser.add_argument('--out', required= True, dest= 'o', 
        help= 'output alignment')
    parser.add_argument('--type', dest= 't', 
        choices=['nucleotide','binary'],default= 'nucleotide',
        help= 'character types (default: %(default))')
    parser.add_argument('--tr', required= False, dest= 'tr', 
        default= None, 
        help='''to consider only a subset of the samples composed of the
                        alignment, the tip labels of this tree is used''')

    args= parser.parse_args()
    main(args)


## detect invariant sites
#detect_result=[]
#for name in aln:
#    for n in range(len(str(aln[name].seq))):
#        char= (str(aln[name].seq[n]).upper() if (not case) else str(aln[name].seq[n]))
#        if len(detect_result) < len(str(aln[name].seq)):
#            detect_result.append({char:1})
#        elif not (char in detect_result[n]):
#            detect_result[n][char]=1
#        else:
#            detect_result[n][char]+=1
#
#variant_counts= [len(aln)-compo[sorted(compo, key= lambda x:-compo[x])[0]] for compo in detect_result]
#variant_loc= [n for n in range(len(variant_counts))if int(variant_counts[n]) >= cutoff]
#
#
#new_recs= []
#for name in aln:
#    new_seq= ''.join([str(aln[name].seq[n]) for n in variant_loc])
#    new_rec= SeqRecord(Seq(new_seq,IUPAC.IUPACAmbiguousDNA), id=name, name= '', description= '')
#    new_recs.append(new_rec)
#
#with open(out_f, 'w') as out_h:
#    SeqIO.write(new_recs, out_h, 'fasta')
