#!/usr/bin/env python3
from Bio import SeqIO
from Bio import AlignIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
#import string
import argparse

def main(args):
    aln_f= args.i
    out_f= args.o
    cutoff= int(args.cutoff_num)
    case= args.case_sensitive
    aln= AlignIO.read(aln_f, 'fasta')
    sample_num= len(aln)
    col_num= aln.get_alignment_length()

    new_aln= aln[:, 0:1] ## removed afterward
    new_aln.sort()
    print([s.id for s in new_aln])
    for n in range(col_num):
        col_str= aln[:, n]
        if not case:
            col_str= col_str.lower()
        dominant_char_num= max(
            [col_str.count(uc) for uc in set(col_str)])
        if (sample_num-dominant_char_num) >= cutoff:
            col_aln= aln[:, n:n+1]
            col_aln.sort()
            new_aln=new_aln+col_aln
    new_aln= new_aln[:, 1:]

    with open(out_f, 'w') as out_fh:
        AlignIO.write(new_aln, out_fh, "fasta")

if __name__== '__main__':
    # read an alignment
    parser= argparse.ArgumentParser()
    parser.add_argument('--in', required= True, dest= 'i', 
        help='input alignment')
    parser.add_argument('--out', required= True, dest= 'o', 
        help= 'output alignment')
    parser.add_argument('--cn', dest= 'cutoff_num', default= 2,
        help= 'least number of variant residue should exist in column')
    parser.add_argument('--s', dest= 'case_sensitive', default= False, 
        action= 'store_true', 
        help= 'case-sensitive (ie. A and a should be viewed differently)')

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
