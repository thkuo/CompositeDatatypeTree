
import argparse
parser= argparse.ArgumentParser()
parser.add_argument('--tr', required= True, dest= 'i', 
    help='input alignment')
parser.add_argument('--out', required= True, dest= 'o', 
    help= 'output table')

args= parser.parse_args()
tr_f= args.i
out_f= args.o
bio_f= ('/net/sgi/metagenomics/data/from_moni/'
        'old.tzuhao/TreePaper/MRKP.v1/data/samples_info/SuppTab1.csv')

## load the tree
from Bio import Phylo
#tr_f= '../../results.5/raxml/RAxML_bestTree.nuc'
tr= Phylo.read(tr_f, 'newick')
tips= [t.name for t in tr.get_terminals()]

## the time table
import os
import pandas as pd
from datetime import datetime
bio_df= pd.read_csv(bio_f, sep= ',', quotechar= '"')
bio_df.set_index('Sequencing Short NO.', inplace= True)

#out_f='../../results.5/transmission_ana/sample_time_subset.tsv'
if not os.path.isdir(os.path.dirname(out_f)):
    os.makedirs(os.path.dirname(out_f))
bio_df.loc[tips, "Specimen date"].to_csv(out_f, sep= '\t', header= False)
