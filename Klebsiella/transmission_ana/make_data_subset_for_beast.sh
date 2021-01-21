#!/usr/bin/env bash

## the time data
#conda activate cdtree_env
#parallel 'python ./prune_time_data.py \
#  --tr ../../results.5/transmission_ana/{1}/clade-wise/{2}/subtree.nwk \
#  --out ../../results.5/transmission_ana/{1}/clade-wise/{2}/time.tab' ::: nuc.15 cd.15 ::: 1 3 4
#conda deactivate

# the alignment subset
conda activate cdtree_env
# nuc
export nuc_aln=../../results.5/alignment/nuc.full.aln
parallel -j 6 'python removeInvariant_for_beast.py \
  --type nucleotide \
  --in $nuc_aln \
  --out ../../results.5/transmission_ana/{1}/clade-wise/{2}/nuc.sub_compressed.aln \
  --tr ../../results.5/transmission_ana/{1}/clade-wise/{2}/subtree.nwk '  ::: nuc.15 cd.15 ::: 1 3 4
## gene presence/absence
#export gpa_aln=../../results.5/gpa/gpa.aln
#parallel 'python removeInvariant_for_beast.py \
#  --type binary \
#  --in $gpa_aln \
#  --out ../../results.5/transmission_ana/{1}/clade-wise/{2}/gpa.sub_compressed.aln \
#  --tr ../../results.5/transmission_ana/{1}/clade-wise/{2}/subtree.nwk '  ::: cd.15 ::: 1 3 4
conda deactivate
