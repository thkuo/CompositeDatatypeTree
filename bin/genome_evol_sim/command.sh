#!/usr/bin/env bash

source activate whichtree_env 
parallel -j 3 \
'python ./genome_evol_sim.py \
  -o \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/set{}/genome_evol_sim \
  -g \
  /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim/bin/re_simulate/which_tree/Streptococcus_pneumoniae_ATCC_700669_v1.db \
  -wd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/set{}/wtrees \
  -c \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/data.v5/genome_evol_sim/alf-input_template.drw' ::: {10..12}
