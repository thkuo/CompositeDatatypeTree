#!/usr/bin/env bash

source activate whichtree_env 
python ./genome_evol_sim.py \
  -o \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v8/data/genome_evol_sim/ \
  -g \
  /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim/bin/re_simulate/which_tree/Streptococcus_pneumoniae_ATCC_700669_v1.db \
  -wd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v8/data/outbreak_sim/wtrees \
  -c \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v8/data/genome_evol_sim/alf-input_template.drw
