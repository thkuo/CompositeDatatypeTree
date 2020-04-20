#!/usr/bin/env bash

source activate whichtree_env
python ./reads_sim.py \
  -o \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v8/data/intergenic_evol_sim\
  -gd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v8/data/genome_evol_sim/ \
  -wd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v8/data/outbreak_sim/wtrees \
  -n 20 
