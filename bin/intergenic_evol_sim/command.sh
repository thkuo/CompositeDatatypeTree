#!/usr/bin/env bash

source activate whichtree_env
python ./intergenic_evol_sim.py \
  -o \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/data/intergenic_evol_sim\
  -gd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/data/genome_evol_sim/ \
  -wd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/data/outbreak_sim/wtrees \
  -c \
  ./dawg-input_template.dawg \
  -co \
  ./intergenic_coordinates.txt
  -n 20 
