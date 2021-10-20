#!/usr/bin/env bash

READS_OUT_DIR="$1"
GENOME_F="$2"
PROJECT_NAME="$3"
CPU_NUM=${4:1}
WHICH_TREE_HOME=/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim/bin/re_simulate/which_tree
PATH=$WHICH_TREE_HOME:$PATH

source activate whichtree_env

# simulate the reads
mkdir -p $READS_OUT_DIR  
pirs simulate -l 100 -x 100 -m 250 -z \
  -o $READS_OUT_DIR/$PROJECT_NAME $GENOME_F
