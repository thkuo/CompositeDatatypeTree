#!/bin/bash

# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

# activate the environment
source activate cdtree_env
# ensure the paths to those prerequisites (see CompositeDatatypeTree/tree_inference/README.md) 
export PATH=$PATH:/path/to/prokka/:/path/to/roary/:/path/to/stampy/

# make the dna list
# ensure the filename is included in the config file
# please refer to ../../data/pseudomonas_aeruginosa.txt, downloading the reads,
# and then generating a reads list accordingly 
export DNA_LIST=./dna_list

# run the workflow
export RESULTS_D=../results
export REF=$( realpath ../../data/references/NC_008463.fasta )
export LOG_F=cdtree.log
export CONFIG=./cdtree_config.yml
# cdtree -h
cdtree --dry --cpu 2 \
  --config $CONFIG \
  $RESULTS_D $DNA_LIST $REF all 2> $LOG_F

