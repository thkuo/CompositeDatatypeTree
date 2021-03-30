#!/usr/bin/env bash

source activate py37
source activate whichtree_env 
export TRSIM_HOME=./
export OUT_DIR=../../data
# genomes
export GD=$OUT_DIR/genome_evol_sim
# intergenic regions
export ID=$OUT_DIR/intergenic_evol_sim
# within-host subtrees
export WD=$OUT_DIR/outbreak_sim/wtrees
# reads
export RD=$OUT_DIR/reads_sim

#python $TRSIM_HOME/genome_evol_sim/genome_evol_sim_v2.py \
python $TRSIM_HOME/genome_evol_sim/genome_evol_sim.py \
  -o $GD \
  -g ../../config/Streptococcus_pneumoniae_ATCC_700669_v1.db \
  -wd $WD \
  -c ../../config/alf-input_template.drw

python $TRSIM_HOME/intergenic_evol_sim/intergenic_evol_sim.py \
  -o $ID \
  -gd $GD \
  -wd $WD \
  -c ../../config/dawg-input_template.dawg \
  -co ../../config/intergenic_coordinates.txt \
  -n 12

python $TRSIM_HOME/reads_sim/reads_sim.py \
  -o $RD \
  -gd $GD \
  -wd $WD \
  -n 12 
