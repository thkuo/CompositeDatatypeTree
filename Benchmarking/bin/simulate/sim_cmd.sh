#!/usr/bin/env bash

source activate py37
source activate whichtree_env 
export current=`pwd`
export TRSIM_HOME=../../../transmission_simulator/
export OUT_DIR=../../data
# genomes
export GD=$OUT_DIR/genome_evol_sim
# intergenic regions
export ID=$OUT_DIR/intergenic_evol_sim
# within-host subtrees
export WD=$OUT_DIR/outbreak_sim/wtrees
# reads
export RD=$OUT_DIR/reads_sim

###
# simulate the evolution
cd $TRSIM_HOME/R/
Rscript sim_phylo_from_transmission.R -s ../../config/outbreak_params.yml

###
# simulate the genomes
cd $current
python $TRSIM_HOME/bin/genome_evol_sim/genome_evol_sim.py \
  -o $GD \
  -g ../../config/Streptococcus_pneumoniae_ATCC_700669_v1.db \
  -wd $WD \
  -c ../../config/alf-input_template.drw

python $TRSIM_HOME/bin/intergenic_evol_sim/intergenic_evol_sim.py \
  -o $ID \
  -gd $GD \
  -wd $WD \
  -c ../../config/dawg-input_template.dawg \
  -co ../../config/intergenic_coordinates.txt \
  -n 12

python $TRSIM_HOME/bin/reads_sim/reads_sim.py \
  -o $RD \
  -gd $GD \
  -wd $WD \
  -n 12 

###
# simulate the gene presence or absence patterns
cd $TRSIM_HOME/R
Rscript ./sim_gpa_v2.R
