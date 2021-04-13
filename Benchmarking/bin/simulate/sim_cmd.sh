#!/usr/bin/env bash

# author: Tzu-Hao Kuo

# If 'source activate' or 'source deactivate' could not work properly, 
# replacing 'source' with 'conda' might be the solution.

source activate py37
export current=`pwd`
export TRSIM_HOME=$( realpath ../../../transmission_simulator/ )
export OUT_DIR=$( realpath ../../data )
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
source activate r35
cd $TRSIM_HOME/R/
Rscript sim_phylo_from_transmission.R -s ../../Benchmarking/config/outbreak_params.yml
source deactivate

###
# simulate the genomes
source activate whichtree_env 
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
source deactivate

###
# simulate the gene presence or absence patterns
source activate r35
cd $TRSIM_HOME/R
Rscript ./sim_gpa_v2.R -t ../../../Benchmarking/data/outbreak_sim/phyloFromPtree.nwk -o $OUT_DIR/gpa.aln
source deactivate
