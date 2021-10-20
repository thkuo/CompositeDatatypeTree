#!/usr/bin/env bash
#$ -V
#$ -l h_core=22
#$ -l h_vmem=100G


source activate whichtree_env
cd /net/metagenomics/data/from_moni/old.tzuhao/transmission_simulator/bin/intergenic_evol_sim
parallel -j 3 --joblog /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/intergenic_evol_sim.log \
'python ./intergenic_evol_sim.py \
  -o \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/set{}/intergenic_evol_sim \
  -gd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/set{}/genome_evol_sim \
  -wd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/set{}/wtrees \
  -c \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/data.v5/intergenic_evol_sim/dawg-input_template.dawg \
  -co \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/data.v5/intergenic_evol_sim/intergenic_coordinates.txt \
  -n 7' ::: {10..12}
