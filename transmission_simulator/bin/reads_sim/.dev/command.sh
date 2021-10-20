#!/usr/bin/env bash
#$ -V
#$ -l h_core=22
#$ -l h_vmem=100G

source activate whichtree_env
cd /net/metagenomics/data/from_moni/old.tzuhao/transmission_simulator/bin/reads_sim/
parallel -j 3 --joblog /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/reads_sim.log \
'python ./reads_sim.py \
  -o \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/set{}/reads_sim \
  -gd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/set{}/genome_evol_sim \
  -wd \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/transmission_simulator.application/set{}/wtrees \
  -n 7 ' ::: {10..12}
