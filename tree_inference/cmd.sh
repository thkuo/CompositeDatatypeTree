#!/bin/bash
#$ -V
#$ -l h_core=26
#$ -l h_vmem=200G

source activate cdtree_env
export CDTREE_SHARED_ENV=/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs
mkdir $CONDA_PREFIX/tmp
export TMPDIR=$CONDA_PREFIX/tmp

cd /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MM_for_paper/tree_inference
#./cdtree --config ./test_config.yml --cpu 25 --dry\
./cdtree --cpu 25 --dry\
  ./test \
  ../../WhichTree_Sim.v7/bin.v5/run_seq2geno/dna_list \
  /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim/data/reference/ATCC_700669.fasta \
  fast_mapping

rm -r $CONDA_PREFIX/tmp
