#!/usr/bin/env bash
#$ -V
#$ -l h_core=16
#$ -l h_vmem=200G

cd /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5.7/ll_cutoff
source activate cdtree_env

#export CONFIG_F=./conc.config.yml
export CONFIG_F=./conc.config.old_gpa.yml
snakemake -p \
 --rerun-incomplete \
 --conda-prefix /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MM_for_paper/tree_inference/shared_envs \
 --use-conda -j 16 \
 --configfile=$CONFIG_F \
 --snakefile=./conc_tree.smk 
source deactivate
