#!/usr/bin/env bash
#$ -V
#$ -l h_core=30
#$ -l h_vmem=300G

cd /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/bin/conc_workflow
source activate sgp_env

snakemake -p \
 --rerun-incomplete \
 --conda-prefix /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs \
 --use-conda -j 28 \
 --configfile=conc.config.yml \
 --snakefile=conc_tree.smk
