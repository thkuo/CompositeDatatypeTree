#!/usr/bin/env bash
#$ -V
#$ -l h_core=24
#$ -l h_vmem=300G

cd /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/bin.5/
source activate cdtree_env
# reuse the same method
export method_dir=/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/bin.4

# snps calling
export PATH=/home/thkuo/seq2geno2pheno/seq2geno/snps/lib:$PATH

## the nucleotide tree
#bash ./nuc_workflow/cmd.sh

snakemake -p \
    --rerun-incomplete --use-conda \
    --conda-prefix /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs \
    --cores 22 \
    --directory=per-site_LogLikelihood/shortcut/ \
    --configfile=per-site_LogLikelihood/shortcut/config.yml \
    --snakefile=per-site_LogLikelihood/shortcut/compute_ll.smk \
    /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/results.5/nuc_col_by_ll.nwk
