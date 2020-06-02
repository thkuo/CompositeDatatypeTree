#!/usr/bin/env bash
#$ -V
#$ -l h_core=30
#$ -l h_vmem=300G

cd /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/bin
source activate sgp_env

## seq2geno-ng
# snps calling
export PATH=/home/thkuo/seq2geno2pheno/seq2geno/snps/lib:$PATH
#
## reuse the precomputed variants 
#snakemake -p --cores 16 --use-conda\
#  --rerun-incomplete \
#  --configfile=snps_config.yml \
#  --directory=snps_workflow/ \
#  --snakefile=snps_workflow/snps.in_one.smk \
#  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.vPre3-3/results/seq2geno/snps/merged_vcf/multisample.snp.vcf.gz

## the nucleotide tree
#snakemake -p --cores 20 --use-conda \
#  --rerun-incomplete \
#  --configfile=nuc_workflow/nuc_config.yml \
#  --conda-prefix=env/ \
#  --directory=nuc_workflow/ \
#  --snakefile=nuc_workflow/nuc_tr.smk \
#  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/results/raxml/RAxML_bipartitions.nuc.bs

## make the gpa alignment


## composite datatype tree
