#!/bin/bash
#$ -V
#$ -l h_core=20
#$ -l h_vmem=300G

source activate cdtree_env
export PATH=/home/thkuo/seq2geno_ng/snps/lib/://net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MM_for_paper/tree_inference/denovo_workflow/lib/Roary/bin/:$PATH
export CDTREE_SHARED_ENV=/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MM_for_paper/tree_inference/shared_envs
#export CDTREE_SHARED_ENV=/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs
mkdir $CONDA_PREFIX/tmp
export TMPDIR=$CONDA_PREFIX/tmp

cd /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MM_for_paper/tree_inference
#./cdtree --config ./test_config.yml --cpu 25 --dry\
./cdtree --cpu 20 --dry\
  /net/flashtest/scratch/thkuo/MM_for_paper/test_klebsiella \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/bin.5/dna_list \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v1/data/reference/GCF_001902475.1_ASM190247v1_genomic.only_chr.fna \
  all

rm -r $CONDA_PREFIX/tmp
