# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

source activate cdtree_env
export CDTREE_SHARED_ENV=/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MM_for_paper/tree_inference/shared_envs
mkdir $CONDA_PREFIX/tmp
export TMPDIR=$CONDA_PREFIX/tmp

output_dir=...
dna_list=...
reference_seq=...
workflow=all
./cdtree --cpu 20  \
  --config ./example_config.yml \
  $test_dir \
  $dna_list \
  $reference_seq \
  $workflow

rm -r $CONDA_PREFIX/tmp
