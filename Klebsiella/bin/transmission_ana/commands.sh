# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

snakemake -p \
  --use-conda --conda-prefix ./env \
  --snakefile ./transmission_ana.smk -j 15 \
  ../../results/transmission_ana/nuc/clade-wise/1/transphylo/nuc_ctree.Rds \
  ../../results/transmission_ana/nuc/clade-wise/2/transphylo/nuc_ctree.Rds \
  ../../results/transmission_ana/nuc/clade-wise/3/transphylo/nuc_ctree.Rds \
  ../../results/transmission_ana/cd/clade-wise/1/transphylo/cd_ctree.Rds \
  ../../results/transmission_ana/cd/clade-wise/2/transphylo/cd_ctree.Rds \
  ../../results/transmission_ana/cd/clade-wise/3/transphylo/cd_ctree.Rds 

