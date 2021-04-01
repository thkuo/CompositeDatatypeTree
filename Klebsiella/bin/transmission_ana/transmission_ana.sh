#!/bin/bash
#$ -V
#$ -l h_core=15
#$ -l h_vmem=300G

###
# adjust the memory allowance in /home/thkuo/miniconda3/envs/beast25_env/bin/beast
#
cd /net/personal-backup/thkuo/MRKP.v2/bin.5/transmission_ana
#source activate beast25_env
#parallel -j 6 --joblog beast.log \
#  'beast -working -seed 1234 -threads 3 /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results.5/transmission_ana/{1}/clade-wise/{2}/v2/beast.xml' \
#  ::: nuc.15 cd.15 ::: 3 4
#source deactivate

# tree dating
#source activate snakemake_unittest_env
#snakemake -p \
#  --use-conda --conda-prefix /net/personal-backup/thkuo/MRKP.v2/bin.5/transmission_ana/env \
#  --snakefile ./transmission_ana.smk -j 15 \
#  /net/personal-backup/thkuo/MRKP.v2/results.5/transmission_ana/cd.15/clade-wise/1/v7/cd.d_tree.nex \
#  /net/personal-backup/thkuo/MRKP.v2/results.5/transmission_ana/cd.15/clade-wise/3/v7/cd.d_tree.nex \
#  /net/personal-backup/thkuo/MRKP.v2/results.5/transmission_ana/cd.15/clade-wise/4/v7/cd.d_tree.nex
#source deactivate

# transmission 
#export data_d=/net/personal-backup/thkuo/MRKP.v2/results.5/transmission_ana/cd.15/clade-wise/
#source activate r35
#parallel -j 3 'Rscript ./transmission_ana.R -t $data_d/{}/v7/cd.d_tree.nex -d $data_d/{}/time.tab -p cd -o $data_d/{}/v7/transphylo' ::: 1 3 4
#for d in $( ls $data_d ); do
#  echo $data_d/$d
#  export d_tree=$data_d/$d/v7/cd.d_tree.nex
#  export time_tab=$data_d/$d/time.tab 
#  export tree_name='cd'
#  export out_d=$data_d/$d/transphylo
#  if [ ! -s $d_tree ];then
#    echo $d_tree' not found'
#    exit
#  fi
#  echo $d_tree
#  Rscript ./transmission_ana.R \
#  -t $d_tree \
#  -d $time_tab \
#  -p $tree_name -o $out_d 
#done
#source deactivate

source activate snakemake_unittest_env
snakemake -p \
  --use-conda --conda-prefix /net/personal-backup/thkuo/MRKP.v2/bin.5/transmission_ana/env \
  --snakefile ./transmission_ana.smk -j 15 \
  /net/personal-backup/thkuo/MRKP.v2/results.5/transmission_ana/cd.15/clade-wise/3/v_neg_control/cd.d_tree.nex \
  /net/personal-backup/thkuo/MRKP.v2/results.5/transmission_ana/nuc.15/clade-wise/3/v_pos_control/nuc.d_tree.nex
source deactivate
