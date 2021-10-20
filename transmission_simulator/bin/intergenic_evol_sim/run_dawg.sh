#!/usr/bin/env bash

IN_SIM_OUT="$1"
GEN_SIM_OUT="$2"
PROJECT="$3"
DAWG_CONFIG_FILE="$4"
INTERGENIC_COOR_FILE="$5"
#echo $IN_SIM_OUT
#echo $GEN_SIM_OUT
#echo $PROJECT
#echo $DAWG_CONFIG_FILE
CPU_NUM="$6"


#WHICH_TREE_HOME=/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim/bin/re_simulate/which_tree
WHICH_TREE_HOME=$( realpath ../../which_tree )
PATH=$WHICH_TREE_HOME:$PATH

#Using dawg for intergenic regions, with parameter file intergenic.dawg
source activate whichtree_env 
PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2/:$PERL5LIB

mkdir -p $IN_SIM_OUT/$PROJECT
cd $IN_SIM_OUT/$PROJECT
dawg --output $IN_SIM_OUT/$PROJECT/intergenic.fa \
  $DAWG_CONFIG_FILE
$WHICH_TREE_HOME/dawg_to_alf.pl $GEN_SIM_OUT/$PROJECT/speciesMapping.txt \
  $INTERGENIC_COOR_FILE \
  $IN_SIM_OUT/$PROJECT/intergenic.fa
mkdir -p $IN_SIM_OUT/$PROJECT/MSA/
mv MSA_i* $IN_SIM_OUT/$PROJECT/MSA/
#combine with the alf results
#creating the genome sequences
#parallel -j $CPU_NUM --joblog alf_df_to_fasta.log "$WHICH_TREE_HOME/alf_db_to_fasta_splice_intergenic.pl $GEN_SIM_OUT/$PROJECT/DB/{}_dna.fa $IN_SIM_OUT/$PROJECT/MSA | sed 's/>.\+\(DB\/SE[0-9]\+_dna.fa\)/>\1/'> $GEN_SIM_OUT/$PROJECT/DB/{}_genome.fa" \
#  ::: `ls $GEN_SIM_OUT/$PROJECT/DB| grep dna.fa| sed 's/_dna.fa//'`

for n in $( ls $GEN_SIM_OUT/$PROJECT/DB| grep dna.fa| sed 's/_dna.fa//' ); do
  $WHICH_TREE_HOME/alf_db_to_fasta_splice_intergenic.pl \
    $GEN_SIM_OUT/$PROJECT/DB/$n'_dna.fa' \
    $IN_SIM_OUT/$PROJECT/MSA |\
    sed 's/>.\+\(DB\/SE[0-9]\+_dna.fa\)/>\1/'> \
    $GEN_SIM_OUT/$PROJECT/DB/$n'_genome.fa'

done
