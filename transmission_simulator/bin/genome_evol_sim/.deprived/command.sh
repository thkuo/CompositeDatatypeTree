#!/bin/bash
#$ -V
#$ -l h_core=18
#$ -l h_vmem=200G

#################
## Preparation ##
#################
#- the tree file
#- the alf parameters
#- the tree, parameters, and coordinates of intergenic regions
CPU_NUM=18
SIM_OUT=/net/metagenomics/data/from_moni/old.tzuhao/transmission_simulator/results/genome_evol_sim
WHICH_TREE_HOME=/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim/bin/re_simulate/which_tree
PATH=$WHICH_TREE_HOME:$PATH
PROJECT=57 

source activate whichtree_env
#cd $SIM_OUT/57

# simulate the evolution
alfsim -o $SIM_OUT ./test.drw

###Using dawg for intergenic regions, with parameter file intergenic.dawg
##PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2/:$PERL5LIB
##dawg --output $SIM_OUT/intergenic_tree_test.fas \
##  $SIM_OUT/intergenic.dawg 
##$WHICH_TREE_HOME/dawg_to_alf.pl $SIM_OUT/$PROJECT/speciesMapping.txt \
##  intergenic_coordinates.txt \
##  $SIM_OUT/intergenic_tree_test.fas
##mkdir -p $SIM_OUT/dawg/MSA/
##mv MSA_i* $SIM_OUT/dawg/MSA/
###combine with the alf results
###creating the genome sequences
##parallel -j $CPU_NUM --joblog alf_df_to_fasta.log "$WHICH_TREE_HOME/alf_db_to_fasta_splice_intergenic.pl $SIM_OUT/$PROJECT/DB/{}_dna.fa $SIM_OUT/dawg/MSA | sed 's/>.\+\(DB\/SE[0-9]\+_dna.fa\)/>\1/'> $SIM_OUT/$PROJECT/DB/{}_genome.fa" \
##  ::: `ls $SIM_OUT/$PROJECT/DB| grep dna.fa| sed 's/_dna.fa//'`

### simulate the reads
##mkdir -p $SIM_OUT/reads/
##parallel -j $CPU_NUM --joblog pirs.log "pirs simulate \
##  -l 100 -x 100 -m 250 \
##  -o $SIM_OUT/reads/{} \
##  $SIM_OUT/$PROJECT/DB/{}_genome.fa " \
##  ::: `ls $SIM_OUT/$PROJECT/DB| grep genome.fa| sed 's/_genome.fa//'`
##parallel -j $CPU_NUM  "gzip -9 {}" \
##  ::: `ls $SIM_OUT/reads/*| grep '.fq$'`
##
### create the annotation files
##source deactivate; source activate prokka_env
##cd $PROJECT 
##mkdir blastn
##parallel -j $CPU_NUM --joblog blastn.log \
##  'blastn -num_threads 5 -query DB/{}_dna.fa -subject DB/{}_genome.fa -evalue 1E-6 -max_hsps 1 -max_target_seqs 2500 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen qframe sframe" | awk "\$4 == \$14"| sort -rk 1 > blastn/{}.tsv' \
##  ::: `cut -f2- speciesMapping.txt`
##source deactivate; source activate py37
##parallel -j $CPU_NUM 'python ../../blastn2gff.py -bl blastn/{}.tsv -genome DB/{}_genome.fa -output blastn/{}.gff' \
##  ::: `cut -f2- speciesMapping.txt`
