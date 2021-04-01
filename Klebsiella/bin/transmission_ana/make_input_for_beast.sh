
##' count the constant sites for the nuc alignment
#export NUC_ALN=/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results.5/alignment/nuc.full.aln
#export OUT_NUC_ALN=/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results.5/alignment/nuc.compressed.aln
#source activate cdtree_env
#./removeInvariant_for_beast.py --in $NUC_ALN --out $OUT_NUC_ALN 
#source deactivate

#' count the constant sites for the gpa alignment
#' edit the datatype in the nexus file manuyally
export GPA_ALN=/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results.5/gpa/gpa.aln
export OUT_GPA_ALN=/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results.5/gpa/gpa.compressed.aln
source activate cdtree_env
./removeInvariant_for_beast.py --in $GPA_ALN --out $OUT_GPA_ALN --type binary
source deactivate

