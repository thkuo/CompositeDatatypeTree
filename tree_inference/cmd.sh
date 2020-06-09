
source activate snakemake_env

#python merge.py nuc_tr \
# --dry \
# --c  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/nuc_workflow/nuc_config.yml \
# --o /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/mapping/raxml/RAxML_bipartitions.nuc.bs

#python merge.py cd_tr \
# --dry \
# --c  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/conc_workflow/conc.config.yml \
# --o '' 

#python merge.py gpa \
# --dry \
# --c  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/gpa_workflow/config.yml \
# --o /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/gpa/gpa.var.aln 

python merge.py mapping \
 --dry \
 --r /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim/data/reference/ATCC_700669.fasta \
 --p /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/ \
 --l /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/run_seq2geno/dna_list

#python merge.py denovo \
# --dry \
# --p /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/seq2geno/denovo \
# --l /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/run_seq2geno/dna_list
