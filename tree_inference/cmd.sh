
source activate snakemake_env
python merge.py mapping \
 --dry \
 --c  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/bin.v5/snps_workflow/config.yml \
 --o /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/mapping/merged_vcf/multisample.snp.vcf.gz

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

#python merge.py denovo \
# --dry \
# --c /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/seq2geno/denovo/denovo_config.yml \
# --o /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/WhichTree_Sim.v7/results.v5/seq2geno/seq2geno/denovo/roary/gene_presence_absence.csv
