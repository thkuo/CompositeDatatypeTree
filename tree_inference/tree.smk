##For tree inference, the composite datatype method includes:
##1. snps_workflow: map the reads for calling variants
##2. nuc_workflow: infer the nucleotide tree
##3. denovo_workflow: compute de novo assembly, detect genes, and cluster the
##   orthologous 
##4. gpa_workflow: convert a gene orthologous matrix (ie the output table of
##   Roary) into gene presence/absence
##   alignment
##5. conc_workflow: concatenate nucleotide and gpa alignments and conduct

subworkflow mapping:
    output:
        (
            '/net/sgi/metagenomics/data/from_moni/'
            'old.tzuhao/TreePaper/Paeru.vPre3-3/'
            'results/seq2geno/snps/merged_vcf/'
            'multisample.snp.vcf.gz')
    workdir:
        "snps_workflow/"
    snakefile:
        "snps_workflow/snps.in_one.smk"
    configfile:
        ""

rule step1:
    input:
        mapping
    output

# mapping
snakemake -p --cores 16 --use-conda\
  --rerun-incomplete \
  --configfile=snps_config.yml \
  --directory= \
  --snakefile= \
  

# the nucleotide tree
snakemake -p --cores 20 --use-conda \
  --rerun-incomplete \
  --configfile=nuc_workflow/nuc_config.yml \
  --conda-prefix=env/ \
  --directory=nuc_workflow/ \
  --snakefile=nuc_workflow/nuc_tr.smk \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/results/raxml/RAxML_bipartitions.nuc.bs

# compute the de novo assemlies and orthologous families
#--> need to ensure <--#
snakemake -p \
 --rerun-incomplete \
 --conda-prefix /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs \
 --use-conda -j 28 \
 --configfile=denovo_workflow/denovo.config.yml \
 --snakefile=denovo.in_one.smk

## make the gpa alignment
#snakemake \
#  --configfile ./config.yml \
#  --snakefile ./gpa.smk \
#  ../../results.v5/gpa/gpa.var.aln

subworkflow gpa:
    workdir:
        "./gpa_workflow"
    snakefile:
        "./gpa_workflow/gpa.smk"
    configfile:
        "./config.yaml"

rule gpa_rule:
    input:
    output:
        gpa_var_aln=gpa(../../results.v5/gpa/gpa.var.aln),
        

## composite datatype tree
subworkflow cd_tree:
    workdir:
        "./conc_workflow"
    snakefile:
        "./conc_workflow/conc_tree.smk"
    configfile:
        "./conc_workflow/conc.config.yml"

rule cd_tree_rule:
    input:
        gpa_var_aln=../../results.v5/gpa/gpa.var.aln,
        nuc_tree=
        
    output:
        ('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.vPre3-3/results/after_nuc_tree/v1/RAxML_bipartitions.conc.bs')

rule all:
    input:
        cd_tree('/net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.vPre3-3/results/after_nuc_tree/v1/RAxML_bipartitions.conc.bs')
    

#snakemake -p \
# --rerun-incomplete \
# --conda-prefix /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs \
# --use-conda -j 28 \
# --configfile=conc.config.yml \
# --snakefile=
