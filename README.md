## Improved tree inference and transmission reconstruction from composite datatype representations of microbial genomes

---
To-include:
- jmodeltest (for determing the model of RAxML)
- outbreak simulation
---
`tree_inference` includes workflows that are involved in the composite datatype method. Those methods were applied to both the clinical and benchmark datasets. `Klebsiella`, `Pseudomonas`, and `outbreak_simulation` include the other relevant analysis also stated in the paper.

For tree inference, the composite datatype method includes:
1. snps_workflow: map the reads for calling variants
2. nuc_workflow: infer the nucleotide tree
3. denovo_workflow: compute de novo assembly, detect genes, and cluster the
   orthologous 
4. gpa_workflow: convert a gene orthologous matrix (ie the output table of
   Roary) into gene presence/absence
   alignment
5. conc_workflow: concatenate nucleotide and gpa alignments and conduct
   composite datatype inference

The substitution models for RAxML (step 2 and 5) are determined the alignments (precisely, the non-redundant alignment) by jmodeltest:
```shell
OUT_F=/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/results/alignment/nuc.var.jmodel.out
ALN_F=/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/results/alignment/nuc.var.aln 
SYM_ALN_F=./tmp_aln
ln -s $ALN_F $SYM_ALN_F 
java -jar jModelTest.jar \
  -f -i -g 4 -tr 13 -AIC \
  -d $SYM_ALN_F \
  -o $OUT_F
```
After the tree inferences, branches are collapsed with branch lengths and
bootstrap cutoffs by the r script:
```shell
call the r script that collapses branches
```

Usage:
```shell
export PATH=/home/thkuo/seq2geno2pheno/seq2geno/snps/lib:$PATH

# mapping
snakemake -p --cores 16 --use-conda\
  --rerun-incomplete \
  --configfile=snps_config.yml \
  --directory=snps_workflow/ \
  --snakefile=snps_workflow/snps.in_one.smk \
  /net/sgi/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.vPre3-3/results/seq2geno/snps/merged_vcf/multisample.snp.vcf.gz

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

# make the gpa alignment
snakemake \
  --configfile ./config.yml \
  --snakefile ./gpa.smk \
  ../../results.v5/gpa/gpa.var.aln


## composite datatype tree
snakemake -p \
 --rerun-incomplete \
 --conda-prefix /net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs \
 --use-conda -j 28 \
 --configfile=conc.config.yml \
 --snakefile=conc_tree.smk
```
