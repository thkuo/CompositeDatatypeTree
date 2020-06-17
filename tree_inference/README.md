## Improved tree inference and transmission reconstruction from composite datatype representations of microbial genomes

`tree_inference` includes workflows that are involved in the composite datatype method. Those methods were applied to both the clinical and benchmark datasets. `Klebsiella`, `Pseudomonas`, and `outbreak_simulation` include the other relevant analysis also stated in the paper.

For tree inference, the composite datatype method includes:
1. snps\_workflow: map the reads for calling variants
2. nuc\_workflow: infer the nucleotide tree
3. denovo\_workflow: compute de novo assembly, detect genes, and cluster the
   orthologous 
4. gpa\_workflow: convert a gene orthologous matrix (ie the output table of
   Roary) into gene presence/absence
   alignment
5. conc\_workflow: concatenate nucleotide and gpa alignments and conduct
   composite datatype inference

# Installation
1. Download the package to proper locations
2. Install the environment using conda (tested version: 4.8.0)
```sh
conda env create -f installation/cdtree_env.yml
```

3. Set up the environmental variable `CDTREE_SHARED_ENV`, which determines where the process-specific environment should be built. This will avoid repetive installation of same environment in the future. For quick and temporary utility, we recommend to try:

```sh
export CDTREE_SHARED_ENV=~/bin/cdtree_sharedEnv/
```

For long-term management, two methods may help:
- adding the above command in files such as `.bashrc`, which sets up the
  variable at OS loggin
- adding the above command following [this tutorial of conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables), which sets up the variable when the environment is activated

# Customization:

This tools is based on snakemake, which means the dependencies among processes
(also called "rules" in snakemake) are determined based on filenames. The user
usually do not need to worry about how to name each of the huge number of
files---from reads mapping to de novo assemblies---but besides the input data,
some intermediate files and parameters can still be customized.

- Path parameters

This enables precomputed files from reads mapping and de novo assembly results; precisely, the multi-sample .vcf.gz file and gene presence/absence table from Roary. They can be specied in the config yaml file:
```
roary_out (string; existing filename)
multisample_vcf (string; existing filename)
```

- Non-path parameters

The parameters (correct datatype needed) below can be specified in the config yaml file:
```
nuc_subs_model (string)
rate_model (string)
outgroup (list of strings)
br_cutoff (float)
bs_cutoff (integer)
```

The two parameters `nuc_subs_model` (default: GTR) and `rate_model` (default: GAMMA) are substition and rate variation models for RAxML. They can be tested using external software jmodeltest.


# Utility



