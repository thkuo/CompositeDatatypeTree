<!--
SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo

SPDX-License-Identifier: GPL-3.0-or-later
-->

# Inference of Composite Datatype Tree 

---
- [Introduction](#introduction)
- [Installation](#installation)
- [Utility and input data formats](#utility) 
---

###<a name="introduction"></a> Introduction
This tool computes composite datatype tree from (1) a list of sequencing reads corresponding to each sample and (2) the reference genome sequence.
Considering the complexity, we partition the whole workflow into five parts:
1. snps\_workflow: map the reads for calling variants (output: multi-sample vcf file)
2. nuc\_workflow: infer the nucleotide tree (output: phylogenetic tree in newick format)
3. denovo\_workflow: compute de novo assembly, detect genes, and cluster the
   orthologous (output: matrix of gene families)
4. gpa\_workflow: convert a gene orthologous matrix (ie the output table of
   Roary) into gene presence/absence alignment (output: alignment of gene presence/absence states)
5. conc\_workflow: concatenate nucleotide and gpa alignments and conduct
   composite datatype inference (output: phylogenetic tree in newick format)

###<a name="installation"></a> Installation
1. Clone this package to proper locations

2. Install the environment using conda (tested version: 4.8.0)
```sh
conda env create -f installation/cdtree_env.yml
```

3. Some external tools still need to be installed additionally:

- stampy: downloadable via the [ofiicial site](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy)
- prokka: it's dependency _tblasn_ regularly expires, so please ensure the [latest version](https://github.com/tseemann/prokka)
- roary: tested version [here](https://github.com/hzi-bifo/Roary)

Except for them, the other process-specific software and environments will be installed (once for the first time) and activated when the process is used. Set up the environmental variable `CDTREE_SHARED_ENV`, which determines where the process-specific environment should be built. This will avoid repetive installation of same environment in the future. For quick and temporary utility, we recommend to try:

```sh
export CDTREE_SHARED_ENV=~/bin/cdtree_sharedEnv/
```

For long-term management, you might want to:
- adding the above command in files such as `.bashrc`, which sets up the
  variable at OS loggin
- adding the above command following [this tutorial of conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables), which sets up the variable when the environment is activated

4. Include the location of this package in the PATH variable. 

###<a name="utility"></a> Utility and input data formats

1. arguments
```
usage: cdtree [-h] [--config CONFIG_F] [--cpu CPU] [--dry]
              project_dir list_f ref
              {mapping,fast_mapping,nuc_tr,col_tr,denovo,gpa,cd_tr,all}

positional arguments:
  project_dir           the directory for project
  list_f                the list of paired-end DNA sequencing reads
  ref                   the reference genome for read mapping
  {mapping,fast_mapping,nuc_tr,col_tr,denovo,gpa,cd_tr,all}
                        the workflow to launch

optional arguments:
  -h, --help            show this help message and exit
  --config CONFIG_F     yaml file to overwrite default parameter settings
  --cpu CPU             cpu number; default: 1
  --dry                 display the processes and exit
```

2. input formats

- list_f

It should be a two-column list, where the first column includes all samples and the second column lists the __paired-end reads files__. The two reads file are separated by a comma. The first line is the first sample.
```
sample01	/paired/end/reads/sample01_1.fastq.gz,/paired/end/reads/sample01_2.fastq.gz
sample02	/paired/end/reads/sample02_1.fastq.gz,/paired/end/reads/sample02_2.fastq.gz
sample03	/paired/end/reads/sample03_1.fastq.gz,/paired/end/reads/sample03_2.fastq.gz
```

- ref

The fasta file of reference genome, which will be used in the mapping-relevant
procedures.

3. customization:

This tool uses Snakemake to resolve the dependencies among processes. Because the process relies on filenames, the intermediate files will be determined automatically. Despite that, there are still some settings that can be customized by setting them in the config file (see *example_config.yml*) with `--config`.

For the models to use in the tree inference, the two parameters `nuc_subs_model` (default: GTR) and `rate_model` (default: GAMMA) are substition and rate variation models for RAxML. To adjust the settings, please use external software [jmodeltest](https://github.com/ddarriba/jmodeltest2) and update the config file.

4. examples

```sh
# activate the environment
source activate cdtree_env
# ensure the location to built and access the process-specific environments
export CDTREE_SHARED_ENV=/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/shared_envs
 
./cdtree --cpu 25 --dry\
  ./test \
  some/dna_reads_list \
  some/reference/genome_seq.fasta \
  all 
```

