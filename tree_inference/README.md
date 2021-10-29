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

### <a name="introduction"></a> Introduction
This tool computes the composite datatype tree from (1) a list of sequencing reads corresponding to each sample and (2) the reference genome sequence.
Considering the complexity, we partition the whole workflow into five parts:
1. snps\_workflow: map the reads for calling variants (output: multi-sample vcf file)
2. nuc\_workflow: infer the nucleotide tree (output: phylogenetic tree in Newick format)
3. denovo\_workflow: compute de novo assembly, detect genes, and cluster the
   orthologous (output: matrix of gene families)
4. gpa\_workflow: convert a gene family table (i.e., the output table of                                                                                                                                           
   Roary) into gene presence/absence alignment (output: alignment of gene presence/absence states)
5. conc\_workflow: concatenate nucleotide and GPA alignments and conduct
   composite datatype inference (output: phylogenetic tree in Newick format)

### <a name="installation"></a> Installation
1. Clone this package to proper locations

2. Install the core environment using conda (tested version: 4.8.0):
```sh
conda env create -f ../envs/cdtree_env.yml
```

3. Some external tools still need to be installed additionally:

- [stampy](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy)
- [roary](https://github.com/hzi-bifo/Roary)
- [prokka](https://github.com/tseemann/prokka): it's dependency _tblasn_ regularly expires. When the annotation procedure has problems, re-installing the package may help.

Except for them, the other process-specific software and environments will be installed (once for the first time) and activated when the process is used. To determine where to create the process-specific environments, please set up the environmental variable `CDTREE_SHARED_ENV`, avoiding installing the same environment in the future:
```sh
export CDTREE_SHARED_ENV=~/bin/cdtree_sharedEnv/
```

- adding the above command in files such as `.bashrc`, which sets up the
  variable at OS login
- adding the above command following [this tutorial of conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables), which sets up the variable when the environment is activated

4. Include the location of this package in the PATH variable. 

### <a name="utility"></a> Utility and input data formats

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
  --config CONFIG_F     YAML file to overwrite default parameter settings
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

3. paths and parameters:

This tool uses Snakemake to resolve the dependencies among processes. Because the process relies on filenames, the intermediate files will be determined automatically. Despite that, there are still some settings that can be customized by setting them in the config YAML file (see *example_config.yml*) with `--config`.

For intermediate file paths, the config file enables setting the multi-sample .vcf.gz file and gene presence/absence table from Roary with the below fields:
```
roary_out (string; existing filename)
multisample_vcf (string; existing filename)
```

The parameters listed below can be specified in the config YAML file:
```
nuc_subs_model (string): the nucleotide substitution model for RAxML
rate_model (string): the rate model for RAxML
outgroup (list of strings): the position of root
cutoff_perc (0.0 to 1.0): the stringency of branch confidence; setting 0.0 will collapse all branches.

```

For the models to use in the tree inference, the two parameters `nuc_subs_model` (default: GTR) and `rate_model` (default: GAMMA) are substition and rate variation models for RAxML. To adjust the settings, please use model test tools such as [jmodeltest](https://github.com/ddarriba/jmodeltest2).

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

