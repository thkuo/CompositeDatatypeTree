# Inference of Composite Datatype Tree 

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

### Installation
1. Clone this package to proper locations

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

4. Include the location of this package in the PATH variable. 

### Customization:

This tools is based on snakemake, which means the dependencies among processes
(also called "rules" in snakemake) are determined based on filenames. The user
usually do not need to worry about how to name each of the huge number of
files---from reads mapping to de novo assemblies---but besides the input data,
some intermediate files and parameters can still be customized.

#### Path parameters

This enables precomputed files from reads mapping and de novo assembly results; precisely, the multi-sample .vcf.gz file and gene presence/absence table from Roary. They can be specied in the config yaml file:
```
roary_out (string; existing filename)
multisample_vcf (string; existing filename)
```

#### Non-path parameters

The parameters (correct datatype needed) below can be specified in the config yaml file:
```
nuc_subs_model (string)
rate_model (string)
outgroup (list of strings)
br_cutoff (float)
bs_cutoff (integer)
```

The two parameters `nuc_subs_model` (default: GTR) and `rate_model` (default: GAMMA) are substition and rate variation models for RAxML. They can be tested using external software jmodeltest.


#### Utility

1. arguments
```
usage: cdtree [-h] [--config CONFIG_F] [--cpu CPU] [--dry]
              project_dir list_f ref {mapping,nuc_tr,cd_tr,gpa,denovo,all}

Tree inference using composite datatype of microbial representations

positional arguments:
  project_dir           the directory for project
  list_f                the list of sequencing reads
  ref                   the reference genome for mapping sequencing reads
  {mapping,nuc_tr,cd_tr,gpa,denovo,all}
                        the function to launch

optional arguments:
  -h, --help            show this help message and exit
  --config CONFIG_F     yaml file to overwrite default parameter settings
  --cpu CPU             cpu number; default: 1
  --dry                 display the processes and exit
```

For the arguments about files:

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

- CONFIG_F

An .yaml file where the customized parameters are listed. For more details about customization,
please refer to the above section. 

2. examples

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


