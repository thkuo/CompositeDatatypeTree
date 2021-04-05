# Improved tree inference and transmission reconstruction from composite datatype representations of microbial genomes

This repository contains the methodology of the phylogenetic tree research:
- `transmission_simulator`: the scripts for simulating the phylogenetic tree, genomes, and gene presence or absence patterns in the scenario of outbreak
- `tree_inference`: workflows of the composite datatype method 
- `envs/`: yaml files for conda environments
- `Benchmarking/`: the parameters for simulating _Streptococcus pneumoniae_ dataset and the results
  - `bin/`:
    - `cutoffs/`: methods and results from the iteration with different cutoffs of log-lieklihood scores
    - `simulate/`: commands for simulating the dataset 
  - `config`: parameters and input settings for the simulation
  - `data/`:
    - `outbreak_sim/`: the tree simulation results using TransPhylo
    - `pseudo_gpa/`: the simulation result of gene presence or absence 
  - `results/`
    - `nuc_tr.nwk`: the nucleotide tree
    - `ll_cutoff/`: the composite datatype trees with different cutoffs of
      log-likelihodd scores
- `Pseudomonas/`: the results and visualizing scripts for the _Pseudomonas aeruginosa_ dataset
    - `paper_figs.Rmd`: the visualizations with R
    - `nuc.var.aln`: the nucleotide alignment
    - `gpa.var.aln`: the gene presence/absence alignment
    - `nuc_tr.nwk`: the nucleotide tree
    - `cd_tr.nwk`: the composite datatype tree
- `Klebsiella`: methods and results of the _Klebsiella pneumoniae_ dataset
  - `results/`: 
    - `nuc.var.aln`: the nucleotide alignment
    - `gpa.var.aln`: the gene presence/absence alignment
    - `nuc_tr.nwk`: the nucleotide tree
    - `cd_tr.nwk`: the composite datatype tree
    - `transmission_ana/`: the BEAST2 and TransPhylo results
    - `nuc_permutation/`: the results from nucleotide replacement test
  - `bin/`:
    - `paper_figs.Rmd`: the visualizations with R
    - `transmission_ana/`: the BEAST2 and TransPhylo methods
    - `nuc_permutation/`: the methods of nucleotide replacement test

# Prerequisites
    - [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) (used version: 4.8.4)
    - file [.condarc](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-channels.html) that includes these channels and is detectable by your conda
      - hzi-bifo
      - conda-forge/label/broken
      - bioconda
      - conda-forge
      - defaults
    - [python](https://www.python.org/downloads/) (used verson: 3.6)
    - [Linux](https://www.cyberciti.biz/faq/find-linux-distribution-name-version-number/) (used version: Debian GNU/Linux 8.8 jessie)
    - [git](https://git-scm.com/downloads) (used version: 2.18)
    - [stampy](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy) (used version: v1.0.23 (r2059))

The other analysis-wise dependencies (that is, the software and environment
required by each analysis) are listed in yaml files and controlled by snakemake
and conda. Those yaml files can be found in `envs/` and used to build
environments with conda ([tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)).
