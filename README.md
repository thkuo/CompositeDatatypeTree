<!--
SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo

SPDX-License-Identifier: GPL-3.0-or-later
-->

# Improved phylogenetic tree inference and transmission reconstruction from composite datatype representations of microbial genomes

---
Overview: 
1. [Contents](#contents)
2. [Prerequisities](#prerequisite)
3. [Contact](#contact)
4. [Licenses](#licenses)
---

### <a name="contents"></a> Contents

This repository contains the methodology and results of this paper. For the results of each application, please visit the `bin/paper_figs.html` file under `Pseudomonas/`, `Klebsiella/`, or `Benchmarking/` (simulated dataset). For the methods, the folder `tree_inference/` includes the setup and usage of the composite datatype inference, and `transmission_simulator/` contains methods for the simulation in an outbreak scenario (please also refer to `Benchmarking/`). 

    .
    ├── tree_inference # workflows of the composite datatype method 
    ├── transmission_simulator #the scripts for simulating the phylogenetic tree, genomes, and gene presence or absence patterns in the scenario of outbreak
    ├── envs/ # yaml files for conda environments
    ├── Pseudomonas/ # methods and results of the Pseudomonas aeruginosa dataset
        ├── bin/
            ├── tree_inference.sh # the commands for tree inference
            └── paper_figs.html # the visualizations with R
        ├── data/ # sampling locations
        └── results/
          ├── nuc.var.aln # the nucleotide alignment
          ├── gpa.var.aln # the gene presence/absence alignment
          ├── nuc_tr.nwk # the nucleotide tree
          └── cd_tr.nwk # the composite datatype tree
    ├── Klebsiella/ # methods and results of the Klebsiella pneumoniae dataset
        ├── bin/
            ├── tree_inference.sh # the commands for tree inference
            ├── paper_figs.html # the visualizations with R
            ├── ana_var/ # early analysis about the divergence among samples
            ├── transmission_ana/ # the BEAST2 and TransPhylo methods
            └── nuc_permutation/ # the methods of nucleotide replacement test
        └── results/
            ├── nuc.var.aln # the nucleotide alignment
            ├── gpa.var.aln # the gene presence/absence alignment
            ├── nuc_tr.nwk # the nucleotide tree
            ├── cd_tr.nwk # the composite datatype tree
            ├── ana_var/ # early analysis about the divergence among samples
            ├── transmission_ana/ # the BEAST2 and TransPhylo results
            └── nuc_permutation/ # the results from nucleotide replacement test
    └── Benchmarking/ # methods and results of simulating and analysing the Streptococcus pneumoniae dataset
        ├── bin/
	    ├──	tree_inference.sh # the commands for tree inference
            ├── paper_figs.html # the visualizations with R
            ├── cutoffs/ # methods and results from the iteration with different cutoffs of log-lieklihood scores
            └── simulate/ # commands for simulating the dataset 
        ├── data/
            ├── outbreak_sim/ # the tree simulation results using TransPhylo
            └── pseudo_gpa/ # the simulation result of gene presence or absence 
        ├── config # parameters and input settings for the simulation
        └── results/
            ├── nuc.var.aln # the nucleotide alignment
            ├── nuc_tr.nwk # the nucleotide tree
            └── ll_cutoff/ # the composite datatype trees with different cutoffs of log-likelihodd scores

### <a name="prerequisite"></a> Prerequisites

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

### <a name="contact"></a> Contact
Tzu-Hao Kuo Tzu-Hao.Kuo@helmholtz-hzi.de


### <a name="licenses"></a> Licenses
Please refer to the copyright header of each file or `LICENSES/`. 
