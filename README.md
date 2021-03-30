# Improved tree inference and transmission reconstruction from composite datatype representations of microbial genomes

This repository contains the methodology of the phylogenetic tree research:
- `tree_inference` includes workflows of the composite datatype method 
- `transmission_ana` includes the methods of transmission analysis
- `outbreak_simulation` includes the methods and parameters for simulating the dataset used in this research
- `Klebsiella` includes the clustering effect analysis the input parameters for transmission analyses, such as
  BEAST2 xml files
- `Pseudomonas`includes the clustering effect analysis 

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
required by each analysis) are listed as yaml files and controlled by snakemake
and conda. 
