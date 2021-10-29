<!--
SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo

SPDX-License-Identifier: GPL-3.0-or-later
-->

The two parts for simulating seqeuncing reads of outbreak
- `R/`: for transmission tree and within-host subtrees simulation
- `bin/`: for genome and reads simulation

Prerequisite:
- r:
    - packages described in ../envs/r35.yml
    - TransPhylo (tested version: 1.2.3)

- python:
    - packages described in ../envs/py37.yml

- [ALF](http://alfsim.org/#index)
- [dawg](https://github.com/reedacartwright/dawg)
- [perl scripts](https://github.com/johnlees/which_tree) for integrating the simulated sequences
- [pIRS](https://github.com/galaxy001/pirs)

Using Conda to install the packages listed in the YAML files can be easier
([tutorial](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)).
