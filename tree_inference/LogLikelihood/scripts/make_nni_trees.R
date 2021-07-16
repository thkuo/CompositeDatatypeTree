# SPDX-FileCopyrightText: 2021 Tzu-Hao Kuo
#
# SPDX-License-Identifier: GPL-3.0-or-later

library(phytools)
library(phangorn)

args<- list('in'= snakemake@input[['nuc_tr']], 'out'= snakemake@output[['rear_trs']])
print(args)
main_tr_f<- args[['in']]
main_tr<- read.newick(main_tr_f) 
out_f<- args[['out']]
rear_trs<- nni(main_tr)
write.tree(c(main_tr, rear_trs), out_f)
