#!/usr/bin/env Rscript
closeAllConnections()
rm(list=ls())

#== NOTE
##Collapse branches of raw output from RAxML
##
##branch length
##bootstrap value
##outgroup???

collapse_branches<- function(tr, target_nodes){
    library(ape)
    br<- tr$edge.length
    tol_value<- min(br[br>0])
    
    new_br<- br
    new_br[tr$edge[,2] %in% target_nodes]<- 0
    col_tr<- di2multi(
        compute.brlen(tr,new_br), 
        tol = tol_value)    
    
    return(col_tr)    
}

find_short_branches<- function(tr, br_cutoff){
    br<- tr$edge.length
    bad_nodes<- tr$edge[br<=br_cutoff, 2]
    return(bad_nodes)
}

find_low_branches<- function(tr, bs_cutoff){
    library(ggtree)
    tr_info<- fortify(tr)
    sub_tr_info<- tr_info[!tr_info$isTip & (tr_info$parent!=tr_info$node), ]
    bad_nodes<- sub_tr_info[as.numeric(sub_tr_info$label) <= bs_cutoff, ]$node
    bad_nodes<- bad_nodes[!is.na(bad_nodes)]
    return(bad_nodes)
}

write.cutoff<-function(cutoff, cutoff_f){
    library(reshape2)
    write.table(melt(as.data.frame(cutoff)), 
                file= cutoff_f, 
                sep = '\t', 
                col.names = F, row.names = F, 
                quote = F)
}

library(phytools)
library(argparse)
w_dir<- getwd()
print(w_dir)
parser <- ArgumentParser(description='Collapse branches with cutoffs of length and bootstrap support')
parser$add_argument('--i', type="character", dest= 'tr_f',
                    help='input tree')
parser$add_argument('--o', type="character", dest= 'output_tr_f', 
                    help='output tree')
parser$add_argument('--br', dest='cutoff.br', default= 1e-4,
		    type= 'double',
                    help='the cutoff for branch length')
parser$add_argument('--bs', dest='cutoff.bs', default= 60,
		    type= 'integer',
                    help='the cutoff for bootstrap support')
parser$add_argument('--og', type="character", nargs= '+', 
		    dest= 'outgroup',
                    help='outgroup for rooting')

args <- parser$parse_args()
tr_f<- args$tr_f
cutoff<-list(br=args$cutoff.br , bs= args$cutoff.bs)
outgroup<- args$outgroup
output_tr_f<- args$output_tr_f
cutoff_f<- sprintf('%s_CUTOFF', output_tr_f)

# nucleotide tree
tr<- unroot(read.newick(tr_f))
tr<- phytools::reroot(tr, 
		      node.number = ifelse(length(outgroup)>1, getMRCA(tr, outgroup), which(tr$tip.label == outgroup[1])),
		      position = 1e-6)

#source('collapse_branches.R')
short_ix<- find_short_branches(tr, br_cutoff = cutoff$br)
low_ix<- find_low_branches(tr, bs_cutoff = cutoff$bs)
col_tr<- collapse_branches(tr= tr, target_nodes = union(short_ix, low_ix))
col_tr<- phytools::reroot(col_tr, 
		 node.number = ifelse(length(outgroup)>1, getMRCA(col_tr, outgroup), which(col_tr$tip.label == outgroup[1])),
		 position = 1e-6)

# write results
dir.create(dirname(output_tr_f), recursive = T)
write.cutoff(cutoff, cutoff_f)
write.tree(col_tr, file = output_tr_f)
