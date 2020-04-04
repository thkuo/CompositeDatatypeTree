library(TransPhylo)
## This is the pipeline of transmission tree simulation
## simulate a transmission tree
set.seed(0)
neg=100/365
off.r=5
w.shape=10
w.scale=0.1
pi=0.75
dateT= 2008
dateStartOutbreak=2005

source('simulateOutbreak_unglued.R')
t_w_trees <- simulateOutbreak_unglued(
    neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
    w.scale=w.scale,dateStartOutbreak=dateStartOutbreak,dateT=dateT)
plotTTree2(list(ttree= t_w_trees$ttree, nam= t_w_trees$nam))

## save the trees
## don't forget the node indices
output_d<- 'outbreak_sim'
dir.create(file.path(output_d), recursive = T)   
## the transmission tree
ttree_f<- sprintf('%s/ttree.mat', output_d)
write.table(t_w_trees$ttree, file = ttree_f, sep = '\t',  col.names = F)
## save the within-host trees
output_wtrees_d<- sprintf('%s/wtrees', output_d)
dir.create(file.path(output_wtrees_d), recursive = T)   
library(ape)
for (patient_ix in 1:length(t_w_trees$wtrees)){
    # patient_ix<- 1
    ptree<- list(ptree= t_w_trees$wtrees[[patient_ix]], nam= t_w_trees$nam)
    dtr<- phyloFromPTree(ptree)
    dtr_f<- sprintf('%s/%s.nwk', output_wtrees_d, patient_ix)
    write.tree(phy= dtr, file = dtr_f)
}