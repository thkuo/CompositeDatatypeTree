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
plotTTree(list(ttree= t_w_trees$ttree, nam= t_w_trees$nam), w.shape = w.shape, w.scale = w.scale)
plotTTree2(list(ttree= t_w_trees$ttree, nam= t_w_trees$nam))
source('TransPhylo/R/glueTrees.R')
source('TransPhylo/R/computeHost.R')
truth<-.glueTrees(t_w_trees$ttree,t_w_trees$wtrees)
# truth[,1]<-truth[,1]+dateStartOutbreak
ctree=list(ctree=truth,nam=mtt$nam,probttree=probttree,probwithin=probwithin)
class(ctree)<-'ctree'
plotCTree(ctree)

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
source('phyloFromWTree.R')
for (patient_ix in 1:length(t_w_trees$wtrees)){
    dtr<- phyloFromWTree(t_w_trees$wtrees[[patient_ix]])
    dtr_f<- sprintf('%s/%s.nwk', output_wtrees_d, patient_ix)
    write.tree(phy= dtr, file = dtr_f)
}
