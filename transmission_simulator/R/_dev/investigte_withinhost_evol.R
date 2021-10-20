closeAllConnections()
rm(list=ls())

## This is the pipeline of transmission tree simulation
## simulate a transmission tree

## param space
## ws.shape= c(10, 1, 20, 15)
## sample.p= c(0.05, 0.5, 0.95, 0.25)
library(TransPhylo)
set.seed(0)
output_d<- '/home/thkuo/bifo-sgi/transmission_simulator.application/set14'

neg=100/365
off.r=2.5
w.shape=1
w.scale=0.01
ws.shape=15
ws.scale=0.01
sample.p=0.95
dateT= 2020
dateStartOutbreak=2019
nSampled = 100

ctree<- simulateOutbreak(nSampled = nSampled, 
    neg=neg,pi=sample.p,off.r=off.r,
    w.shape=w.shape,w.scale=w.scale,
    ws.shape = ws.shape, ws.scale = ws.scale,
    dateStartOutbreak=dateStartOutbreak,dateT=dateT)
##plotCTree(ctree)
##source('wtreesFromCTree.R')
##wtrees<- wtreesFromCTree(ctree)
##
##source('phyloFromWTree_v2.R')
##w_phylos<-lapply(wtrees, function(wtree) phyloFromWTree(wtree))
##
### source('simulateOutbreak_unglued.R')
### t_w_trees <- simulateOutbreak_unglued(
###     neg=neg,pi=pi,off.r=off.r,w.shape=w.shape,
###     w.scale=w.scale,dateStartOutbreak=dateStartOutbreak,dateT=dateT)
### plotTTree(list(ttree= t_w_trees$ttree, nam= t_w_trees$nam), w.shape = w.shape, w.scale = w.scale)
### plotTTree2(list(ttree= t_w_trees$ttree, nam= t_w_trees$nam))
### source('TransPhylo/R/glueTrees.R')
### source('TransPhylo/R/computeHost.R')
### truth<-.glueTrees(t_w_trees$ttree,t_w_trees$wtrees)
### # truth[,1]<-truth[,1]+dateStartOutbreak
### ctree=list(ctree=truth,nam=mtt$nam,probttree=probttree,probwithin=probwithin)
### class(ctree)<-'ctree'
### plotCTree(ctree)
##
#### save the trees
##dir.create(file.path(output_d), recursive = T)   
##library(ape)
#### don't forget the node indices
#### save the parameters
##params_var<- list(neg= neg,
##                  off.r=off.r,
##                  w.shape=w.shape, 
##                  w.scale=w.scale,
##                  ws.shape= ws.shape,
##                  ws.scale= ws.scale,
##                  pi=sample.p, 
##                  dateT= dateT, 
##                  dateStartOutbreak=dateStartOutbreak,
##                  nSampled = nSampled)
##params_f<- sprintf('%s/transphylo_params', output_d)
##write.table(t(data.frame(params_var)), file = params_f, 
##            col.names = F, row.names = T,
##            quote = F, 
##            sep= '\t')
#### print the ctree
##ctree_p_f<- sprintf('%s/ctree.pdf', output_d)
##pdf(ctree_p_f)
##plotCTree(ctree)
##dev.off()
#### the transmission tree
##ttree_f<- sprintf('%s/ttree.mat', output_d)
##ttree<- extractTTree(ctree)
##write.table(ttree$ttree, file = ttree_f, sep = '\t',  col.names = F)
##ttree_p_f<- sprintf('%s/ttree.pdf', output_d)
##pdf(ttree_p_f)
##plotTTree2(ttree)
##dev.off()
#### the (time-measured) phylogeny
##ptree_f<- sprintf('%s/phyloFromPtree.nwk', output_d)
##write.tree(phy = phyloFromPTree(extractPTree(ctree)), file = ptree_f)
#### the within-host trees
##output_wtrees_d<- sprintf('%s/wtrees', output_d)
##dir.create(file.path(output_wtrees_d), recursive = T)   
#### ori index as the tree names
#### so that the trees are traceable (i.e. tip of one to the origin of another)
##for (i in 1:length(w_phylos)){
##    dtr<- w_phylos[[i]]$phylo
##    dtree.name<- w_phylos[[i]]$ori
##    dtr_f<- sprintf('%s/%s.nwk', output_wtrees_d, dtree.name)
##    write.tree(phy= dtr, file = dtr_f)
##}
