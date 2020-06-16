
closeAllConnections()
rm(list=ls())

## This is the pipeline of transmission tree simulation
## simulate a transmission tree
library(argparse)
library(TransPhylo)
set.seed(0)
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
    dest="verbose", help="Print little output")
parser$add_argument("-s", "--simulation_parameters", dest= 'yaml', nargs= 1,
    help="yaml file describing parameters for simulation using TransPhylo")
parser$add_argument("-o", "--output_dir", dest= 'out',
    help="output directory for the output trees")
args <- parser$parse_args()
output_d<- args$out
transphylo_params<- list()
if (!file.exists(args$yaml)){
    stop(sprintf("Specified file ( %s ) does not exist", args$yaml))
}else{
    transphylo_params<- read_yaml(args$yaml)
}
ctree<- simulateOutbreak(transphylo_params)
# extract the within-host subtrees
source('wtreesFromCTree.R')
wtrees<- wtreesFromCTree(ctree)
# convert the transmission trees to phylo
source('phyloFromWTree_v2.R')
w_phylos<-lapply(wtrees, function(wtree) phyloFromWTree(wtree))

## save the trees
dir.create(file.path(output_d), recursive = T)   
library(ape)
## don't forget the node indices
## save the parameters
params_f<- sprintf('%s/transphylo_params', output_d)
write.table(t(data.frame(transphylo_params)), file = params_f, 
            col.names = F, row.names = T,
            quote = F, 
            sep= '\t')
## print the ctree
ctree_p_f<- sprintf('%s/ctree.pdf', output_d)
pdf(ctree_p_f)
plotCTree(ctree)
dev.off()
## the transmission tree
ttree_f<- sprintf('%s/ttree.mat', output_d)
ttree<- extractTTree(ctree)
write.table(ttree$ttree, file = ttree_f, sep = '\t',  col.names = F)
ttree_p_f<- sprintf('%s/ttree.pdf', output_d)
pdf(ttree_p_f)
plotTTree2(ttree)
dev.off()
## the (time-measured) phylogeny
ptree_f<- sprintf('%s/phyloFromPtree.nwk', output_d)
write.tree(phy = phyloFromPTree(extractPTree(ctree)), file = ptree_f)
## the within-host trees
output_wtrees_d<- sprintf('%s/wtrees', output_d)
dir.create(file.path(output_wtrees_d), recursive = T)   
## ori index as the tree names
## so that the trees are traceable (i.e. tip of one to the origin of another)
for (i in 1:length(w_phylos)){
    dtr<- w_phylos[[i]]$phylo
    dtree.name<- w_phylos[[i]]$ori
    dtr_f<- sprintf('%s/%s.nwk', output_wtrees_d, dtree.name)
    write.tree(phy= dtr, file = dtr_f)
}
