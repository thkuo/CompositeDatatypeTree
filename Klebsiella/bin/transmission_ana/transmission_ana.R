library(ape)
library(TransPhylo)
library(argparse)

parser <- ArgumentParser(description='Run transmission analysis with TransPhylo')
parser$add_argument('-t', type="character",
                   help="the filename of dated tree")
parser$add_argument('-d', type="character", 
                   help="the filename of dates")
parser$add_argument('-p', type="character", 
                   help="project name that will be used as the prefix of output data")
parser$add_argument('-o', type="character",
                   help="output directory")

args<- parser$parse_args()
#' input/output
d_tree_f<- args$t
sample_date_f<- args$d
out_d<- args$o
#d_tree_f<- '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results.5/transmission_ana/nuc.15/nuc_dtree.nex'
#sample_date_f<- '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results.5/transmission_ana/sample_time_subset.tsv'
#out_d<- '../../results.5/transmission_ana/nuc.15/transphylo'
tr_name<- args$p
#tr_name<- 'nuc'
dir.create(out_d, recursive = TRUE)
#' the intermediate data of RData
#' that are needed for optimizing the visualization
ttree_rds_f<- sprintf('%s/%s_ttree.Rds', out_d, tr_name)
ctree_rds_f<- sprintf('%s/%s_ctree.Rds', out_d, tr_name)

#' read the dated tree
tr<- read.nexus(d_tree_f)

#' read the time information
t_df<- read.csv(sample_date_f, sep= '\t', header = F)
colnames(t_df)<- c('name', 'time')
#' numerize the dates (units: years)
dates<- apply(t_df, 1, function(r) strptime(r['time'],  format = "%d/%m/%Y"))
dates_num<- lapply(dates, function(d) 1900+d$year+d$yday/365)
last_sample_date<- max(unlist(dates_num))

#' the input format for TransPhylo
set.seed(123)
ptree<- ptreeFromPhylo(tr, dateLastSample=last_sample_date)

#' the transmission tree (MCMC procedure)
w.shape<- 1
w.scale<- 0.5
dateT<- Inf
mcmc_length<- 1e5
ttree<- inferTTree(ptree,
    mcmcIterations=mcmc_length,
    w.shape=w.shape,
    w.scale=w.scale,
    dateT=dateT)
saveRDS(ttree, file = ttree_rds_f)

##' compute the ESS
#library(coda)
#effectiveSize(convertToCoda(ttree))

#' compute the consensus transmission tree
ctree<- consTTree(ttree)
saveRDS(ctree, file = ctree_rds_f)

#' visualize the results
## the clinical records
bio_f<- '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v1/data/samples_info/SuppTab1.csv'
bio_df<- read.table(bio_f, sep = ',', header= T, stringsAsFactors = F)
bio_df$Ward<- as.character(bio_df$Ward)
bio_df$W<- as.character(bio_df$Ward)
rownames(bio_df)<- bio_df$Sequencing.Short.NO.

#' save the results
source('../../../MM_for_paper/Klebsiella/transmission_reconstruction/plotTTree3.R')
library(jcolors)
make_symbols<- function(bio_df){
    #' determine the appearance of each node
    available_cols<- as.character(jcolors('pal2')[1:5])
    available_shapes<- c( 8, 13, 15, 19)
    all_locs<- unique(bio_df$Ward)
    all_locs<- all_locs[order(all_locs)]

    #' the symbols assigned to each
    cols<- available_cols[(1:length(all_locs) %% length(available_cols)) + 1]
    names(cols)<- all_locs 
    shapes<- available_shapes[(1:length(all_locs) %% length(available_shapes)) + 1]
    names(shapes)<- all_locs 
    return(list(cols= cols, shapes= shapes))
}

#' the colors and shapes of nodes
d<-sprintf('%s/%s', out_d, tr_name) 
dir.create(d, recursive= T)

#' determine the color and shape
symbols<- make_symbols(bio_df)
sample.col<- as.character(symbols$cols[bio_df[ctree$nam, 'Ward']])
sample.shape<- as.numeric(symbols$shapes[bio_df[ctree$nam, 'Ward']] )
epi_tr_f<- file.path(d, paste0(tr_name, '_', 'epi_tree.pdf'))
pdf(epi_tr_f)
plotTTree3(ctree,
	   sample.col= sample.col, sample.shape= sample.shape,
	   col.legend= symbols$cols, shape.legend= symbols$shapes)
dev.off()
# low memory method
# library(coda)
# cons_trees<- list()
# for (x in names(ptrees)){
#     print(x)
#     ptree<- ptrees[[x]]
#     ttree<- inferTTree(ptree,
#                        mcmcIterations=mcmc_length,
#                        w.shape=w.shape,
#                        w.scale=w.scale,
#                        dateT=dateT)
#     cons_trees[[x]]<- consTTree(ttree)
#     ess<- effectiveSize(convertToCoda(ttree))
#     print(ess)
# }

####
## TransPhylo functions
#library(viridis)
#wiw_probs<- lapply(records, function(r)computeMatWIW(r))
#wiwp_mat<- wiw_probs$clade3_md
#strain_order<- order(bio_df[colnames(wiwp_mat), 'Ward'])
#wiwp_mat<- wiwp_mat[strain_order, strain_order]
#pheatmap(wiwp_mat, cluster_rows = F, cluster_cols = F,
#         annotation_row = bio_df[, c('Ward', 'Campus')], 
#         color = viridis(20, option = 'C'))
#for (x in 1:ncol(wiwp_mat)-1){
#  for (y in (x+1):nrow(wiwp_mat)){
#    	strain1<- colnames(wiwp_mat[x])
#    	strain2<- colnames(wiwp_mat[y])
#	class<- ifelse(bio_df[strain1, 'Ward'] == bio_df[strain2, 'Ward'], 'within', 'inter')
#	p<- wiwp_mat[strain1, strain2]
#  }
#}

#####
### Analyze the transmission tree like a graph instead of a phylo object
#ctree_to_edge<- function(cons_tree){
#    library(igraph)
#    edge<- matrix(c(cons_tree$ttree[,3], 1:nrow(cons_tree$ttree)), ncol= 2)
#    edge[which(edge[,1] == 0), 1]<- edge[which(edge[,1] == 0), 2]
#    edge
#}
#ctree_to_graph<- function(cons_tree){
#    edge<- ctree_to_edge(cons_tree)
#    # add the weights, which are the time interval of each infection time
#    g<- graph_from_edgelist(edge, directed= F)
#    t<- cons_tree$ttree[,1]
#    E(g)$weight<- t[edge[,2]]  - t[edge[,1]]    
#    g
#}
#edges<- lapply(cons_trees, function(ctree) ctree_to_edge(ctree))

## compute the pairwise shortest path
#d<- shortest.paths(g)
#colnames(d) <- c(cons_tree$nam, as.character((length(cons_tree$nam)+1):ncol(d)))
#rownames(d) <- c(cons_tree$nam, as.character((length(cons_tree$nam)+1):ncol(d)))
#library(pheatmap)
#library(viridis)
#pheatmap(d[1: length(cons_tree$nam), 1: length(cons_tree$nam)], 
#         annotation_col = bio_df[, c('Ward', 'Campus')], 
#         color = viridis(20),
#         treeheight_row = 0, treeheight_col = 0)

## the contact pairs
#within_ward_infection<- function(edge, cons_tree){
#    wiw_mat<- t(apply(edge, 1, 
#                      function(l) cons_tree$nam[l]))
#    wiw_filter<- apply(wiw_mat, 1, function(l) ! any(is.na(l)))
#    wiw_mat_subset<- wiw_mat[wiw_filter, ]
#    # if the direct transmission pair is inter-ward or within-pair
#    summary(bio_df[wiw_mat_subset[,2], 'Ward'] == bio_df[wiw_mat_subset[,1], 'Ward'])
#}
#print(sapply(names(edges), function(x) within_ward_infection(edges[[x]], cons_trees[[x]])))
