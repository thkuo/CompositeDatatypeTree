d_tree_files<- list()
sample_time_files<- list()
for (clade in c('clade1', 'clade2', 'clade3')){
    for (datatype in c('nuc', 'md')){
        x<- paste(clade, datatype, sep =  '_')
        d_tree_files[[x]]<- paste('/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results/after_nuc_tree/v1', 
                                  clade, 
                                  datatype, 
                                  'model_test.best', 
                                  paste(datatype, 'd_tree', sep= '.'), 
                                  sep= '/')
        sample_time_files[[x]]<- paste('/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v2/results/after_nuc_tree/v1', 
                                       clade,
                                       'sample_time_subset.tsv',
                                       sep= '/')
    }
}                  

library(ape)
library(TransPhylo)
# tr<- read.nexus(d_tree_f)
# tr<- compute.brlen(tr, tr$edge.length/365)
trs<- lapply(d_tree_files, function (f) read.nexus(f))
trs<- lapply(trs, function(tr) compute.brlen(tr, tr$edge.length/365))

compute_ptree<- function(t_f, tr){
  print('compute ptree')
  t_df<- read.csv(t_f, sep= '\t', header = F)
  colnames(t_df)<- c('name', 'time')
  t_base<- 2013
  t_df$time<- t_df$time/365 +t_base
  last_sample_date<- max(t_df$time)

  ptreeFromPhylo(tr, dateLastSample=last_sample_date)
}
# ptree<- compute_ptree(sample_time_f, tr)
ptrees<- list()
for (x in names(trs)){
    ptrees[[x]]<- compute_ptree(sample_time_files[[x]], 
                                trs[[x]])
}

w.shape<- 1
#w.scale<- 0.5^(-1)
w.scale<- 0.5
dateT<- Inf
#mcmc_length<- 1e6
# mcmc_length<- 2e5
mcmc_length<- 4e5
# ttree<- inferTTree(ptree,
#               mcmcIterations=mcmc_length,
#               w.shape=w.shape,
#               w.scale=w.scale,
#               dateT=dateT)
# cons_tree<- consTTree(ttree)
set.seed(123)
library(parallel)
records<- mclapply(ptrees, function(ptree) inferTTree(ptree,
                                                 mcmcIterations=mcmc_length,
                                                 w.shape=w.shape,
                                                 w.scale=w.scale,
                                                 dateT=dateT), mc.cores = 6)
library(coda)
print(sapply(records, function(r) effectiveSize(convertToCoda(r))))
cons_trees<- mclapply(records, 
		      function(rec) consTTree(rec), 
		      mc.cores = 6)

## load the clinical records
bio_f<- '/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/MRKP.v1/data/samples_info/SuppTab1.csv'
bio_df<- read.table(bio_f, sep = ',', header= T, stringsAsFactors = F)
bio_df$Ward<- as.character(bio_df$Ward)
bio_df$W<- as.character(bio_df$Ward)
rownames(bio_df)<- bio_df$Sequencing.Short.NO.

####
## save the results
iteratively_create_d<- function(d){
   if (dir.exists(dirname(d))){
     dir.create(d)
   }else{
    iteratively_create_d(dirname(d))
    dir.create(d)
  }
}
source('./plotTTree3.R')
library(jcolors)
out_d<- 'transmission_results.v7/'
# determine the appearance of each node
available_cols<- as.character(jcolors('pal2')[1:5])
available_shapes<- c( 8, 13, 15, 19)
all_hospitals<- unique(bio_df$Ward)
all_hospitals<- all_hospitals[order(all_hospitals)]
cols<- available_cols[(1:length(all_hospitals) %% length(available_cols)) + 1]
names(cols)<- all_hospitals 
shapes<- available_shapes[(1:length(all_hospitals) %% length(available_shapes)) + 1]
names(shapes)<- all_hospitals 
for (tr_name in names(cons_trees)){
    cons_tree<- cons_trees[[tr_name]]

    d<- paste0(out_d, tr_name)
    iteratively_create_d(d)   

    # determine the color and shape
    sample.col<- as.character(cols[bio_df[cons_tree$nam, 'Ward']])
    sample.shape<- as.numeric(shapes[bio_df[cons_tree$nam, 'Ward']] )
    epi_tr_f<- file.path(d, paste0(tr_name, '_', 'epi_tree.pdf'))
    pdf(epi_tr_f)
    plotTTree3(cons_tree, 
	       sample.col= sample.col, sample.shape= sample.shape,
	       col.legend= cols, shape.legend= shapes)
    dev.off()
}
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

####
## Analyze the transmission tree like a graph instead of a phylo object
ctree_to_edge<- function(cons_tree){
    library(igraph)
    edge<- matrix(c(cons_tree$ttree[,3], 1:nrow(cons_tree$ttree)), ncol= 2)
    edge[which(edge[,1] == 0), 1]<- edge[which(edge[,1] == 0), 2]
    edge
}
ctree_to_graph<- function(cons_tree){
    edge<- ctree_to_edge(cons_tree)
    # add the weights, which are the time interval of each infection time
    g<- graph_from_edgelist(edge, directed= F)
    t<- cons_tree$ttree[,1]
    E(g)$weight<- t[edge[,2]]  - t[edge[,1]]    
    g
}
edges<- lapply(cons_trees, function(ctree) ctree_to_edge(ctree))

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

# the contact pairs
within_ward_infection<- function(edge, cons_tree){
    wiw_mat<- t(apply(edge, 1, 
                      function(l) cons_tree$nam[l]))
    wiw_filter<- apply(wiw_mat, 1, function(l) ! any(is.na(l)))
    wiw_mat_subset<- wiw_mat[wiw_filter, ]
    # if the direct transmission pair is inter-ward or within-pair
    summary(bio_df[wiw_mat_subset[,2], 'Ward'] == bio_df[wiw_mat_subset[,1], 'Ward'])
}
print(sapply(names(edges), function(x) within_ward_infection(edges[[x]], cons_trees[[x]])))
