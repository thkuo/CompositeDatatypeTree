closeAllConnections()
rm(list=ls())
source('/home/thkuo/projects/TreePaper/bin/collapse_branches.R')
source('entropy_ana.functions.R')
library(ggtree)
library(phytools)
library(entropy)
library(reshape2)

# load data
## the experimental data
source('load_biological_data.R')
bio_df<- load_biological_data()

# and the trees
source('load_outgroup.R')

tr_files<- list('nuc'= '/home/thkuo/bifo-sgi/TreePaper/Paeru.v4/results/raxml/RAxML_bipartitions.nuc.bs',
                'cd'= '/home/thkuo/bifo-sgi/TreePaper/Paeru.v4/results/cd/cutoff.v1/RAxML_bipartitions.conc.bs')
trs<- lapply(tr_files, function(f) read.newick(f))

source('/home/thkuo/projects/TreePaper/bin/collapse_branches.R')
cutoff<- list(br= 5e-4, bs= 60)
outgroup<- c('CH4433')
col_trs<- list()
for (x in names(trs)){
    tr<- trs[[x]]
    # reroot
    tr<- phytools::reroot(tr, node.number = which(tr$tip.label == outgroup[1]), position = 1e-6)
    trs[[x]]<- tr
    # collapse
    short_ix<- find_short_branches(tr, br_cutoff = cutoff$br)
    low_ix<- find_low_branches(tr, bs_cutoff = cutoff$bs)
    col_tr<- collapse_branches(tr= tr, target_nodes = union(short_ix, low_ix))
    col_tr<- phytools::reroot(col_tr, node.number = which(col_tr$tip.label == outgroup[1]), position = 1e-6)
    col_trs[[x]]<- col_tr
}
# map the branches
# distinguish new branches and their most recent, reproduced ancetor
node_match_mat<- matchNodes(col_trs$cd, col_trs$nuc)
newly_resolved_ix<- node_match_mat[is.na(node_match_mat[, 2]),1]
#! only include the new clades that do not overlap with each other 
#! (i.e. those not having another newly resolved nodes under it)
ix_filter<- sapply(newly_resolved_ix, 
                   function(x) sum(getDescendants(col_trs$cd, node = x) %in% newly_resolved_ix) == 0)
nonoverlap_newly_resolved_ix<- newly_resolved_ix[ix_filter]
## trace the source
# the parent of each new node
trace_resolved_origin<- function(tr, newly_resolved_nodes){
    origins<- c()
    for (ix in newly_resolved_nodes){
        # the parental node
        p_ix<- getParent(tr, ix)
        ## find the most recent ancestor that could be mapped to the nuc tree
        while (p_ix %in% newly_resolved_nodes){
            p_ix<- getParent(tr, p_ix)
        }
        origins<- c(origins, p_ix)
    }
    data.frame(node= newly_resolved_nodes, origin= origins)
}
resolved_origins<- trace_resolved_origin(col_trs$cd, newly_resolved_ix)

# compute the entropies
source('entropy_ana.functions.R')
bio_list<- as.list(bio_df)
for (n in 1:length(bio_list)){
    names(bio_list[[n]])<- rownames(bio_df)
}
traits<- bio_list$`Geographic Origin`
entropies_outputs<- node_entropies(tr= col_trs$cd, traits = traits)

# make the dataframe for ggplot
fig_df<- resolved_origins
fig_df$node<- as.character(fig_df$node)
fig_df$origin<- as.character(fig_df$origin)
fig_df$entropies<- entropies_outputs[match(fig_df$node, entropies_outputs$node_ix), ]$entropies
fig_df$origin_entropies<- entropies_outputs[match(fig_df$origin, entropies_outputs$node_ix), ]$entropies
fig_df$sizes<- entropies_outputs[match(fig_df$node, entropies_outputs$node_ix), ]$clade_size
fig_df$origin_sizes<- entropies_outputs[match(fig_df$origin, entropies_outputs$node_ix), ]$clade_size

library('ggrepel')
sub_fig_df<- fig_df[fig_df$node %in% nonoverlap_newly_resolved_ix,]
ggplot(sub_fig_df)+
    geom_abline(intercept = 0, slope = 1, color= 'grey70')+
    geom_point(aes(x= entropies, y= origin_entropies, size= sizes), alpha= 0.5, color= 'dodgerblue3')+
    geom_text_repel(aes(x= entropies, y= origin_entropies, label = node)) +
    # scale_x_continuous(expand = c(0,0))+
    # scale_y_continuous(expand = c(0,0))+
    xlab('entropy of new branches')+ylab('entropy of MRRAs')+
    # geom_point(aes(x= entropies, y= node), color= 'blue')+
    # geom_point(aes(x= origin_entropies, y= node), color= 'red')+
    theme_classic()
# ggsave('entropy.pdf', height = 6, width= 6.5)

# resolution in perfect clusters 
in_pure<- sum((sub_fig_df$origin_entropies == 0)&(sub_fig_df$entropies == 0))
# increased resolution in mixed clades
improved<- sum(sub_fig_df$entropies < sub_fig_df$origin_entropies) 
# not improved cases
nonimproved<- sum(sub_fig_df$entropies >= sub_fig_df$origin_entropies) - in_pure
