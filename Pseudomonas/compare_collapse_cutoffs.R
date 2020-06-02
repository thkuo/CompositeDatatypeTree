## Visualize the cutoffs in the br-bs plot
library(phytools)
# First, reproduce the plot with nucleotide tree 
nuc_tr_f<- '/home/thkuo/bifo-sgi/TreePaper/Paeru.v4/results/raxml/RAxML_bipartitions.nuc.bs'
nuc_tr<- read.newick(nuc_tr_f)

outgroup<- c('CH4433')
nuc_tr<- phytools::reroot(tr, node.number = which(tr$tip.label == outgroup[1]), position = 1e-6)

# Read the cutoffs
read_cutoff_files<- function(f){
    cutoff_df<- read.table(f, sep= '\t', header = F)
    cutoff<- cutoff_df[, 2]
    names(cutoff)<- cutoff_df[, 1]
    as.list(cutoff)
}
cutoff_files<- sapply(1:4, function(n) sprintf('/home/thkuo/bifo-sgi/TreePaper/Paeru.v4/results/cd/cutoff.v%d/cutoffs.txt',n))
names(cutoff_files)<- sapply(1:4, function(n) sprintf('v%d', n))
cutoff_set<- lapply(cutoff_files, function(f) read_cutoff_files(f))
brs<- c()
bss<- c()
for (x in names(cutoff_set)){
    brs<- c(brs, cutoff_set[[x]][['br']])
    bss<- c(bss, cutoff_set[[x]][['bs']])
}
cutoff_set_df<- data.frame(list(br= brs, bs= bss), row.names = names(cutoff_set))

# the br-bs plot
library(ggplot2)
library(viridis)
library(ggtree)
## all nodes
tr_info<- fortify(nuc_tr)
plot_df<- tr_info[!tr_info$isTip & tr_info$label != 'Root', c('branch.length', 'label', 'isTip')]
plot_df$clade_size<- apply(tr_info[!tr_info$isTip & tr_info$label != 'Root', ], 
                           1, function(r) length(get.offspring.tip(tr, as.numeric(r['node']))))
plot_df$label<- as.numeric(plot_df$label)
plot_df<- plot_df[! is.na(plot_df$label), ]
sc_p<- ggplot(plot_df)+
    geom_point(aes(x= branch.length, y= label, color= clade_size), size = 2, alpha= .3)+
    geom_point(data= cutoff_set_df, color= '#F2293A', size= 2, shape= 4, inherit.aes = F,
               aes(x= br, y= bs))+
    scale_x_continuous(trans= 'log10', na.value = 0)+
    # scale_color_brewer(type = 'seq', palette = 'YlGnBu')+
    scale_color_viridis(option = 'A')+
    # scale_color_gradient()+
    xlab('branch length')+ylab('bootstrap support')+
    theme_classic()

hist_right<- ggplot(plot_df)+
    geom_histogram(aes(x= label,y=cumsum(..count..)), bins = 40, boundary= 0)+
    xlab('')+ylab('')+
    coord_flip()+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
hist_top<- ggplot(plot_df)+
    geom_histogram(aes(x= branch.length, y=cumsum(..count..)), bins = 40, boundary= 0)+
    scale_x_continuous(trans= 'log10', na.value = 0)+
    xlab('')+ylab('')+
    theme_classic()
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
    theme(axis.ticks=element_blank(),
          panel.background=element_blank(),
          axis.text.x=element_blank(), axis.text.y=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank())
library(ggpubr)
p<- ggarrange(hist_top, empty, sc_p, hist_right,
          ncol=2, nrow=2, widths=c(5, 1), heights=c(1, 5), common.legend = T, legend = 'bottom')
# ggsave(p, filename = 'cutoffs_br-bs.pdf')
# compare the improvement of resolution under different cutoff values
## For each, load the trees
source('/home/thkuo/projects/TreePaper/bin/collapse_branches.R')
library(phytools)
col_trs<- list()
result_dir<- '/home/thkuo/bifo-sgi/TreePaper/Paeru.v4/results/cd/'
col_trs<- list()
for (n in 1:4){
    set_name<- sprintf('v%d', n)
    print(set_name)
    cd_tr_f<- sprintf('%s/cutoff.%s/RAxML_bipartitions.conc.bs', result_dir, set_name)
    cd_tr<- read.newick(cd_tr_f)
    trs<- list(nuc= nuc_tr, cd= cd_tr)
    col_tr_pair<- list()
    cutoff<- cutoff_set[[set_name]]
    for (x in names(trs)){
        tr<- trs[[x]]
        # collapse
        short_ix<- find_short_branches(tr, br_cutoff = cutoff$br)
        low_ix<- find_low_branches(tr, bs_cutoff = cutoff$bs)
        col_tr<- collapse_branches(tr= tr, target_nodes = union(short_ix, low_ix))
        col_tr<- phytools::reroot(col_tr, node.number = which(col_tr$tip.label %in% outgroup), position = 1e-6)
        col_tr_pair[[x]]<- col_tr
    }
    
    col_trs[[set_name]]<- col_tr_pair
}
# for (set_name in names(col_trs)){
#     rerooted_trs<- lapply(col_trs[[set_name]], 
#                           function(tr) root(phy=tr, outgroup = load_outgroup(), resolve.root = T))
#     col_trs[[set_name]]<-rerooted_trs
# }
# check if they are collapsed
# tr_files<- list(
#     nuc= '/home/thkuo/bifo-sgi/TreePaper/Paeru.vPre3-3/results/after_nuc_tree/v4/nuc.col.nwk',
#     md= '/home/thkuo/bifo-sgi/TreePaper/Paeru.vPre3-3/results/after_nuc_tree/v4/md.col.nwk')
# trs<- lapply(tr_files, function(f) phytools::read.newick(f))
# lapply(trs, function(tr) count_unresolved(tr)$summary[1, 1])

source('/home/thkuo/projects/TreePaper/bin/count_unresolved.v3.R')
source('/home/thkuo/projects/TreePaper/bin/dRF.R')
#sapply(col_trs, function(tr_set) lapply(tr_set, function(tr) is.binary(tr)))
# the resolution levels (number of resolved nodes)
# resolved_counts_df<- t(sapply(col_trs, function(tr_set) lapply(tr_set, function(tr) count_unresolved(tr)$summary[1, 1])))
count_multifurcating_nodes<- function(tr){
    is_multi<- table(tr$edge[,1]) > 2
    is_multi<- as.integer(names(is_multi[is_multi]))
    length(is_multi)
}
count_binary_nodes<- function(tr){
    Nnode(tr)-count_multifurcating_nodes(tr)
}

total_nodes<- t(sapply(col_trs, function(tr_set) lapply(tr_set, function(tr) Nnode(tr))))
bi_nodes<- t(sapply(col_trs, function(tr_set) lapply(tr_set, function(tr) count_binary_nodes(tr))))
# t(sapply(col_trs, function(tr_set) lapply(tr_set, function(tr) count_binary_nodes(tr)/Nnode(tr))))
multi_nodes<- t(sapply(col_trs, function(tr_set) lapply(tr_set, function(tr) count_multifurcating_nodes(tr))))
# cbind(cutoff_set_df, resolved_counts_df)
# apply(resolved_counts_df, 1, function(r) r[['md']]/r[['nuc']])
