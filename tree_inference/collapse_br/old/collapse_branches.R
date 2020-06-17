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
