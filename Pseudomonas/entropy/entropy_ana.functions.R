
compute_br_cutoff<- function(nuc_tr){
    br_l<- nuc_tr$edge.length
    br_cutoff<- as.numeric(quantile(br_l, probs= 0.5)-(quantile(br_l, probs= 0.5)-quantile(br_l, probs= 0.25))*1.5)
    return(br_cutoff)
}

compute_bs_cutoff<- function(nuc_tr){
    return(30)
}

find_root_and_first_children<- function(tr){
    # find the index of root node
    tr_info<- fortify(tr)
    root_ix<- as.character(tr_info[tr_info$parent == tr_info$node, ]$node)
    # the first level of internal nodes
    first_child_ixs<- as.character(
        tr_info[(tr_info$parent == root_ix) & (!tr_info$isTip) & (tr_info$node != root_ix), ]$node)
    return(list(root= root_ix, first_children= first_child_ixs))
}


node_entropies_multiple_trees<- function(trs, traits){
    # for a list of trees
    # calculate the entropy values of nodes in each tree
    entropies_df<- data.frame()
    if (length(unique(traits)) < 2){
        for (tr_name in names(trs)){
            entropies<- rep(0, Nnode(trs[[tr_name]]))
            datatype<- rep(tr_name, length(entropies))
            node_ix<- (Ntip(trs[[tr_name]])+1):(Nnode(trs[[tr_name]])+Ntip(trs[[tr_name]]))
            tmp_entropies_df<- data.frame(entropies, datatype, node_ix)
            if (nrow(entropies_df) == 0){
                entropies_df<- tmp_entropies_df
            }else{
                entropies_df<- rbind(entropies_df, tmp_entropies_df)
            }
        }
    }else{
        print(traits)
        recons_outputs<- lapply(trs, 
                                function (tr) rerootingMethod(tree = tr,x = traits, model = 'ER'))
        print(dim(recons_outputs))
        for (tr_name in names(recons_outputs)){
            ma<- recons_outputs[[tr_name]]$marginal.anc
            entropies<- apply(ma, 1, function(r) entropy.empirical(r, unit= 'log2'))
            datatype<- rep(tr_name, length(entropies))
            node_ix<- names(entropies)
            entropies<- as.numeric(entropies)
            tmp_entropies_df<- data.frame(entropies, datatype, node_ix)
            if (nrow(entropies_df) == 0){
                entropies_df<- tmp_entropies_df
            }else{
                entropies_df<- rbind(entropies_df, tmp_entropies_df)
            }
        }   
    }
    return(entropies_df)
}


node_entropies<- function(tr, traits, node_ix= NULL){
    # calculate the entropy values based only on tips
    # freee of ancestral reconstruction
    # for one tree each call
    library(ggtree)
    library(entropy)
    entropies_df<- data.frame()
    tr_info<- fortify(tr)
    tr_info$entropy<- 0
    
    if (is.null(node_ix)){
        node_ix<- tr_info$node    
    }
    incl_tips<- sapply(node_ix, function(n)get.offspring.tip(tr, node = n))
    incl_traits<- lapply(incl_tips, function(t) traits[t])
    entropies<- sapply(incl_traits, function(t) entropy.empirical(table(t)/length(t)))
    clade_size<- sapply(incl_traits, function(t) length(t))
    entropies_df<- data.frame(entropies, node_ix,clade_size)

    return(entropies_df)
}


root_and_children_entropy_change<- function(trs, entropies_outputs, tr_name, trait_name){
    entropies<- entropies_outputs[[trait_name]][entropies_outputs[[trait_name]]$datatype == tr_name, ]$entropies
    names(entropies)<- entropies_outputs[[trait_name]][entropies_outputs[[trait_name]]$datatype == tr_name, ]$node_ix
    tr<- trs[[tr_name]]
    root_and_children<- find_root_and_first_children(tr)
    sum(entropies[root_and_children$first_children] - entropies[root_and_children$root])/
        length(root_and_children)
}
