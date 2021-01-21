
extract_all_bipart_tips_unrooted<-function(target_tr){
    # compute all bipartitions
    library(phytools)
    library(ape)
    bipart_strs<- c()
    for (x in target_tr$edge[,2]){
        if (x <= Ntip(target_tr)){
            bipart_strs<- c(bipart_strs, NA)
        }else{
            s1<- ape::extract.clade(phy = target_tr, node= x)
            if (Ntip(s1) < Ntip(target_tr)){
                s2<- ape::drop.tip(phy= target_tr, tip= s1$tip.label)
                tips1<- s1$tip.label
                tips1<- tips1[order(tips1)]
                tips1_str<- paste0(tips1, collapse = ':')
                tips2<- s2$tip.label
                tips2<- tips2[order(tips2)]
                tips2_str<- paste0(tips2, collapse = ':')
                tips_strs<- c(tips1_str, tips2_str)
                tips_strs<- tips_strs[order(tips_strs)]
                bipart_str<- paste0(tips_strs, collapse = '::')
                bipart_strs<- c(bipart_strs, bipart_str)
            }   
        }
    }
    return(bipart_strs)
}

matchBranches_unrooted<-function(target_tr, ref){
    # match the branches for a pair of trees
    library(phytools)
    clades<- lapply(list(target= target_tr, ref= ref), function(tr) extract_all_bipart_tips_unrooted(tr))
    correctness<- match(clades$target, clades$ref)
    correctness[is.na(clades$target)]<- NA
    match_mat<- cbind(target_tr$edge[, 2], sapply(correctness, function(x) ifelse(is.na(x), NA, ref$edge[x, 2])))
    match_mat[match_mat[, 1] > Ntip(target_tr),]
}

read.ll_mat<- function(ll_mat_f, all_rear_trs_f, cpus_for_reading_table= 1){
    # read and calculate the PSLLs for each branch
    ll_mat<- data.frame()
    if (cpus_for_reading_table>1){
        library(data.table)
        ll_mat_data.table<- fread(ll_mat_f, header= F, skip= 1, 
				  sep= ' ', nThread= cpus_for_reading_table, showProgress= T)
        ll_mat<- as.data.frame(ll_mat_data.table)
        rownames(ll_mat)<- sapply(ll_mat[, 1], function(x) unlist(strsplit(x = x,split = '\t' ))[1])
        ll_mat[, 1]<- sapply(ll_mat[, 1], function(x) as.numeric(unlist(strsplit(x = x,split = '\t' ))[2]))
        colnames(ll_mat)<- as.character(1:ncol(ll_mat))
    }else{
        ll_mat<- read.table(ll_mat_f, header= F, skip = 1)
        rownames(ll_mat)<- ll_mat[, 1]
        ll_mat<- ll_mat[, 2:ncol(ll_mat)]
        colnames(ll_mat)<- as.character(1:ncol(ll_mat))
    }
    dim(ll_mat)
    #' higher the value, less supporting to the branch composing the current nucleotide tree
    reference_ll<- ll_mat[1, ]
    ll_gap_mat<- t(sapply(1:nrow(ll_mat), 
                          function(x) ll_mat[x, ]-reference_ll))
    rownames(ll_gap_mat)<- rownames(ll_mat)
    ll_gap_mat<- ll_gap_mat[2:nrow(ll_gap_mat),]
    
    #' Which branch was rearranged?
    library(phytools)
#    source('./tree_correctness_unrooted.R')
    # 
    all_rear_trs<- read.newick(all_rear_trs_f)
    rear_br_ix_list<- list()
    for (x in 1:length(all_rear_trs)){
        print(x)
        br_match_mat<- matchBranches_unrooted(all_rear_trs[[1]], all_rear_trs[[x]])
        rear_br_ix<- br_match_mat[which(is.na(br_match_mat[, 2])), 1]
        rear_br_ix_list[[x]]<- rear_br_ix
    }
    names(rear_br_ix_list)<- sprintf('tr%d', 1:length(all_rear_trs))
    
    #' merge tress that were made by reaaranging the same branch
    merged_ll_gap_mat<- rep(0, ncol(ll_gap_mat))
    for (node_ix in unique(unlist(rear_br_ix_list))){
        # node_ix<- unique(unlist(rear_br_ix_list))[1]
        tr_ixs<- names(rear_br_ix_list)[which(rear_br_ix_list == node_ix)]
        merged_ll_gap_mat<- rbind(merged_ll_gap_mat, 
                                  apply(ll_gap_mat[tr_ixs, ], 2, function(c) mean(unlist(c))))
    }
    merged_ll_gap_mat<- merged_ll_gap_mat[2:nrow(merged_ll_gap_mat), ]
    rownames(merged_ll_gap_mat)<- as.character(unique(unlist(rear_br_ix_list)))
    merged_ll_gap_mat<- merged_ll_gap_mat[order(rowSums(merged_ll_gap_mat)), ]
    
    list('branch.merged'= merged_ll_gap_mat, 'reference_ll'= reference_ll, 
         'ori'= ll_gap_mat, 'rear.trs'= all_rear_trs)
}

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


nuc_psll_f<-snakemake@input[['nuc_pslls']]
rear_trs_f<-snakemake@input[['rear_trs']]
psll_sums_out_f<- snakemake@output[['nuc_psll_sum_per_br']]
col_cutoff_pr<- snakemake@params[['col_cutoff_pr']]
col_tr_out_f<- snakemake@output[['col_nuc_tr']]

#' compute the sums 
nuc_psll<- read.ll_mat(ll_mat_f = nuc_psll_f, 
	       all_rear_trs_f = rear_trs_f, 
	       cpus_for_reading_table= snakemake@threads[1])
nuc_psll_sum_per_branch<- rowSums(nuc_psll$branch.merged)
write.table(as.data.frame(nuc_psll_sum_per_branch), 
	    psll_sums_out_f, sep= '\t')

#' determine the branches to collapse
nuc_psll_cutoff<- quantile(nuc_psll_sum_per_branch, probs = col_cutoff_pr)
br.ix_to_col<- names(nuc_psll_sum_per_branch[which(nuc_psll_sum_per_branch >= nuc_psll_cutoff)])
br.ix_to_col<- as.numeric(br.ix_to_col)

if (! dir.exists(dirname(col_tr_out_f))){
    dir.create(dirname(col_tr_out_f))
}
nuc_tr<- nuc_psll$rear.trs[[1]]
write.tree(
    phy= unroot(collapse_branches(tr = nuc_tr, target_nodes = br.ix_to_col)),
    file = col_tr_out_f)

