
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
source('/net/metagenomics/data/from_moni/old.tzuhao/TreePaper/Paeru.v4/bin.3/per-site_LogLikelihood/visualize/pslltools/read_ll_mat.R')
nuc_psll_f<-snakemake@input[['nuc_pslls']]
rear_trs_f<-snakemake@input[['rear_trs']]
col_tr_out_f<- snakemake@output[['col_nuc_tr']]
cpu_num<- snakemake@threads
cutoff_perc<- as.numeric(snakemake@params[['ll_perc']] )

#' compute the sums 
nuc_psll<- read.ll_mat(ll_mat_f = nuc_psll_f, 
	       all_rear_trs_f = rear_trs_f, 
	       large_table= T, cpu_num=cpu_num)
nuc_psll_sum_per_branch<- rowSums(nuc_psll$branch.merged)

#' determine the branches to collapse
nuc_psll_cutoff<- quantile(nuc_psll_sum_per_branch, probs = cutoff_perc)
br.ix_to_col<- names(nuc_psll_sum_per_branch[which(nuc_psll_sum_per_branch >= nuc_psll_cutoff)])
br.ix_to_col<- as.numeric(br.ix_to_col)

if (! dir.exists(dirname(col_tr_out_f))){
    dir.create(dirname(col_tr_out_f))
}
nuc_tr<- nuc_psll$rear.trs[[1]]
write.tree(
    phy= unroot(collapse_branches(tr = nuc_tr, target_nodes = br.ix_to_col)),
    file = col_tr_out_f)

