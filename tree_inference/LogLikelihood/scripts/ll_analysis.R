
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

determine_rearranged_br<-function(reference_tr,tr){
  br_match_mat<- matchBranches_unrooted(reference_tr, tr)
  rear_br_ix<- br_match_mat[which(is.na(br_match_mat[, 2])), 1]
  return(rear_br_ix)
}

parse.ll<- function(ll_f, all_rear_trs_f, cpus_for_reading_table= 1){
  library(stringr)
  if (!require('future.apply')){
    install.packages('future.apply', repos= 'https://cloud.r-project.org')
  }
  con<- file(ll_f, 'r')
  lines<- readLines(con)
  close(con)
  pat<- 'Tree ([0-9]+) Likelihood (-[0-9\\.]+)'
  tr.ixs<- c()
  lls<- c()
  for (l in lines){
    if (grepl(pat, l)){
      re_out<- str_match( string= l,pattern= pat)
      print(re_out)
      tr.ix<- as.integer(re_out[2])
      ll<- as.numeric(re_out[3])
      tr.ixs<- c(tr.ixs, tr.ix)
      lls<- c(lls, ll)
    }
  }
  print(length(lls))
  ll_gap<- lls-lls[1]
  stopifnot(length(ll_gap) == length(lls))
  print(length(ll_gap))
  print(length(sprintf('tr%d',  1:length(ll_gap))))
  names(ll_gap)<-sprintf('tr%d',  1:length(ll_gap))

  #' Which branch was rearranged?
  library(phytools)
  library(future.apply)
  plan(multiprocess, workers = cpus_for_reading_table)
  all_rear_trs<- read.newick(all_rear_trs_f)
  rear_br_ix_list<- list()
  rear_br_ix_list<- future_lapply(all_rear_trs, function(tr) determine_rearranged_br(all_rear_trs[[1]], tr))
  names(rear_br_ix_list)<- sprintf('tr%d', 1:length(all_rear_trs))

  merged_ll_gap<-list()
  for (node_ix in unique(unlist(rear_br_ix_list))){
    # node_ix<- unique(unlist(rear_br_ix_list))[1]
    tr_ixs<- names(rear_br_ix_list)[which(rear_br_ix_list == node_ix)]
    merged_ll_gap[[as.character(node_ix)]]<- mean(ll_gap[tr_ixs])
  }
    
  list('branch.merged'= merged_ll_gap, 'reference_ll'= ll[1], 
       'rear.trs'= all_rear_trs)
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


nuc_lls_f<- snakemake@input[['nuc_lls']]
rear_trs_f<-snakemake@input[['rear_trs']]
ll_sums_out_f<- snakemake@output[['nuc_ll_sum_per_br']]
col_tr_out_f<- snakemake@output[['col_nuc_tr']]
cutoff_perc<- as.numeric(snakemake@params[['cutoff_perc']])

#' compute the sums 
nuc_ll<- parse.ll(ll_f = nuc_lls_f, 
	       all_rear_trs_f = rear_trs_f, 
	       cpus_for_reading_table= snakemake@threads[1])
write.table(t(data.frame(nuc_ll$branch.merged,check.names = F )), 
	    ll_sums_out_f, sep= '\t')
nuc_ll_sum_per_branch<- unlist(nuc_ll$branch.merged)
print(nuc_ll_sum_per_branch)

#' determine the branches to collapse
nuc_ll_cutoff<- quantile(nuc_ll_sum_per_branch, probs = cutoff_perc)
br.ix_to_col<- names(nuc_ll_sum_per_branch[which(nuc_ll_sum_per_branch >= nuc_ll_cutoff)])
br.ix_to_col<- as.numeric(br.ix_to_col)

if (! dir.exists(dirname(col_tr_out_f))){
    dir.create(dirname(col_tr_out_f))
}
nuc_tr<- nuc_ll$rear.trs[[1]]
write.tree(
    phy= unroot(collapse_branches(tr = nuc_tr, target_nodes = br.ix_to_col)),
    file = col_tr_out_f)

