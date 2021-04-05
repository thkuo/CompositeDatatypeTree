

extract_all_clades_tips<-function(tr){
    library(ape)
    tr_all_clades<- c()
    for (n in ((Ntip(tr)+1):(Ntip(tr)+Nnode(tr)))){
        sub<- extract.clade(tr, node = n)
        tips<- sub$tip.label
        tips<- tips[order(tips)]
        tips_str<- paste0(tips, collapse = ':')
        tr_all_clades<- c(tr_all_clades, tips_str)
    }
    return(tr_all_clades)
}

tree_correctness<- function(target_tr, ref){
    clades<- lapply(list(target= target_tr, ref= ref), function(tr) extract_all_clades_tips(tr))
    all_clades<- union(clades$target, clades$ref)
    clades_pa<- sapply(clades, function(c) ifelse(all_clades %in% c, 1, 0))
    
    library(limma)
    pa_counts<- vennCounts(clades_pa)
    correctness<- pa_counts[, 'Counts']
    names(correctness)<- c('TN', 'FN', 'FP', 'TP')
    correctness['recall']<- correctness['TP']/(correctness['TP']+correctness['FN'])
    correctness['precision']<- correctness['TP']/(correctness['TP']+correctness['FP'])
    correctness['f1']<- 2*(correctness['precision']*correctness['recall'])/(correctness['precision']+correctness['recall'])
    return(correctness)
}

make_three_venn_diagram<- function(a, b, c, names= NA){
    all_elements<- Reduce(union, list(a,b,c)) 
    bin_mat<- sapply(list(a,b,c), function(x) ifelse(all_elements %in% x, 1, 0))
    library(limma)
    v_counts<- vennCounts(bin_mat)
    cir_colors<- c('#023059', '#858C4D', '#BF6B04')
    if (any(is.na(names))){
        vennDiagram(v_counts,circle.col=cir_colors)
    }else{
        vennDiagram(v_counts,circle.col=cir_colors, 
                    names = names)
    }
    
}

tips_str<- function(tr, n){
    tips<- extract.clade(tr, node= n)$tip.label
    paste0(tips[order(tips)], collapse = ':')    
}

three_trees_venn<- function(nuc_tr, cd_tr, ref_tr){
    trs_to_compare<- list(nuc= nuc_tr, cd= cd_tr, ref= ref_tr)
    tips_strs<- lapply(trs_to_compare, function(tr)
        sapply((Ntip(tr)+1):(Ntip(tr)+Nnode(tr)) , function(n) tips_str(tr, n)))
    make_three_venn_disgram(a= tips_strs$nuc, b= tips_strs$cd, c= tips_strs$ref, 
                            names = c('nucleotide', 'composite\ndatatype', 'reference'))
}
