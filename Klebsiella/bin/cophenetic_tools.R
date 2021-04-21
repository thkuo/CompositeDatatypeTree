library(reshape2)
library(viridis)
library(lemon)
library(ggplot2)
norm_tr<- function(tr){
    new_br<- tr$edge.length
    new_br<- new_br/sum(new_br)
    new_tr<- compute.brlen(tr, new_br)
    return(new_tr)
}
mat2df<- function(m, traits){
    df<- melt(m, stringAsFactor= F)
    df[, 1]<- as.character(df[, 1])
    df[, 2]<- as.character(df[, 2])
    df<- df[df[, 1] > df[, 2], ]
    df<- df[order(traits[df[, 1]],traits[df[, 2]]), ]
    for (n in 1:nrow(df)){
        sample_pair<- c(df[n, 1], df[n, 2])
        traits_pair<- c(traits[df[n, 1]], traits[df[n, 2]])
        sample_pair<- sample_pair[order(traits_pair)]
        df[n, 1]<-sample_pair[1]
        df[n, 2]<-sample_pair[2]
    }
    df
}
conc_locs<- function(locs){
    ordered_locs<- locs[order(locs)]
    paste0(ordered_locs, collapse = '-')
}
compute_cophe_df<- function(trs, traits){
    if (length(trs) > 2){
        print('Only the first two trees are considered')
        trs<- trs[[1:2]]
    }
    #' un-weigthed
    cophe_mats<- lapply(trs, function(t)cophenetic.phylo(compute.brlen(t, 1)))
    #' #' weighted
    #' norm_trs<-lapply(trs, function(t) norm_tr(t))
    # cophe_mats<- lapply(norm_trs, function(t)cophenetic.phylo(t))
    cophe_dfs<- lapply(cophe_mats, function(m) mat2df(m, traits))
    cophe_one_df<- merge(cophe_dfs[[1]], cophe_dfs[[2]], by=c('Var1', 'Var2')) 
    colnames(cophe_one_df)<- c('sample1', 'sample2', 'dist1', 'dist2')
    cophe_one_df$loc1<-  as.character(traits[cophe_one_df$sample1])
    cophe_one_df$loc2<-  as.character(traits[cophe_one_df$sample2])
    cophe_one_df$loc_pair<- apply(cophe_one_df, 1, function(r) conc_locs(c(r['loc1'], r['loc2'])))
    cophe_one_df
}

determine_sig<- function(count_df, h1, h2, alternative= NA){
    print(sprintf('%s-%s', h1, h2))
    if (is.na(alternative)){
        alternative<-'less'
        if (h1 == h2){
            alternative= 'greater'
        }    
    }
    
    sig<- 1
    t_out<- t.test(count_df[(count_df$loc1 == h1) & (count_df$loc2 == h2), 1], 
                   count_df[(count_df$loc1 == h1) & (count_df$loc2 == h2), 2], 
                   paired = TRUE, alternative = alternative)      
    sig<- t_out$p.value
    sig
}

sig_for_all_loc_pairs<- function(cophe_one_df){
    all_loc_pairs<- unique(cophe_one_df[, c('loc1', 'loc2')])
    print(all_loc_pairs)
    loc1s<- c()
    loc2s<- c()
    sigs<- c()
    alternative_hs<- c()
    pairs_num<- c()
    for (n in 1:nrow(all_loc_pairs)){
        h1<- all_loc_pairs[n, 'loc1']
        h2<- all_loc_pairs[n, 'loc2']
        print(h1)
        print(h2)
        loc_filter<- (cophe_one_df$loc1 == h1) & (cophe_one_df$loc2 == h2)
        if (sum(loc_filter) == 0){
            next
        }
        count_df<- cophe_one_df[loc_filter, 
                                c('dist1', 'dist2', 'loc1', 'loc2')]
        
        # print(table(count_df$loc1, count_df$loc2))
        loc1s<- c(loc1s, h1)
        loc2s<- c(loc2s, h2)
        sigs<- c(sigs, determine_sig(count_df, h1, h2))
        print(sigs[length(sigs)])
        pairs_num<- c(pairs_num, nrow(count_df))
        alternative_hs<-c(alternative_hs, ifelse(h1==h2, 'greater', 'less'))
    }
    data.frame(loc1s, loc2s, sigs, alternative_hs, pairs_num)   
}
show_violin<- function(cophe_one_df){
    plot_max_val<- max(abs(cophe_one_df$dist2-cophe_one_df$dist1))
    ggplot(cophe_one_df, aes(x= loc_pair))+
        geom_violin(aes(y= dist2-dist1, fill= loc1==loc2))+
        geom_hline(yintercept = 0, color= 'grey90')+
        scale_y_continuous(limits= c(-plot_max_val, plot_max_val), 
                           breaks= seq(from= -plot_max_val,to= plot_max_val,by=1 ))+
        xlab('locational pairs')+
        theme_classic()
}
show_density<- function(cophe_one_df, dist1_name= '', dist2_name= ''){
    binning_breaks<- unique(c(cophe_one_df$dist1, cophe_one_df$dist2))
    binning_breaks<-binning_breaks[order(binning_breaks)]
    ggplot(cophe_one_df)+
        # geom_point(aes(x= nuc_dist, y= cd_dist), alpha= .3, color= 'orange')+
        geom_bin2d(aes(x= dist1, y= dist2), breaks=binning_breaks) +
        scale_fill_viridis_c(direction = -1)+
        # scale_fill_gradient(low ='#BAE4BC', high = 'darkorange2', name= 'sample pairs')+
        geom_abline(intercept = 0, slope = 1, color= 'grey50')+
        xlab(dist1_name)+
        ylab(dist2_name)+
        coord_fixed(ratio= 1)+
        theme_classic()
}
show_pair_counts<- function(cophe_one_df, loc_pairs){
    sub_cophe_one_stat_df<- dcast(
        as.data.frame(
            table(ifelse(cophe_one_df$dist2 > cophe_one_df$dist1, 'increased', 'decreased'),
                  cophe_one_df$loc_pair)), Var2~Var1)
    colnames(sub_cophe_one_stat_df)<-c('hospital pairs', 'decreased', 'increased')
    rownames(sub_cophe_one_stat_df)<- sub_cophe_one_stat_df$`hospital pairs`
    sub_cophe_one_stat_df$`hospital pairs`<- NULL
    sub_cophe_one_perc_df<- sub_cophe_one_stat_df/rowSums(sub_cophe_one_stat_df) * 100
    colnames(sub_cophe_one_perc_df)<- c( 'decreased', 'increased')
    
    sub_cophe_one_stat_df$locs<- rownames(sub_cophe_one_stat_df)
    sub_cophe_one_perc_df$locs<- rownames(sub_cophe_one_perc_df)
    sub_cophe_one_stat_long<- melt(sub_cophe_one_stat_df)
    sub_cophe_one_stat_long<- sub_cophe_one_stat_long[order(sub_cophe_one_stat_long$locs, sub_cophe_one_stat_long$variable), ]
    sub_cophe_one_perc_long<- melt(sub_cophe_one_perc_df)
    sub_cophe_one_perc_long<- sub_cophe_one_perc_long[order(sub_cophe_one_perc_long$locs, sub_cophe_one_perc_long$variable), ]
    
    plot_df<- cbind(sub_cophe_one_stat_long, sub_cophe_one_perc_long$value)
    plot_df$locs<- factor(plot_df$locs, levels= loc_pairs)
    
    colnames(plot_df)<- c('locs', 'change', 'count', 'percentage')
    ggplot(plot_df, aes(x= locs, y= percentage,fill= change))+
        geom_bar(stat="identity")+
        geom_text(aes(label=count, color= change), position = position_stack(vjust= .5))+
        scale_fill_manual(values=c('decreased'='yellowgreen', 'increased'='dodgerblue4'))+
        scale_color_manual(values= c('decreased'= 'black', 'increased'= 'grey'))+
        xlab('locational pairs')+
        theme_classic()
}
