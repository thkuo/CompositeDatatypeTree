#' Parse within-patient trees from the ctree
#' Not extracting them before ttree and wtrees are merged
#' because the node indices were not yet determined, causing difficultied in 
#' deterning the connection between wtrees

wtreesFromCTree<- function(ctree){
    wtrees<-c()
    patients_ics<- unique(ctree$ctree[, 4])
    patients_ics<- patients_ics[order(patients_ics)]
    patients_ics<- patients_ics[patients_ics> 0]
    
    for (patient_ix in patients_ics){
        print(patient_ix)
        sub_ctree<- ctree$ctree[which(ctree$ctree[,4] == patient_ix), ]
        if (sum(ctree$ctree[,4] == patient_ix)==1){
            sub_ctree<- t(as.data.frame(sub_ctree))
        }
        node_ics<- which(ctree$ctree[,4] == patient_ix)
        # find the infected time
        ori_time<- 0
        ori<- which(ctree$ctree[,2] == node_ics[which(sub_ctree[,1] == min(sub_ctree[,1]))])
        if (length(ori) == 0){
            ori<- which(ctree$ctree[,3] == sub_ctree[which(sub_ctree[,1] == min(sub_ctree[,1])), 5])
        }
        ori_time<- ctree$ctree[ori, 1]
        root.edge.length<-sub_ctree[which(sub_ctree[,1] == min(sub_ctree[,1])), 1]-ori_time
        sub_ctree<-rbind(sub_ctree,  ctree$ctree[ori, ])
        for (n in 2:4){
            sub_ctree[, n]<- as.integer(sub_ctree[, n])
        }
        
        sub_ctree<- cbind(sub_ctree, c(node_ics, ori))
        ## reformate the wtree
        ## move the tips to the top
        is_tips<- (!sub_ctree[,2] %in% sub_ctree[,5]) & (!sub_ctree[,3] %in% sub_ctree[,5])
        sub_ctree<- sub_ctree[order(!is_tips, -sub_ctree[,1]),]
        wtrees<- c(wtrees, list(sub_ctree))
    }
    wtrees
}
