#' Parse within-patient trees from the ctree
#' Not extracting them before ttree and wtrees are merged
#' because the node indices were not yet determined, causing difficultied in 
#' deterning the connection between wtrees
t<- ctree$ctree
patients_ics<- unique(t[, 4])
patients_ics<- patients_ics[order(patients_ics)]

for (patient_ix in patients_ics){
    sub_ctree<- ctree$ctree[which(ctree$ctree[,4] == patient_ix), ]
    node_ics<- which(ctree$ctree[,4] == patient_ix)
    # find the infected time
    ori_time<- 0
    ori<- which(ctree$ctree[,2] == node_ics[which(sub_ctree[,1] == min(sub_ctree[,1]))])
    if (length(ori) == 0){
        ori<- which(ctree$ctree[,3] == sub_ctree[which(sub_ctree[,1] == min(sub_ctree[,1])), 5])
    }
    ori_time<- ctree$ctree[ori, 1]
    root.edge.length<-sub_ctree[which(sub_ctree[,1] == min(sub_ctree[,1])), 1]-ori_time
    
}
