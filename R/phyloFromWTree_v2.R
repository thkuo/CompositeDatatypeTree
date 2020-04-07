#' Converts a phylogenetic tree into an ape phylo object
#' revised from TransPhylo::phyloFromPtree
#' @param wtree_from_ctree phylogenetic tree; col1: time; col2-3: indices of offsprings; col4: patien index; col5: node/tip indices
#' @return phylo object
#' @export
phyloFromWTree <- function(wtree_from_ctree) {
    wtree<- wtree_from_ctree
    # nam=wtree$nam
    # wtree=wtree$wtree
    # 
    # n<-ceiling(nrow(wtree)/2)
    # internal_ics<- which((wtree[1:(nrow(wtree)-1),2]!=0 | wtree[1:(nrow(wtree)-1),3]!=0))
    # tip_ics<- setdiff(wtree[1:(nrow(wtree)-1),5], internal_ics)
    
    tip_ics<-which((!wtree[,2] %in% wtree[,5]) & (!wtree[,3] %in% wtree[,5]))
    internal_ics<- setdiff(1:(nrow(wtree)-1), tip_ics)
    nam<- as.character(wtree[tip_ics, 5])
    node.nam<- as.character(wtree[internal_ics, 5])
    n<- length(internal_ics)
    
    root<-which(wtree[1:(nrow(wtree)-1),1]==min(wtree[1:(nrow(wtree)-1),1]))
    root.edge.length<- wtree[root,1] - min(wtree[nrow(wtree),1])
    
    # if (n==1) return(ape::read.tree(text='(1);'))
    if (length(tip_ics)==1) {
        ori<- wtree[nrow(wtree), 5]
        ori.time<- wtree[nrow(wtree), 1]
        root.edge.length<- wtree[root,1]- ori.time
        tr<- ape::read.tree(text=sprintf('(%s:0.0,pseudo:0.0);', wtree[1, 5] ))
        tr$root.edge<- root.edge.length
        return(list(phylo= tr, ori= ori, ori.time= ori.time))
    }
    
    tr<-list()
    # tr$Nnode<-n-1
    tr$Nnode<-n
    tr$abel<-nam
    # tr$edge<-matrix(0,nrow(wtree)-1,2)
    # except for the root
    tr$edge<-matrix(0,length(c(tip_ics, setdiff(internal_ics, root))),2)
    # tr$edge.length<-rep(0,nrow(wtree)-1)
    tr$edge.length<-rep(0,nrow(tr$edge))
    # tr$edge<-matrix(0,n*2-2,2)
    # tr$edge.length<-rep(0,n*2-2)
    
    iedge<-1
    # root<-which(wtree[,1]==min(wtree[,1]))
    tra<-c(1:length(tip_ics),root,setdiff(internal_ics,root))
    tra2<-1:length(tra)
    tra[tra]<-tra2
    # for (i in (n+1):(2*n-1)) {
    # for (i in internal_ics){
    for (i in internal_ics){
        cur_ix<- internal_ics[i]
        for (offspring_ix in wtree[i, c(2,3)]){
            if (offspring_ix != 0){
                # tr$edge[iedge,]<-c(tra[i],tra[offspring_ix])
                tr$edge[iedge,]<-c(tra[i],tra[which(wtree[,5] == offspring_ix)])
                tr$edge.length[iedge]<-wtree[which(wtree[,5] == offspring_ix),1]-wtree[i,1]
                iedge<-iedge+1
            }
        } 
        # tr$edge[iedge,]<-c(tra[i],tra[wtree[i,2]])
        # tr$edge.length[iedge]<-wtree[wtree[i,2],1]-wtree[i,1]
        # iedge<-iedge+1
        # tr$edge[iedge,]<-c(tra[i],tra[wtree[i,3]])
        # tr$edge.length[iedge]<-wtree[wtree[i,3],1]-wtree[i,1]
        # iedge<-iedge+1
    }
    # tr$root.time=min(wtree[,1])
    tr$root.edge<- root.edge.length
    class(tr)<-'phylo'
    tr$tip.label<- wtree[tip_ics, 5]
    tr$node.label<- wtree[internal_ics, 5]
    
    ori<- wtree[nrow(wtree), 5]
    ori.time<- wtree[nrow(wtree), 1]
    
    return(list(phylo= tr, ori= ori, ori.time= ori.time))
}
