#' Converts a phylogenetic tree into an ape phylo object
#' revised from TransPhylo::phyloFromPtree
#' @param wtree_from_ctree phylogenetic tree; col1: time; col2-3: indices of offsprings; col4: patien index; col5: node/tip indices
#' @return phylo object
#' @export
phyloFromWTree <- function(wtree_from_ctree, root, root.edge.length) {
  # nam=wtree$nam
  # wtree=wtree$wtree
  # 
  # root<-which(wtree[1:(nrow(wtree)-1),1]==min(wtree[1:(nrow(wtree)-1),1]))
  # root.edge.length<- wtree[root,1] - min(wtree[,1])
  
  # n<-ceiling(nrow(wtree)/2)
  internal_ics<- which((wtree[1:(nrow(wtree)-1),2]!=0 | wtree[1:(nrow(wtree)-1),3]!=0))
  tip_ics<- setdiff(1:(nrow(wtree)-1), internal_ics)
  nam<- as.character(tip_ics)
  n<- length(internal_ics)
  
  # if (n==1) return(ape::read.tree(text='(1);'))
  if (length(tip_ics)==1) return(ape::read.tree(text=sprintf('(1:%f);', root.edge.length)))
  tr<-list()
  # tr$Nnode<-n-1
  tr$Nnode<-n
  tr$tip.label<-nam
  # tr$edge<-matrix(0,nrow(wtree)-1,2)
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
  for (i in internal_ics){
      for (offspring_ix in wtree[i,c(2,3)]){
          if (offspring_ix != 0){
              tr$edge[iedge,]<-c(tra[i],tra[offspring_ix])
              tr$edge.length[iedge]<-wtree[offspring_ix,1]-wtree[i,1]
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
  tr$root.time=min(wtree[,1])
  tr$root.edge<- root.edge.length
  class(tr)<-'phylo'
  tr$node.label<- ((Ntip(tr)+1):(Ntip(tr)+Nnode(tr)))
  return(tr)
}
