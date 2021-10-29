#' given a number of offsprings and effective sizes, 
#' this script simulate back
#' For the subsequent simulation that considers LGT, 
#' 

# Simulates the within-host coalescent model
# @param times times at which N samples are taken(counted forward in time from infection time) 
# @param neg is the product of the within-host effective population size and the generation duration in days 
# @return array of size(2N)*3 where each row is a node,the first column indicate the date of the node and the last two columns indicate the two children. This array has internal nodes sorted in order of most recent to most ancient node(and remains so during the algorithm). The last node corresponds to infection time and only has one child 
# withinhost = function(times,neg)  {
wtree<- wtrees[[1]]
nodes<- wtree 
nodes[,2]<- match(nodes[,2], nodes[,5]);nodes[is.na(nodes[,2]),2]<- 0;
nodes[,3]<- match(nodes[,3], nodes[,5]);nodes[is.na(nodes[,3]),3]<- 0;
nodes<- nodes[, 1:3]

# start_date<- wtree[nrow(wtree),1]
# wtree[,1]<- wtree[,1]-start_date
# tips
tip_ics<-which((!wtree[,2] %in% wtree[,5]) & (!wtree[,3] %in% wtree[,5]))
tips<- wtree[tip_ics, ]
times<- tips[,1]
pseudo_samples_size<- 3
times<- c(times, rep(min(times)-1e-6, pseudo_samples_size))

# internal nodes
internal_ics<- setdiff(1:(nrow(wtree)-1), tip_ics)
internal_nodes<- matrix()
if (length(internal_ics)==1){
    internal_nodes<- matrix(wtree[internal_ics, ], nrow=1)
}else{
    internal_nodes<- wtree[internal_ics, ]
}
internal_nodes[,2]<- match(internal_nodes[,2], tips[,5])
internal_nodes[,3]<- match(internal_nodes[,3], tips[,5])
# ori
ori<- matrix(wtree[nrow(wtree),], nrow= 1)
ori[,2]<- match(ori[,2], wtree[, 5])

# initial state
nodes<- rbind(internal_nodes[, 1:3], ori[, 1:3])
starting_time<- nodes[nrow(nodes), 1]
nodes[nrow(nodes), 1]<- 0
node_ix_top<- 749
# wtree<-wtrees[[1]]
# sample.times<- wtree[1:(nrow(wtree)-1), 1]
neg<- 100/365
# times<- rep(sample.times, effective_size) # except for the one that was already included
# nodes<- wtree[1:(nrow(wtree)-1),1:3]
# nodes[,1]<- 0
#  
  prob <- 0 
  ind=order(times,decreasing=T);tim=times[ind]
  n <- length(tim)
  # n <- length(tim) + nrow(nodes)
  # nodes <- cbind(0,ind[1],0);#Start with one node at time 0 and with the first isolate connected to it
  # i <- 2
  i<- nrow(nodes)+1
  while (i <= n) {#Graft branches one by one 
    curt <- tim[i]
    fi <- which( nodes[ ,1] < curt );fi<-fi[1]
    trunc=0
    for (j in (fi:nrow(nodes)))  {
        trunc=trunc+(curt-nodes[j,1]) * (i-j) 
        curt <- nodes[j,1] 
    }
    r=rexpT(1/neg,trunc)
    prob=prob+r$prob
    r=r$x
    curt <- tim[i]#Current time:start with date of isolate and go back in time until coalescence happens 
    fi <- which( nodes[ ,1] < curt );fi<-fi[1]
    for (j in (fi:nrow(nodes)))  {
      if (r > (curt-nodes[j,1]) * (i-j))  { 
        r <- r-(curt-nodes[j,1]) * (i-j) 
        curt <- nodes[j,1] 
      } else { 
        curt <- curt-r/(i-j)#Found the time for grafting
        r <- 0 
        break 
      } 
    } 
    if (r>0) print('error: this should not happen')
    #Create new node 
    a <- nodes[ ,2:3];a[a >= j + n] <- a[a >= j + n] + 1;nodes[ ,2:3] <- a;#Renumbering according to table insertion in next line 
    nodes <- rbind(nodes[seqML(1,j-1), ],
                   c(curt,ind[i],0),
                   nodes[seqML(j,nrow(nodes)),]) 
    #Now choose on which branch to regraft amongst the branches alive at time curt 
    no <- j 
    side <- 2 
    #prob <- prob + log(1/(nrow(nodes)-j))
    w <- 1 + floor(runif(1) * (nrow(nodes)-j)) 
    while (w > 0)  { 
      no <- no + side-1 
      side <- 3-side 
      if (nodes[no,side + 1] <= n ||(nodes[no,side + 1] > n && nodes[nodes[no,side + 1]-n,1] > curt))  { 
        w <- w-1 
      } 
    } 
    nodes[j,3] <- nodes[no,side + 1] 
    nodes[no,side + 1] <- n + j 
    i <- i + 1 
  } 
  # nodes<-cbind(nodes, wtree[match(nodes[,1], wtree[,1]), c(4,5)])
  # nodes[which(is.na(nodes[,4])), 4]<- wtree[1,4] # patient id
  # pseudo tips
  pseudo_tips<- matrix(0, nrow=pseudo_samples_size, ncol= 3)
      
      matrix(c(times, rep(min(times)-1e-6, pseudo_samples_size))
  p_tips<- tips[match(times[(nrow(tips)+1):length(times)], tips[,1]), ]
  p_tips[,2]<- 0
  p_tips[,3]<- 0
  p_tips[,5]<- NA
  
  nodes<- rbind(tips, p_tips, nodes)
 # change id to transmission tree node ids
  # global_ix<- nodes[,5]
  # global_ix[which(is.na(global_ix))]<- sapply(which(is.na(global_ix)), function(x) paste0('pseudo', as.character(x)))
  nodes[which(is.na( nodes[,5])),5]<- sapply(which(is.na( nodes[,5])), function(x) node_ix_top+x)
  # tips were done
  # so just update the internal ones
  # offspring1<- global_ix[nodes[(length(times)+1):nrow(nodes), 2]]
  # offspring2<- global_ix[nodes[(length(times)+1):(nrow(nodes)-1), 3]]
  nodes[(length(times)+1):nrow(nodes), 2]<- nodes[nodes[(length(times)+1):nrow(nodes), 2],5]
  nodes[(length(times)+1):(nrow(nodes)-1), 3]<- nodes[nodes[(length(times)+1):(nrow(nodes)-1), 3],5]
  # nodes[,5]<- global_ix
  # nodes[(length(times)+1):nrow(nodes), 2]<- offspring1
  # nodes[(length(times)+1):(nrow(nodes)-1), 3]<- offspring2
  # nodes_df<- as.data.frame(nodes, stringsAsFactors = F)
  # nodes_df[,1]<- as.numeric(nodes_df[,1])
  # nodes[nrow(nodes),]<- wtree[nrow(wtree),]
  
  return(list(nodes = nodes,prob = prob))
# } 

  
  #truncated exponential distribution
  rexpT = function(rate,trunc) {
      x=-log(1-runif(1)*(1-exp(-rate*trunc)))/rate
      prob=log(rate)-rate*x-log(1-exp(-rate*trunc))
      return(list(x=x,prob=prob))
  }
  