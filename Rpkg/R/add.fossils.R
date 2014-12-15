add.fossils <- function(tree,foss.rate,node.ages=NULL,tp=-Inf,t0=NULL) {
  if (is.null(node.ages)) node.ages <- node.age(tree,tp,t0)
  fossils <- apply(node.ages,1,function(x,foss.rate) {
    foss.time <- c()
    next.t <- x[6]
    # add fossils between ancestor and next node
    while (1) {
      next.t <- next.t+log(1/runif(1))/foss.rate
      if (next.t > x[3]) break
      foss.time <- c(foss.time,next.t)
    }
    foss.time
  },foss.rate=foss.rate)
  foss.edge <- lapply(1:nrow(node.ages),function(i) {
    if (length(fossils[[i]]) > 0) {
      cbind(node.ages[i,2],fossils[[i]],node.ages[i,7],i)
    } else {
      NULL
    }
  })
  do.call(rbind,foss.edge)
}

