get.tip.coords <- function(node.ages) {
  get.coords <- function(id,next.tip) {
    node <- which(node.ages[,2] == id)
    off <- which(node.ages[,1] == id)
    if (length(off) == 0) {
      # node is a tip
      node.ages[node,7] <- next.tip
      return(list(next.tip=next.tip+1,next.tip,cbind(id=id,coord=next.tip)))
    } else {
      off.ids <- node.ages[off,2]
      coords <- c()
      oc <- c()
      for (of in off.ids) {
        off.coord <- get.coords(of,next.tip)
        next.tip <- off.coord[[1]]
        coords <- rbind(coords,off.coord[[3]])
        oc <- c(oc,off.coord[[2]])
      }
      coords <- rbind(coords,cbind(id,mean(oc)))
      node.ages[node,7] <- mean(oc)
      return(list(next.tip=next.tip,mean(oc),coords))
    }
  }
  gc <- get.coords(node.ages[which.min(node.ages[,3]),2],0)
  ii <- as.integer(sapply(gc[[3]][,1],function(x) which(node.ages[,2]==x)))
  node.ages[ii,7] <- as.numeric(gc[[3]][,2])
  list(node.ages,gc)
}

