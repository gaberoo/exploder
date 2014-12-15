get.parents <- function(node.ages,code) {
  i <- which(node.ages[,2]==code)
  parent <- node.ages[i,1]
  if (! is.na(parent)) {
    return(rbind(cbind(i,parent),get.parents(node.ages,parent)))
  } else {
    return(cbind(i,parent))
  }
}

