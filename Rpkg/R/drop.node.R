drop.node <- function(node.ages,to.drop) {
  for (node in to.drop) {
    i <- which(node.ages[,2]==node)
    parent.code <- node.ages[i,1]
    parent.i <- which(node.ages[,2]==parent.code)
    child.i <- which(node.ages[,1]==parent.code)
    sister.i <- child.i[!(child.i %in% i)]
    if (length(sister.i) == 1) {
      # kill parent
      node.ages[sister.i,1] <- node.ages[parent.i,1]
      node.ages[sister.i,6] <- node.ages[parent.i,6]
      node.ages <- node.ages[-c(parent.i,i),]
    }
  }
  node.ages
}

