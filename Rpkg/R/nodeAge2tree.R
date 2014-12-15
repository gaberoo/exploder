nodeAge2tree <- function(node.ages) {
  root.id <- which(node.ages[,1]==NA)
  root.edge <- NULL
  if (length(root.id) == 1) {
    root.edge <- edge.len[rood.id]
  }
  edge.list    <- node.ages[-root.id,1:2]
  edge.lenngth <- node.ages[-root.id,3]-node.ages[-root.id,6]
  tip.labels   <- as.character(node.ages[,2])
  # relabel edges
  tips  <- node.ages[node.ages[,4]==0,2]
  nodes <- node.ages[node.ages[,4]>0 ,2]
  # NOT COMPLETE -- WILL NOT WORK !!
  tree <- list(edge=edge.list,edge.length=edge.len,
               tip.label=tip.labels,
               Nnode=max(edge.list[[1]])-max(epi$id),root.edge=root.edge)
  class(tree) <- "phylo"
  tree
}
