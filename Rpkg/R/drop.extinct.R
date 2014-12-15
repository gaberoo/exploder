drop.extinct <- function(tree,node.ages=NULL,tp=-Inf,t0=NULL) {
  if (is.null(node.ages)) node.ages <- node.age(tree,tp,t0)
  to.drop <- as.integer(node.ages[node.ages[,4]==0 & node.ages[,5]==0,2])
  if (length(to.drop) > 0) {
    return(drop.tip(tree,to.drop))
  } else {
    return(tree)
  }
}

