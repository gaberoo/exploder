epi2igraph <- function(epi) {
  require(igraph)
  edges <- cbind(epi$parent,epi$id,epi$itimes-epi$dtimes)
  graph.edgelist(edges[-1,1:2]+1)
}

