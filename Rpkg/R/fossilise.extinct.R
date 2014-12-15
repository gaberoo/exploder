fossilise.extinct <- function(tree,fossils=NULL,foss.rate=NA,
                              node.ages=NULL,tp=-Inf,t0=NULL) 
{
  if (is.null(fossils) & is.na(foss.rate)) {
    stop("Please supply either fossils or fossilization rate.")
  }
  if (is.null(node.ages)) {
    node.ages <- node.age(tree,tp,t0)
  }
  if (is.null(fossils)) {
    fossils <- add.fossils(tree,foss.rate,node.ages,tp,t0)
  }
  # order fossils by fossil age
  foss.o <- order(fossils[,2],decreasing=TRUE)
  fossils <- cbind(fossils,NA)
  for (j in foss.o) {
    f <- fossils[j,]
    # get node ID
    code <- f[1]
    # get row in node.ages
    i <- which(node.ages[,2]==code)
    # check if node is in an extinct clade
    if (node.ages[i,5] == 0) {
      to.remove <- c()
      # check if it has decendants
      if (node.ages[i,4] > 0) {
        # remove all descendants of this node
        to.remove <- rbind(to.remove,get.tips(node.ages,code))
        # set number of descendent to zero
        node.ages[i,4] <- 0
      }
      # set age to fossilization time
      node.ages[i,3] <- f[2]
      # update parent extant subclade count
      parents <- get.parents(node.ages,code)
      add.extinct <- 1*(node.ages[parents[-nrow(parents),1],5]==0)
      node.ages[parents[-1,1],5] <- node.ages[parents[-1,1],5] + add.extinct
      # set clade to surviving
      node.ages[i,5] <- 1
      # remove extinct subclade
      if (length(to.remove) > 0) {
        node.ages <- node.ages[-to.remove[,1],]
      }
      fossils[j,5] <- 2
    } else {
      fossils[j,5] <- 3
    }
  }
  rownames(fossils) <- NULL
  colnames(fossils) <- c("code","foss.time","coords","row.id","foss.type")
  list(fossils=fossils,node.ages=node.ages)
}

