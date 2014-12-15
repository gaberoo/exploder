ape2time <- function(tree=NULL,node.ages=NULL,eps=0.0,tp=-Inf,t0=NULL) {
  if (is.null(node.ages)) node.ages <- node.age(tree,tp,t0)
  
  bev <- node.ages[,4]>0 & node.ages[,5]>1    # branching events
  sev <- node.ages[,4]==0 & node.ages[,5]>0   # sampling events

  btimes <- as.numeric(node.ages[bev,3])
  stimes <- as.numeric(node.ages[sev,3])

  btimes.id <- as.numeric(node.ages[bev,2])
  stimes.id <- as.numeric(node.ages[sev,2])

  max.time <- max(btimes,stimes)
  btimes <- max.time - btimes
  stimes <- max.time - stimes

  bt <- cbind(btimes,1,btimes.id)
  st <- cbind(stimes,0,stimes.id)

  bt.children <- t(sapply(btimes.id,function(code) {
    x <- node.ages[node.ages[,"parent"]==code,"code"]
    x[!is.na(x)]
  }))

  bt <- cbind(bt,bt.children)
  st <- cbind(st,-1,-1)

  root <- which(is.na(node.ages[,1]))
  if (length(root) > 0) {
    root.code <- node.ages[root,"code"]
    root.age <- max.time+node.ages[root,3]-node.ages[root,6]
    bt <- rbind(bt,c(root.age,1,0,root.code,-1))
  }

  out <- rbind(bt,st)
  out <- out[order(out[,1]),]

  extant <- out[,1] <= eps
  out[!extant,]
}

