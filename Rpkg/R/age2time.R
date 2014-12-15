age2time <- function(node.ages,eps=-Inf) {
  btimes <- node.ages[node.ages[,4]>0,3]
  stimes <- node.ages[node.ages[,4]==0,3]
  out <- rbind(cbind(btimes,1),cbind(stimes,0))
  max.time <- max(out[,1])
  out[,1] <- max.time - out[,1]
  out <- out[out[,1]>eps,]
  srt <- order(out[,1])
  out[srt,]
}

