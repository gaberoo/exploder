plot.edgePhylo <- function(tree=NULL,node.ages=NULL,tp=-Inf,t0=NULL,
                           cols=brewer.pal(3,"Set1"),lwd=1,
                           fossils=NULL,axes=TRUE,xlab="Time")
{
  require(RColorBrewer)
  if (is.null(node.ages) & is.null(tree)) {
    stop("Supply either a phylo or nodeAge object.")
  }
  if (is.null(node.ages)) {
    nedge <- nrow(tree$edge)+1
    node.ages <- node.age(tree,tp,t0)
  } else {
    nedge <- nrow(node.ages)
  }
  edge.col <- rep(cols[1],nedge)
  edge.col[node.ages[,5]>0] <- cols[2]
  tips <- node.ages[,4]==0
  parent.coords <- c(NA,as.numeric(sapply(node.ages[-1,1],
    function(x) node.ages[node.ages[,2]==x,7])))
  xlim <- range(node.ages[,c(3,6)])
  ylim <- range(node.ages[,7])
  plot(1,type="n",xlim=xlim,ylim=ylim,axes=FALSE,xlab=xlab,ylab=NA)
  segments(node.ages[,3],node.ages[,7],
           node.ages[,6],node.ages[,7],col=edge.col)
  segments(node.ages[,6],node.ages[,7],
           node.ages[,6],parent.coords,col=edge.col)
  if (! is.null(fossils)) {
    points(fossils[,2],fossils[,3],pch=16)
  }
  if (axes) axis(1,tcl=-.3)
  if (is.finite(tp)) abline(v=max(node.ages[,3])+tp,lwd=1,lty=2,col=rgb(.5,.5,.5))
  invisible(cbind(node.ages,parent.coords))
}

