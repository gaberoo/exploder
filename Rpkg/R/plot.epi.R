plot.epi <- function(epi) {
  coords <- epi.coords(epi)[[2]]
  coord.ids <- sapply(coords[,1],function(x) which(epi$id==x))
  parent.r <- sapply(epi$parent[coord.ids],function(x) {
    y <- which(coords[,1]==x)
    if (length(y) == 0) return(NA)
    else return(y)
  })
  coords <- cbind(coords,coord.ids,parent.r)
  colnames(coords) <- c("id","y","arr.id","parent.row")
  xlim <- range(-epi$times,-epi$dtimes)
  ylim <- range(coords[,2])
  plot(1,type="n",xlim=xlim,ylim=ylim,axes=FALSE)
  segments(-epi$itimes[coords[,3]],coords[,2],
           -epi$dtimes[coords[,3]],coords[,2])
  segments(-epi$itimes[coords[,3]],coords[,2],
           -epi$itimes[coords[,3]],coords[coords[,4],2])
  axis(1,tcl=-.3)
  invisible(coords)
}

