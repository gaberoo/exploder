det.lineages <- function(forest,pars,extant.num,extant.root.time=NULL,
                         ode.times=NULL,length.out=100,scond=FALSE,
                         root.lineages=1) {
  require(deSolve)
  if (pars[5] > 0.0) {
    n <- pars[5]*extant.num
  } else {
    n <- mean(sapply(forest,function(tree) sum(tree[,2])))
  }
  r <- pars[2]-(pars[3]+pars[4])
  if (is.null(extant.root.time)) {
    if (pars[4] > 0.0) {
      if (pars[1] > 0) {
        bn <- pars[2]/pars[1]
        extant.root.time <- log((exp(bn*n/pars[4])-1)*r/bn)/r
      } else {
        extant.root.time <- log(n*r/pars[4]+1)/r
      }
    } else {
      extant.root.time <- mean(sapply(forest,function(tree) max(tree[,1])))
    }
  }
  if (is.null(ode.times) & length.out > 0)
    ode.times <- seq(0,extant.root.time,length.out=length.out)
  #if (pars[5] == 0) {
    Y <- c(L=extant.num,S=0,C=0)
    ode(Y,ode.times,func="derivs",
        parms=c(unlist(pars),extant.root.time,1.0*scond,1.0*root.lineages),
        method="ode45",jacfunc="jac",dllname="expoTree",initfunc="initmod",
        nout=3,outnames=c("samp","coal","It"))
  #} else {
  #  Y <- c(L=extant.num,S=0,C=0,I=extant.num/rho)
  #  ode(Y,ode.times,func="derivsli",
  #      parms=c(unlist(pars),extant.root.time,1.0*scond,1.0*root.lineages),
  #      method="ode45",dllname="expoTree",initfunc="initmod",
  #      nout=2,outnames=c("samp","coal"))
  #}
}


