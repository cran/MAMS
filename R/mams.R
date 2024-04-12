mams <- function(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.75, p0=0.5, delta=NULL, delta0=NULL, sd=NULL,
                 ushape="obf", lshape="fixed", ufix=NULL, lfix=0, nstart=1, nstop=NULL, sample.size=TRUE, N=20,
                 type="normal", parallel=TRUE, print=TRUE){

  #require(mvtnorm) ## the function pmvnorm is required to evaluate multivariate normal probabilities

  ############################################################################################
  ## 'mesh' creates the points and respective weights to use in the outer quadrature integrals
  ## you provide one-dimensional set of points (x) and weights (w) (e.g. midpoint rule, gaussian quadrature, etc)
  ## and the number of dimensions (d) and it gives d-dimensional collection of points and weights
  ###############################################################################################

  mesh<-function(x,d,w=1/length(x)+x*0){
    n<-length(x)
    W<-X<-matrix(0,n^d,d)
    for (i in 1:d){
      X[,i]<-x
      W[,i]<-w
      x<-rep(x,rep(n,length(x)))
      w<-rep(w,rep(n,length(w)))
    }
    w<-exp(rowSums(log(W)))
    list(X=X,w=w)
  }


  R_prodsum1 <-function(x,l,u,r,r0,r0diff,J,K,Sigma){
    ########################################################################################################
    ## x is vector of dummy variables ( the t_j in gdt paper ), l and u are boundary vectors
    ########################################################################################################
    .Call("C_prodsum1",
       x2 = as.double(x), l2 = as.double(l), u2 = as.double(u),
       r2 = as.double(r), r02 = as.double(r0), r0diff2 = as.double(r0diff),
       J2 = as.integer(J), K2 = as.integer(K), Sigma2 = Sigma,
       maxpts2 = as.integer(25000), releps2 = as.double(0),
       abseps2 = as.double(0.001), tol2 = as.double(0.0001))
  }

  #############################################################################################################
  ##  'typeI' performs the outer quadrature integral in type I error equation using 'mesh' and 'prodsum'
  ##  and calculates the difference with the nominal alpha.
  ##  The accuracy of the quadrature integral will depend on the choice of points and weights.
  ##  Here, the number of points in each dimension is an input, N.
  ##  The midpoint rule is used with a range of -6 to 6 in each dimension.
  #############################################################################################################

  typeI <-function(C,alpha,N,r,r0,r0diff,J,K,Sigma,mmp,
                   ushape,lshape,lfix=NULL,ufix=NULL,parallel=parallel,print=print){

    ########################################################################
    ## the form of the boundary constraints are determined as functions of C.
    ########################################################################
    if(!is.function(ushape)){
      if (ushape=='obf'){
        u<-C*sqrt(r[J]/r)
      }
      else if (ushape=='pocock'){
        u<-rep(C,J)
      }
      else if (ushape=='fixed'){
        u<-c(rep(ufix,J-1),C)
      }
      else if (ushape=='triangular') {
        u<-C*(1+r/r[J])/sqrt(r)
      }
    }else{
      u <- C*ushape(J)
    }

    if(!is.function(lshape)){
      if (lshape=='obf'){
        l<- c(-C*sqrt(r[J]/r[1:(J-1)]),u[J])
      }
      else if (lshape=='pocock'){
        l<-c(rep(-C,J-1),u[J])
      }
      else if (lshape=='fixed'){
        l<-c(rep(lfix,J-1),u[J])
      } else if (lshape=='triangular') {
        if(ushape=="triangular"){
          l<--C*(1-3*r/r[J])/sqrt(r)
        }else{
          l<--C*(1-3*r/r[J])/sqrt(r)/(-1*(1-3)/sqrt(J))
        }
      }
    }else{
      l <- c(C*lshape(J)[1:(J-1)],u[J])
    }


    if(parallel){
        evs <- future.apply::future_apply(mmp$X,1,R_prodsum1,
                    l=l,u=u,r=r,r0=r0,r0diff=r0diff,J=J,K=K,Sigma=Sigma,
                    future.seed=TRUE,future.packages="MAMS")
    }else{
        evs <- apply(mmp$X,1,R_prodsum1,
                    l=l,u=u,r=r,r0=r0,r0diff=r0diff,J=J,K=K,Sigma=Sigma)
    }
    if(print){
        message(".",appendLF=FALSE)
    }
    truealpha<-1-mmp$w%*%evs
    return(truealpha-alpha)
  }
  if(print){message("\n",appendLF=FALSE)}


  ####################################################################################
  ## 'R_prodsum2' evaluates the integrand of Pi_1 according to  gdt paper:
  #####################################################################################

  R_prodsum2 <-function(x,r,r0, u,K,delta,delta0,n,sig){
    .Call("C_prodsum2",
          x2 = as.double(x), r2 = as.double(r), r02 = as.double(r0),
          K2 = as.double(K), u2 = as.double(u),
          delta2 = as.double(delta), delta02 = as.double(delta0),
          n2 = as.double(n), sig2 = as.double(sig))
  }


  ############################################################################################
  ## 'R_prodsum3' evaluates the integrand of Pi_j for j>1.
  ############################################################################################

  R_prodsum3 <-function(x,l,u,r,r0,r0diff,J,K,delta,delta0,
                        n,sig,Sigma,SigmaJ){
    .Call("C_prodsum3",
          x2 = as.double(x), l2 = as.double(l), u2 = as.double(u),
          r2 = as.double(r), r02 = as.double(r0), r0diff2 = as.double(r0diff),
          Jfull2 = as.integer(J), K2 = as.integer(K),
          delta2 = as.double(delta), delta02 = as.double(delta0),
          n2 = as.double(n), sig2 = as.double(sig),
          Sigma2 = Sigma, SigmaJ2 = SigmaJ,
          maxpts2 = as.integer(25000), releps2 = as.double(0),
          abseps2 = as.double(0.001), tol2 = as.double(0.001))
  }

  ######################################################################################################
  ##  'typeII' performs the outer quadrature integrals of Pi_j j=1,...,J  using 'mesh' (stored in mmp_j),
  ##  'prodsum2' and 'prodsum3' as well as summing the Pi_1,...,Pi_J and calculates the difference
  ##  with the nominal power.
  ##  The accuracy of the quadrature integral again depends on the choice of points and weights.
  ##  Here, the number of points in each dimension is an input, N.
  ##  The midpoint rule is used with a range of -6 to 6 in each dimension.
  ########################################################################################################

  typeII <-function(n,beta,l,u,N,r,r0,r0diff,J,K,delta,delta0,sig,Sigma,mmp_j,parallel=parallel){

    evs <- apply(mmp_j[[1]]$X,1,R_prodsum2,
                 r=r,r0=r0,K=K,u=u,delta=delta,delta0=delta0,n=n,sig=sig)
    pi<-mmp_j[[1]]$w%*%evs

    if(J>1){
      for (j in 2:J){
        A<-diag(sqrt(r[j]/(r[j]-r[1:(j-1)])),ncol=j-1)
        SigmaJ<-A%*%(Sigma[1:(j-1),1:(j-1)]-Sigma[1:(j-1),j]%*%t(Sigma[1:(j-1),j]))%*%A
        if(parallel){
            evs <- future.apply::future_apply(mmp_j[[j]]$X,1,R_prodsum3,
                     l=l,u=u,r=r,r0=r0,r0diff=r0diff,J=j,K=K,delta=delta,delta0=delta0,
                     n=n,sig=sig,Sigma=Sigma[1:j,1:j],SigmaJ=SigmaJ,
                     future.seed=TRUE, future.packages="MAMS")
        }else{
            evs <- apply(mmp_j[[j]]$X,1,R_prodsum3,
                     l=l,u=u,r=r,r0=r0,r0diff=r0diff,J=j,K=K,delta=delta,delta0=delta0,
                     n=n,sig=sig,Sigma=Sigma[1:j,1:j],SigmaJ=SigmaJ)
        }
        pi<-pi+mmp_j[[j]]$w%*%evs
      }
    }
    return(1-beta-pi)
  }

  ## checking input parameters
  if(K%%1>0 | J%%1>0){stop("K and J need to be integers.")}
  if(K<1 | J<1){stop("The number of stages and treatments must be at least 1.")}
  if(N<=3){stop("Number of points for integration by quadrature to small or negative.")}
  if(N>3 & N<=10){warning("Number of points for integration by quadrature is small which may result in inaccurate solutions.")}
  if(alpha<0 | alpha>1 | power<0 | power>1){stop("Error rate or power not between 0 and 1.")}
  if(length(r)!=length(r0)){stop("Different length of allocation ratios on control and experimental treatments.")}
  if(length(r)!=J){stop("Length of allocation ratios does not match number of stages.")}

  if(is.numeric(p) & is.numeric(p0) & is.numeric(delta) & is.numeric(delta0) & is.numeric(sd)){
    stop("Specify the effect sizes either via (p, p0) or via (delta, delta0, sd) and set the other parameters to NULL.")
  }

  if(is.numeric(p) & is.numeric(p0)){
    if(p<0 | p>1 | p0<0 | p0>1){stop("Treatment effect parameter not within 0 and 1.")}
    if(p<=p0){stop("Interesting treatment effect must be larger than uninteresting effect.")}
    if(p0<0.5){warning("Uninteresting treatment effect less than 0.5 which implies that reductions in effect over placebo are interesting.")}
  }else{
    if(is.numeric(delta) & is.numeric(delta0) & is.numeric(sd)){
      if(sd<=0){stop("Standard deviation must be positive.")}
    }else{
      stop("Specify the effect sizes either via (p, p0) or via (delta, delta0, sd).")
    }
  }

  if(is.function(ushape) & is.function(lshape)){
    warning("You have specified your own functions for both the lower and upper boundary. Please check carefully whether the resulting boundaries are sensible.")
  }

  if(!is.function(ushape)){
    if(!ushape%in%c("pocock","obf","triangular","fixed")){stop("Upper boundary does not match the available options.")}
    if(ushape=="fixed" & is.null(ufix)){stop("ufix required when using a fixed upper boundary shape.")}
  }else{
    b <- ushape(J)
    if(!all(sort(b,decreasing=TRUE)==b)){stop("Upper boundary shape is increasing.")}
  }
  if(!is.function(lshape)){
   if(!lshape%in%c("pocock","obf","triangular","fixed")){stop("Lower boundary does not match the available options.")}
   if(lshape=="fixed" & is.null(lfix)){stop("lfix required when using a fixed lower boundary shape.")}
  }else{
    b <- lshape(J)
    if(!all(sort(b,decreasing=FALSE)==b)){stop("Lower boundary shape is decreasing.")}
  }

  ############################################################################
  ## Convert treatment effects into absolute effects with standard deviation 1:
  ############################################################################

  if(is.numeric(p) & is.numeric(p0)){
    delta <- sqrt(2) * qnorm(p)
    delta0 <- sqrt(2) * qnorm(p0)
    sig <- 1
  }else{
    delta <- delta
    delta0 <- delta0
    p0 <- pnorm(delta0/sqrt(2 * sd^2)) # for subsequent if(J==1 & p0==0.5)
    sig <- sd
  }

  ############################################################################
  ## gaussian quadrature's grid and weights for stages 1:J
  ############################################################################

  mmp_j = as.list(rep(NA,J))
  for(j in 1:J){
      mmp_j[[j]] = mesh(x=(1:N-.5)/N*12-6,j,w=rep(12/N,N))
  }


  ############################################################################
  ## Ensure equivalent allocation ratios yield same sample size
  ############################################################################

  h <- min(c(r,r0))
  r <- r/h
  r0 <- r0/h

  ####################################################################
  ## Create the variance covariance matrix from allocation proportions:
  ####################################################################

  bottom<-matrix(r,J,J)
  top<-matrix(rep(r,rep(J,J)),J,J)
  top[upper.tri(top)]<-t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)]<-t(bottom)[upper.tri(bottom)]
  Sigma<-sqrt(top/bottom)

  ###############################################################################
  ## Create r0diff: the proportion of patients allocated to each particular stage
  ###############################################################################

  r0lag1<-c(0,r0[1:J-1])
  r0diff<-r0-r0lag1


  ################################
  ## Find boundaries using 'typeI'
  ################################

  if(print){
     message("   i) find lower and upper boundaries\n      ",appendLF=FALSE)
  }

  # Quick & dirty fix to enable single-stage design with specification lshape="obf"
  if(!is.function(lshape)){
    if(J==1 & lshape=="obf"){
      lshape <- "pocock"
    }
  }

  uJ<-NULL
  ## making sure that lfix is not larger then uJ
  try(uJ<-uniroot(typeI,c(qnorm(1-alpha)/2,5),alpha=alpha,N=N,r=r,r0=r0,r0diff=r0diff,
                  J=J,K=K,Sigma=Sigma,mmp=mmp_j[[J]],
                  ushape=ushape,lshape=lshape,lfix=lfix,ufix=ufix,parallel=parallel,print=print,
                  tol=0.001)$root, silent=TRUE)
  #if(is.null(uJ)){stop("Lower boundary (lfix) is too large.")}
  if(is.null(uJ)){stop("No boundaries can be found.")}

  if(!is.function(ushape)){
    if (ushape=='obf'){
      u<-uJ*sqrt(r[J]/r)
    }
    else if (ushape=='pocock'){
      u<-rep(uJ,J)
    }
    else if (ushape=='fixed'){
      u<-c(rep(ufix,J-1),uJ)
    } else if (ushape=='triangular') {
      u<-uJ*(1+r/r[J])/sqrt(r)
    }
  }else{
    u <- uJ*ushape(J)
  }

  if(!is.function(lshape)){
    if (lshape=='obf'){
      l<- c(-uJ*sqrt(r[J]/r[1:(J-1)]),u[J])
    }
    else if (lshape=='pocock'){
      l<-c(rep(-uJ,J-1),u[J])
    }
    else if (lshape=='fixed'){
      l<-c(rep(lfix,J-1),u[J])
    } else if (lshape=='triangular') {
       if(ushape=="triangular"){
          l<--uJ*(1-3*r/r[J])/sqrt(r)
       }else{
          l<--uJ*(1-3*r/r[J])/sqrt(r)/(-1*(1-3)/sqrt(J))
       }
    }
  }else{
    l <- c(uJ*lshape(J)[1:(J-1)],u[J])
  }

  #########################################################
  ## Find alpha.star
  #########################################################

  if(print){
     message("\n  ii) define alpha star\n",appendLF=FALSE)
  }

  alpha.star <- numeric(J)
  alpha.star[1] <- typeI(u[1],
                   alpha = 0, N = N, r = r[1], r0 = r0[1],r0diff = r0diff[1],
                   J = 1, K = K, Sigma = Sigma, mmp = mmp_j[[1]],
                   ushape = "fixed",lshape = "fixed",lfix = NULL, ufix = NULL,
                   parallel = parallel, print = FALSE
                   )
  if (J > 1){
      for (j in 2:J){
          alpha.star[j] <- typeI(u[j],
                           alpha = 0, N = N, r = r[1:j], r0 = r0[1:j],r0diff = r0diff[1:j],
                           J = j, K = K, Sigma = Sigma, mmp = mmp_j[[j]],
                           ushape = "fixed",lshape = "fixed",lfix = l[1:(j - 1)],ufix = u[1:(j - 1)],
                           parallel = parallel, print = FALSE
                           )
      }
  }

  #############################################################
  ##  Now find samplesize for arm 1 stage 1 (n)  using 'typeII'.
  ##  Sample sizes for all stages are then determined by
  ##  r*n and r0*n.
  #############################################################

  if(J==1 & p0==0.5) {
    if(r0>r){
      r <- r/r0
      r0 <- r0/r0
    }
    rho <- r / (r + r0)
    corr <- matrix(rho, K, K) + diag(1 - rho, K)
    if (K==1){
      quan <- qmvnorm(1-alpha, mean=rep(0, K), sigma=1)$quantile
    } else {
      quan <- qmvnorm(1-alpha, mean=rep(0, K), corr=corr)$quantile
    }
    n <- ((quan + qnorm(power)) / (qnorm(p)*sqrt(2)))^2 *(1+1/r)

  }else{

    n <- nstart
    ###################################################################################################
    ## program could be very slow starting at n=0, may want to start at more sensible lower bound on n
    ## unlike type I error, power equation does not neccessarily have unique solution n therefore search
    ## for smallest solution:
    ####################################################################################################

    pow <- 0
    if(sample.size){
      if(print){
         message(" iii) perform sample size calculation\n",appendLF=FALSE)
      }
      if(is.null(nstop)){
        nx <- nstart
        po <- 0
        while (po==0){
          nx <- nx + 1
          po <- (typeII(nx, beta=1 - power, l=l, u=u, N=N, r=r, r0=r0, r0diff=r0diff, J=1,
                        K=K, delta=delta, delta0=delta0, sig=sig, Sigma=Sigma, mmp_j = mmp_j,
                        parallel = parallel
                        ) < 0)
        }
        nstop <- 3 * nx
      }
      if(print){
         message(paste0("      (maximum iteration number = ",nstop-nstart+1,")\n      "),appendLF=FALSE)
      }
      while (pow==0 & n<=nstop){
        n <- n + 1
        pow <- (typeII(n, beta=1 - power, l=l, u=u, N=N, r=r, r0=r0, r0diff=r0diff, J=J,
                      K=K, delta=delta, delta0=delta0, sig=sig, Sigma=Sigma, mmp_j = mmp_j,
                      parallel = parallel
                      ) < 0)
        if(print){
            if(any(seq(0,nstop,50)==n)){message(n,"\n      ",appendLF=FALSE)}else{message(".",appendLF=FALSE)}
        }
      }
      if(print){
         message("\n",appendLF=FALSE)
      }
      if((n - 1)==nstop){warning("The sample size search was stopped because the maximum sample size (nstop, default: 3 times the sample size of a fixed sample design) was reached.\n")}

    }else{
      n <- NULL
    }
  }

  res <- NULL
  res$l <- l
  res$u <- u
  res$n <- n

  res$rMat <- rbind(r0,matrix(r,ncol=J,nrow=K,byrow=TRUE)) ## allocation ratios
  res$N <- sum(ceiling(res$rMat[,J]*res$n)) ## maximum total sample sizeres$N <- K*r[J]*n+r0[J]*n ## maximum total sample size

  res$K <- K
  res$J <- J
  res$alpha <- alpha
  res$alpha.star <- alpha.star
  if(sample.size){
    res$power <- power
  }else{
    res$power <- NA
  }

  res$type <- type

  class(res)<-"MAMS"

  return(res)

}


