#' tfica with multicore computation to save time
#'
#' @param Xin: M x T data input.
#' @param M: number of independent components
#' @param Win: initial value for unmixing matrix
#' @param tol (=0.001): convergence criteria
#' @param maxit: maximum iteration
#' @param nmaxit: maximum number of restart

#' @param sopoc : Sub-octaves per octave calculated.
#' @param ncores (=4)
#' @import coloredICA
#' @import polspline
#' @import dplR
#' @import seewave
#' @return
#' @export
#'
#' @examples
tfica.multicore<-function (Xin,
                   M = dim(Xin)[1],
                   Win = diag(M),
                   tol = 1e-04,
                   maxit = 20,
                   nmaxit=5,
                   sopoc=0.1,
                   ncores=4,
                   maxknots=10)

  #sopoc :
{
  p = dim(Xin)[1]
  if (M > p) {
    stop("Number of sources must be less or equal than number \n  of variables")
  }
  # selected frequency index

  ### pre-whitening
  T = ncol(Xin)
  X.c = scale(t(Xin), center = TRUE, scale = FALSE) # T x M matrix
  svdcovmat = eigen(cov(X.c))#svd(Xc/sqrt(T))
  K = svdcovmat$vector %*% diag(1/sqrt(svdcovmat$value)) %*% t(svdcovmat$vector)
 # K = K[1:M, ]
  Xc = t(X.c %*% K)

  ### Initialization
  Wold = Wnew = Wtmp = Win
  X.fft = mvfft(t(Xc))

#  k0=6
#  fourier_factor <- (4 * pi)/(k0 + sqrt(2 + k0^2))
  X.morlet = mclapply(1:nrow(Xc),function(idx){morlet(y1 = Xc[idx,], x1 = 1:T, dj = sopoc, siglvl = 0.95)}, mc.cores=ncores)

  perid.x = X.morlet[[1]]$period

  S.mat = Wold %*% Xc
  #S.mat.ppt = Wold %*% t(X.fft)
  freqlength=T#floor(T/2)
  # nperid = length(perid.x)
  fft.indx=2:(freqlength)
  lspec.ests = do.call(rbind,
                       mclapply(1:nrow(S.mat),function(jjj){
                         obj=S.mat[jjj,]
                         mm = length(obj)
                         tmpfit<-c()
                         try(tmpfit<-lspec(obj,maxknots=maxknots))
                         if(is.null(tmpfit)){
                           tmpfit<-lspec(obj + rnorm(mm)*0.000001,maxknots=maxknots)
                         }#lspec(period=Mod((obj)[2:(mm/2)])^2/(2*pi*mm))
                         pred = dlspec( 1/perid.x* 2*pi,tmpfit) #dlspec( 1/perid.x *2*pi,tmpfit)
                         pred$m[is.na(pred$m)]<-0
                         lspec.est = (pred$d + pred$m) ## need to straightedn up.
                         return(lspec.est)}, mc.cores=ncores))
  print(apply(lspec.ests, 1, function(ss)sum(ss==0)))
#  plot(1/perid.x[241:1], lspec.ests[1,241:1])
#  plot(1/perid.x[241:1], lspec.ests[2,241:1])
#  plot(1/perid.x[241:1], lspec.ests[3,241:1])

  perall = mclapply(1:length(perid.x),
         function(xx){
           ttt=apply(do.call(cbind, lapply(X.morlet, function(ss)ss$wave[,xx])), 1,
                     function(xxx)as.vector(Re(xxx) %*% t(Re(xxx)) + Im(xxx) %*% t(Im(xxx))))
         }, mc.cores=ncores) ## we don't need to repeat this calculation

  Ixf = lapply(1:M, function(m){
    matrix(apply(do.call(cbind,
                         lapply(1:length(perid.x), function(jjj) apply(perall[[jjj]]/lspec.ests[m,jjj] + log(lspec.ests[m,jjj]),1,mean))),1,mean),
           M,M)
  })
#  Ixf = lapply(split(tmp,f=rep(1:M, each=M*M)), function(xx)Re(matrix(xx,M,M)))

  ### Optimization
  lim=1
  iter=0
  wlik = -Inf
  NInv = 0
  Ndec = 0

  while (lim > tol & iter < maxit & NInv < nmaxit) {
    iter = iter + 1
    taucount = 1
    err = 1
    orthoerror = 1
    tau = 0.5
    eigenval = rep(0, M)
    Wold <- Wnew
    ## update W
    while (taucount < 60 & err > 1e-05 & orthoerror > 1e-05) {
      for (j in 1:M) {
        Gam = 0
        if (j > 1) {
          for (k in 1:(j - 1)) {
            nu = matrix(Wnew[k, ], M, 1) %*% matrix(Wnew[k,], 1, M)
            Gam = Gam + nu
          }
        }

        tmpV = Ixf[[j]] + tau * Gam
        eigenv = eigen(tmpV)
        eigenval[j] = eigenv$values[M]
        Wnew[j, ] = eigenv$vectors[,M]
      }
      orthoerror = (sum(sum((Wnew %*% t(Wnew) - diag(rep(1, M)))^2)))
      err = amari_distance(rerow(Wnew), rerow(Wold))
      taucount = taucount + 1
      tau = 2 * tau
    }

    ## Update spectral density (when iter=1, it is initialization)
    S.mat = Wnew %*% Xc
    #S.mat.ppt = Wnew %*% t(X.fft)
    freqlength=T#floor(T/2)
    # nperid = length(perid.x)
    lspec.ests = do.call(rbind,
                         mclapply(1:M,function(jjj){
                           obj=S.mat[jjj,]
                           mm = length(obj)
                           tmpfit<-c()
                           try(tmpfit<-lspec(obj,maxknots=maxknots))
                           if(is.null(tmpfit)){
                             tmpfit<-lspec(obj + rnorm(mm)*0.000001,maxknots=maxknots)
                             }#lspec(period=Mod((obj)[2:(mm/2)])^2/(2*pi*mm))
                           pred = dlspec( 1/perid.x* 2*pi,tmpfit) #dlspec( 1/perid.x *2*pi,tmpfit)
                           pred$m[is.na(pred$m)]<-0
                           lspec.est = (pred$d + pred$m) ## need to straightedn up.
                           return(lspec.est)}, mc.cores=ncores))


    Ixf = lapply(1:M, function(m){
      matrix(apply(do.call(cbind,
                           lapply(1:length(perid.x), function(jjj) apply(perall[[jjj]]/lspec.ests[m,jjj] + log(lspec.ests[m,jjj]),1,mean))),1,mean),
             M,M)
    })

    ## Whittle likelihood
    wlik2 = -1 * sum(eigenval) - 1 * mean(log(lspec.ests))+ log(abs(det(Wnew)))

    lim = err
    print(paste("cICA-tf - Iteration ", iter, ": error is equal to ",lim, sep = ""))


    if (wlik < wlik2 ) {
      #       print('Whittle likelihood increased.')
      wlik = wlik2
      Wtmp = Wnew
    }else{print(paste("cICA-tf - Iteration ", iter,
                      ": current Whittle likelihood(", wlik2, ") is smaller than previous one (",
                      wlik, ")."))
      Wtmp = Wold ## Wtmp is the temporal saving of the local maxima
      Ndec = Ndec+1 ## we will allow 3 consecutive decrease in the whittle likelihood.
    }

    if (( (iter == maxit ) & NInv < nmaxit)) {
      print("Color ICA: iteration reaches to maximum. Start new iteration.")
      Wnew = qr.Q(qr(matrix(rnorm(M * M), M, M)))
      iter = 0
      Ndec = 0
      NInv = NInv + 1
      lim=1
    }

    if (NInv == nmaxit) {
      print("Color ICA: no convergence")
    }

    ## update old W to new
  }

  wt = Wtmp %*% K
  result = list()
  result$W = Wtmp
  result$K = K
  result$A = solve(wt)#t(wt) %*% solve(wt %*% t(wt))
  result$S = wt %*% Xin
  result$X = Xin
  result$iter = iter
  result$NInv = NInv
  return(result)

}
