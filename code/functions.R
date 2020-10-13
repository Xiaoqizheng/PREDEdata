rm(list = ls())

## library packages
library(matrixStats)
library(quadprog)
library(e1071)
library(corrplot)
library(RColorBrewer)
library(gplots)
library(gtools)
library(vioplot)

## reference-based 
## Y: the methylation profiles of tumor samples
## W: the methylation profiles of reference cell lines
## output: the proportions of constituent cell types in each input sample

RB <- function(Y,W){ 
  ## this function is the same for both GE and ME data
  H.pred = t(QPfunction(Y = Y, Xmat = W, sumLessThanOne=TRUE, nonNeg=!sumLessThanOne))
  return(list(H = H.pred))
}

## reference-free
## Y: the methylation profiles of tumor samples
## W: the methylation profiles of cell types
## K: the number of cell types
## type: "GE" for gene; "ME" for methylation
## output: the methylation profiles of cell types and the proportions of cell types in each input sample

RF <- function(Y,W = NULL,K,type = "GE",iters = 500,rssDiffStop=1e-10){
  mu0 <- RefFreeCellMixInitialize(Y,K=K,method="ward.D2") ## initial profiles of cell types
  mu0 <- matrix(mu0,ncol = K)
  rss0 = 0
  for(i in 1:iters){
    flag <- !apply(is.na(mu0),1,any)
    omega <- QPfunction(Y[flag,],mu0[flag,],sumLessThanOne=TRUE)
    if (type == "GE"){
      mu <- QPfunction(t(Y), omega, sumLessThanOne=FALSE,W_lessThanOne = F)
    } else if (type == "ME"){
      mu <- QPfunction(t(Y), omega, sumLessThanOne=FALSE,W_lessThanOne = T)
    } else {
      stop("Please input the correct data type!")
    }
    
    rss.new = norm(Y - mu%*%t(omega),"F")^2
    ## check convergence
    dd = abs(rss.new-rss0)
    if(dd < rssDiffStop)
      break
    
    ## updata
    mu0 = mu; rss0 = rss.new
  }
  
  W.pred = mu
  H.pred = t(omega)
  
  # AIC 
  rss = norm(Y - W.pred %*% H.pred,type = "F")^2
  nSamples = ncol(Y)*nrow(Y) ### sample size * number of CpG
  nParam = K*(nrow(Y)+ncol(Y)) ### total number of parameters
  aic = nSamples*log(rss/nSamples)+ 2*nParam + (2*nParam*(nParam+1))/(nSamples-nParam-1) 
  
  if(is.null(W)){
    out = list(W=W.pred, H = H.pred, aic = aic)
  } else {
    # Adjust W.pred and H.pred using W
    out = adjustWH(W,W.pred,H.pred)
    out$aic = aic
  }
  out
}



## partial reference-based 
## Y: the methylation profiles of tumor samples
## W: the methylation profiles of cell types
## W1: the methylation profiles of partial reference cell lines
## type: "GE" for gene; "ME" for methylation
## K: the number of cell types

PR <- function(Y,W = NULL,W1,type = "GE",K,iters = 500,rssDiffStop=1e-10){
  if (is.null(W1)){
    ### ref-free
    print("This is ref-free case")
    out.RF = RF(Y,W = W, K = K,type = type)
    W.pred = out.RF$W
    H.pred = out.RF$H
    
    if(is.null(W)){ 
      ## output result directly
      out = list(W=W.pred, H = H.pred)
    } else {
      # Adjust W.pred and H.pred using W
      out = adjustWH(W,W.pred,H.pred)
    }
    return(out)
  } else {
    if (ncol(W1) == K){
      #### ref-based
      print("This is ref-based case")
      return(RB(Y,W))
    } else if (ncol(W1) < K){
      ### partial ref
      mu02 <- RefFreeCellMixInitialize(Y,K=K-ncol(W1),method="ward.D2")
      mu02 <- as.matrix(mu02)
      mu0 <-cbind(as.matrix(W1),mu02)
      rss0 = 0
      
      for(i in 1:iters){
        flag <- !apply(is.na(mu0),1,any)
        omega <- QPfunction(Y[flag,],mu0[flag,])
        omega1 <- omega[,1:ncol(W1)]
        omega2 <- omega[,(ncol(W1)+1):ncol(mu0)]

        if (type == "GE"){
          mu1 <- QPfunction(t(Y-W1%*%t(omega1)), omega2, sumLessThanOne=FALSE)
        } else if (type == "ME"){
          mu1 <- QPfunction(t(Y-W1%*%t(omega1)), omega2, sumLessThanOne=FALSE,W_lessThanOne = T)
        } else {
          stop("Please input the correct data type!")
        }
        mu = cbind(W1,mu1)
        rss.new = norm(Y - mu%*%t(omega),"F")^2
        ## check convergence
        dd = abs(rss.new-rss0)
        if(dd < rssDiffStop)
          break
        
        ## updata
        mu0 = cbind(W1,mu1); rss0 = rss.new
      }
      
      W.pred = mu
      H.pred = t(omega)
      
      # AIC
      rss = norm(Y - W.pred %*% H.pred,type = "F")^2
      nSamples = ncol(Y)*nrow(Y) ### sample size * number of CpG
      nParam = K*(nrow(Y)+ncol(Y)) - nrow(W1)*ncol(W1) ### total number of parameters
      aic = nSamples*log(rss/nSamples)+ 2*nParam + (2*nParam*(nParam+1))/(nSamples-nParam-1) 
      
      if(is.null(W)){
        out = list(W=W.pred, H = H.pred, aic = aic)
      } else {
        # Adjust W.pred and H.pred using W
        out = adjustWH(W,W.pred,H.pred)
        out$aic = aic
      }

      return(out)
    } else {
      stop("ncol(W1) should be less than or equal to ncol(W)!")
    }
  }
}



QPfunction <- function(Y, Xmat, sumLessThanOne=TRUE, nonNeg=!sumLessThanOne, W_lessThanOne = F){
  Xmat = as.matrix(Xmat)
  nCol = dim(Xmat)[2]   # the number of cell types
  nSubj = dim(Y)[2]     # sample size
  mixCoef = matrix(0, nSubj, nCol) ## the proportion of cell types in each samples
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(Xmat)
  
  # calculate H matrix from W
  if(sumLessThanOne){
    Amat = cbind(rep(-1,nCol), diag(nCol))
    b0vec = c(-1,rep(0,nCol))
    
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,]) %*% Xmat[obs,]
      dvec = t(Xmat[obs,])%*%Y[obs,i]
      meq = 1
      mixCoef[i,] = solve.QP(Dmat=Dmat, dvec=dvec, Amat=Amat, bvec=b0vec, meq=meq)$sol
    }
  }
  
  
  # calculate W matrix from H
  else if(nonNeg){  
    
    if(W_lessThanOne == F){ ## only x >= 0, for gene expression
      Amat = cbind(diag(nCol))
      b0vec = c(rep(0,nCol))
    } else if(W_lessThanOne == T){ ## 0 <= x <= 1, for DNA methylation
      Amat = cbind(-diag(nCol), diag(nCol))
      b0vec = c(rep(-1,nCol),rep(0,nCol))
    }
    
    for(i in 1:nSubj){
      obs = which(!is.na(Y[,i]))
      Dmat = t(Xmat[obs,])%*%Xmat[obs,]
      mixCoef[i,] = solve.QP(Dmat, t(Xmat[obs,])%*%Y[obs,i], Amat, b0vec)$sol
    }
  }
  
  return(mixCoef)
}



adjustWH <- function(W,W.pred,H.pred){
  MAT = cor(W,W.pred)
  n = min(dim(MAT))
  a.index = 1:n
  
  for(i in 1:n) { 
    z = which(MAT == max(MAT),arr.ind=T)
    a.index[z[1,1]] = z[1,2] 
    MAT[,z[1,2]] = -1
    MAT[z[1,1],] = -1 
  }
  
  W.adj = W.pred[,a.index]
  H.adj = H.pred[a.index,]
  return(list(W = W.adj,H = H.adj))
}

RefFreeCellMixInitialize <- function(Y,K=2,Y.Distance=NULL, Y.Cluster=NULL, 
                                     largeOK=FALSE, dist.method = "euclidean", ...){
  
  if(!is.matrix(Y) | !is.numeric(Y)){
    stop("Y is not a numeric matrix\n")
  }
  n <- dim(Y)[2]
  
  if(is.null(Y.Cluster)){
    if(is.null(Y.Distance)){
      if(n>2500 & !largeOK){
        stop("Y has a large number of subjects!  If this is what you really want, change 'largeOK' to TRUE\n")
      }
      Y.Distance <- dist(t(Y),method=dist.method)
    }
    Y.Cluster <- hclust(Y.Distance,...)
  } 
  
  classes <- cutree(Y.Cluster, K)
  s <- split(1:n,classes)
  
  sapply(s, function(u) apply(Y[,u,drop=FALSE],1,mean,na.rm=TRUE))
}

###obtain M,H,W
generate_bulk <- function(celllines,nSample =20,csd = 0.1){
  k = ncol(celllines)
  m = nrow(celllines)
  ## generate composition matrix H
  
  H.dat = matrix(runif(k*nSample,min = 0,max = 1),nrow = k)
  H = apply(H.dat,2, function(x) x/sum(x))
  
  Y = as.matrix(celllines) %*% H
  ## add noise
  sd = as.vector(Y)*csd 
  noise = matrix(rnorm(n = m*nSample,mean = 0,sd = sd),nrow = m)
  Y = Y + noise
  
  Y[Y<0] = 0
  
  return(list(Y = Y,H = H,W = celllines))
}


DoCBS <- function(Y, W, nu.v = c(0.25, 0.5, 0.75)) {
  require(e1071)
  
  est.lm <- list()
  
  nui <- 1
  
  for (nu in nu.v) {
    est.m <- matrix(nrow = ncol(Y), ncol = ncol(W))
    
    for (s in 1:ncol(Y)) {
      svm.o <-
        svm(
          x = W,
          y = Y[, s],
          scale = TRUE,
          type = "nu-regression",
          kernel = "linear",
          nu = nu
        )
      
      coef.v <- t(svm.o$coefs) %*% svm.o$SV
      
      coef.v[which(coef.v < 0)] <- 0
      
      total <- sum(coef.v)
      
      coef.v <- coef.v / total
      
      est.m[s, ] <- coef.v
      
    }
    est.lm[[nui]] <- est.m
    
    nui <- nui + 1
    
  }
  
  H = matrix(nrow = ncol(Y), ncol = ncol(W))
  #### select best nu
  rmse.m <- matrix(NA, nrow = ncol(Y), ncol = length(nu.v))
  
  for (nui in 1:length(nu.v)) {
    reconst.m <- W %*% t(est.lm[[nui]])
    
    for (s in 1:ncol(Y)) {
      rmse.m[s, nui] <- sqrt(mean((Y[, s] - reconst.m[, s]) ^ 2))
      
    }
  }
  nu.idx <- apply(rmse.m, 1, which.min)
  
  H <- est.m
  
  for (s in 1:ncol(Y)) {
    H[s, ] <- est.lm[[nu.idx[s]]][s, ]
    
  }
  return(H)
  
}

### simulation function 
simuAll <- function(W,W1.index = 1:4,csd = 0.1,nSample = 100,type = "GE",fsmethod = "cv",method = "all",nrep = 20){
  result = list()
  
  W1 = W[,W1.index]
  W2.index = setdiff(1:ncol(W),W1.index)
  
  if(method == "all" | method == "RB"){
    cat("Running RB ...\n")
    RBout = c()
    for (rep in 1:nrep){ ## repeat 20 times
      cat(rep,"\t")
      Bulk= generate_bulk(celllines = W,nSample = nSample,csd = csd) ## generate 20 samples 
      
      ### select feature 
      feat = select_feature(mat = Bulk$Y,method = fsmethod,nmarker = 1000,startn = 0)
      Bulk$Y = Bulk$Y[feat,]
      W1 = W[feat,W1.index]
      
      ### run deconvolution
      H1 = Bulk$H[W1.index,]
      out <- RB(Y = Bulk$Y, W = W1)
      H_AbsBias = mean(colSums(abs(H1 - out$H)) / nrow(H1))
      H_corr = mean(diag(cor(t(H1),t(out$H))))
      RBout = rbind(RBout,c(H_AbsBias,H_corr))
    }
    colnames(RBout) = c("H_AbsBias","H_corr")
    result$RBout = RBout
  }
  
  if (method == "all" | method == "RF"){
    cat("\nRunning RF ...\n")
    RFout = c()
    for (rep in 1:nrep){ ## repeat 20 times
      cat(rep,"\t")
      Bulk= generate_bulk(celllines = W,nSample = nSample,csd = csd) ## generate 20 samples
      
      ### select feature 
      feat = select_feature(mat = Bulk$Y,method = fsmethod,nmarker = 1000,startn = 0)
      Bulk$Y = Bulk$Y[feat,]
      W.new = W[feat,]
      
      ### run deconvolution
      out <- RF(Y = Bulk$Y, W = W.new,type = type,K = ncol(W.new),iters = 100)
      H_AbsBias = mean(colSums(abs(Bulk$H - out$H)) / nrow(Bulk$H))
      H_corr = mean(diag(cor(t(Bulk$H),t(out$H))))
      W_corr = mean(diag(cor(W.new,out$W))) ### all W
      RFout = rbind(RFout,c(H_AbsBias,H_corr,W_corr))
    }
    colnames(RFout) = c("H_AbsBias","H_corr","W_corr")
    result$RFout = RFout
  }
  if(method == "all" | method == "PR"){
    cat("\nRunning PR ...\n")
    PRout = c()
    for (rep in 1:nrep){ ## repeat 20 times
      cat(rep,"\t")
      #browser()
      Bulk= generate_bulk(celllines = W,nSample = nSample,csd = csd) ## generate 20 samples 
      
      ### select feature 
      feat = select_feature(mat = Bulk$Y,method = fsmethod,nmarker = 1000,startn = 0)
      Bulk$Y = Bulk$Y[feat,]
      W.new = W[feat,]
      W1 = W[feat,W1.index]
      
      ### run deconvolution
      
      out <- PR(Y = Bulk$Y, W = as.matrix(W.new), W1 = as.matrix(W1),type = type,K = ncol(W.new),iters = 100)
      H_AbsBias = mean(colSums(abs(Bulk$H - out$H)) / nrow(Bulk$H))
      H_corr = mean(diag(cor(t(Bulk$H),t(out$H))))
      W_corr = mean(diag(as.matrix(cor(W.new[,W2.index],out$W[,W2.index])))) ## only W2
      PRout = rbind(PRout,c(H_AbsBias,H_corr,W_corr))
    }
    colnames(PRout) = c("H_AbsBias","H_corr","W_corr")
    result$PRout = PRout
  }
  if(method == "all" | method == "CBS"){
    cat("\nRunning CBS ...\n")
    CBSout = c()
    for (rep in 1:nrep){ ## repeat 20 times
      cat(rep,"\t")
      Bulk= generate_bulk(celllines = W,nSample = nSample,csd = csd) ## generate 20 samples 
      
      ### select feature 
      feat = select_feature(mat = Bulk$Y,method = fsmethod,nmarker = 1000,startn = 0)
      Bulk$Y = Bulk$Y[feat,]
      W1 = W[feat,W1.index]
      ### run deconvolution
      
      H1 = Bulk$H[W1.index,]
      out <- DoCBS(Y=Bulk$Y, W=W1) 
      H_AbsBias = mean(colSums(abs(H1 - t(out))) / nrow(H1))
      H_corr = mean(diag(cor(t(H1),out)))
      CBSout = rbind(CBSout,c(H_AbsBias,H_corr))
    }
    colnames(CBSout) = c("H_AbsBias","H_corr")
    result$CBSout = CBSout
  }
  return(result)
}

### feature selection 

select_feature <- function(mat,method = "cv",nmarker=1000, startn=0){

  if(method == "topvar"){ ## sort sd
    Var = rowVars(mat)
    difgene = order(Var,decreasing = T)[startn + 1:nmarker]     #choose the top 1000 features
    index.select = rownames(mat)[difgene]
  }
  if(method == "random"){ ## random
    index.select = sample(rownames(mat),size = nmarker)
  }
  if(method == "cv"){ ## sort sd/mean 
    mm = rowMeans(mat)
    vv = rowVars(mat)
    cv = sqrt(vv) / (mm + 1)
    #cv = sqrt(vv)
    cv[is.na(cv)] = 0
    ix = sort(cv, dec=TRUE, index=TRUE)$ix
    index.select = rownames(mat)[ix[startn + 1:nmarker]]
  }
  return(index.select)
}




