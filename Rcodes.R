# Rcodes for calculating robust covariance matrix and presicion matrix via gamma divergence

#####################gam.cov############################
# The function "gam.cov" is for calculating robust covariance matrix via gamma divergence.
# Need to define the function "gam.meansd", "gam.cov" and "obj.rho" in advance.
# Requires
#     fdata:    Input matrix of size n (samples) times p (variables).
#     gamma:    Tuning parameter of gamma divergence. Default 0.3.
#     method:   Optimization algorithm to be employed. Default "projection"
#               "projection" is corresponding to the projected gradient in the article.
#               If one wants to use "optimize" in R, select "optim".


#####################gam.graph###########################
# The function "gam.graph" is to output the graph with the tuning parameter selected from data.
# This requires QUIC, huge, clime packages
# Need to define "gam.cov" and "pd" in advance
# Requires
#     precision:      Whether "glasso" (defalt), "nodewise" or "clime" shoul be employed.
#     nlambda:        Number of candidate lambda. Default 10.
#     lambda.ratio:   This value determines min(lambda) as lambda.ratio * max(lambda).
#                     The max(lambda) is automatically generated from the data.
#     node.beta:      Only applicable when method="nodewise". The cut off value in StARS selection. Default 0.2.
#     repeats:        Number of replications in 2-fold cross validation or StARS. Default 1.
#                     If one selects "nodewise", this value should be larger than 10.
# Outputs
#     If "glasso" or "clime" is selected, the output is a precision matrix.
#     Otherwise, the output is an adjacency matrix.

# Sample codes
#     fdata <- matrix(rnorm(100*20),100,20)
#     Sigma <- gam.cov(fdata)
#     QUIC(Sigma,rho=0.3,msg=0)$X
#     gam.graph(fdata)



###################################################################
library(QUIC)
library(huge)
library(clime)

#gamma divergence for mean and standard deviation
#input is n * 1 vector
gam.meansd <- function(vector,gamma){ ## m_0 and s_0 is a constant
  count <- 0; epsilon <- 10
  m_0 <- median(vector)
  s_0 <- (qnorm(0.75,0,1)^(-1))*median(abs(vector-m_0))
  while(count <= 500 && epsilon >= 0.0001){
    z <- vector - m_0
    w <- exp((-gamma/(2*s_0))*(z^2))/sum(exp((-gamma/(2*s_0))*(z^2)))
    m_1 <- sum(vector*w)
    z_1 <- vector - m_1
    s_1 <- (1+gamma)*sum((z_1^2)*w)
    count <- count + 1
    epsilon <- (max(abs(m_1-m_0)) + max(abs(s_1-s_0)))/2
    m_0 <- m_1; s_0 <- s_1
  }
  return(c(m_0,s_0))
}


#gamma divergence for correlation
#input is n * 2 matrix standardized by gam.meansd
proj.rho <- function(data,gamma,rho0=0,phi0=0.01,R=0.99){
  count <- 0; epsilon <- 100
  rho <- rho0
  phi <- phi0
  x1 <- data[,1]; x2 <- data[,2]
  topt <- 1; opt <- 0
  while(count <= 500 && epsilon >= 0.0001){
    repeat{	
      if(topt > opt){
        tmp1 <- x1^2 + x2^2 - 2*rho*x1*x2
        tmp2 <- gamma*((1+rho^2)*x1*x2 - rho*(x1^2 + x2^2))/((1-rho^2)^2)
        obj <- (-1/gamma)*log(sum(exp((-0.5*gamma*tmp1/(1-rho^2))))) + (0.5*log(1-rho^2))/(1+gamma)
        weight <- exp((-0.5*gamma*tmp1/(1-rho^2)))/(sum(exp((-0.5*gamma*tmp1/(1-rho^2)))))
        grad <- (-1/gamma)*sum(weight * tmp2) -rho/((1+gamma)*(1-rho^2))
        rho_bar <- rho - (1/phi)*grad
        if(rho_bar^2 <= R){rho_hat <- rho_bar} else{rho_hat <- sign(rho_bar)*R}
        tmp3 <- x1^2 + x2^2 - 2*rho_hat*x1*x2
        topt <- (-1/gamma)*log(sum(exp((-0.5*gamma*tmp3/(1-rho_hat^2))))) + (0.5*log(1-rho_hat^2))/(1+gamma)
        opt <- obj + grad * (rho_hat - rho) + (phi/2)*(rho_hat - rho)^2
        phi <- 2*phi
      }
      else break
    }
    epsilon <- abs(rho-rho_hat)
    rho <- rho_hat
    phi <- max(phi0,phi/2)
    count <- count + 1
  }
  return(rho)
}



#objective function to use "optimize" function in R
obj.rho <- function(rho,data,gamma){
  x1 <- data[,1]; x2 <- data[,2]
  tmp1 <- x1^2 + x2^2 - 2*rho*x1*x2
  (-1/gamma)*log(sum(exp((-0.5*gamma*tmp1/(1-rho^2))))) + (0.5*log(1-rho^2))/(1+gamma)
}



#robust covariance matrix via gamma divergence
gam.cov <- function(fdata,gamma=0.3,method="projection"){
  n <- dim(fdata)[1]; p <- dim(fdata)[2]
  cor <- matrix(0,p,p); cov <- matrix(0,p,p)
  tmp <- apply(fdata,2,function(x)gam.meansd(x,gamma))
  means <- tmp[1,]; vars <- tmp[2,]
  zz <- sweep(fdata,2,means,"-"); zz <- sweep(zz,2,sqrt(vars),"/")
  for(j in 1:(p-1)){
    tmp_zz <- as.matrix(zz[,-(1:j)])
    tolist <- lapply(1:ncol(tmp_zz),function(i)tmp_zz[,i])
    tolist2 <- lapply(tolist,function(x)t(rbind(zz[,j],as.vector(x))))
    if(method=="projection"){gam <- lapply(tolist2,function(x)proj.rho(x,gamma))}
    if(method=="optim"){gam <- lapply(tolist2,function(x)optimize(obj.rho,c(-0.99,0.99),data=x,gamma=gamma)$minimum)}
    cor[j,(j+1):p] <- unlist(gam)
    cov[j,(j+1):p] <- sqrt(vars[j])*sqrt(vars[(j+1):p])*cor[j,(j+1):p]
  }
  cov <- cov + t(cov); diag(cov) <- vars 
  return(cov)
}



##positive (semi)definite projection
pd <- function(S,delta=0){
  eigs <- eigen(S); eigs$values[which(eigs$values <= delta)] <- delta
  pdS <- eigs$vectors %*% diag(eigs$values) %*% t(eigs$vectors)
  return(pdS)
}



##The estimated graph via gamma divergence.
gam.graph <- function(fdata, gamma=0.3, method="projection", precision="glasso", 
                           nlambda=10, lambda.ratio=0.05, node.beta=0.2, repeats=1){
  
  p <- dim(fdata)[2]
  cvs <- matrix(0,repeats,nlambda)
  star.set <- matrix(0,nlambda,p*(p-1)/2)
  
  Sigma <- pd(gam.cov(fdata,gamma=gamma,method=method))
  lmax <- max(abs(Sigma - diag(diag(Sigma)))); lmin <- lambda.ratio*lmax
  lam <- exp(seq(log(lmin),log(lmax),length=nlambda))
  
  if(precision=="glasso"){all.path <- lapply(lam,function(x)round(QUIC(Sigma,rho=x,msg=0)$X,5))}
  if(precision=="clime"){all.path <- lapply(lam,function(x)round(clime(Sigma,lambda=x,sigma=T,perturb=T)$Omega[[1]],5))}
  if(precision=="nodewise"){all.path <- huge(x=Sigma,lambda=lam,method="mb",verbose=F)$path}
  
  for(j in 1:repeats){ ##repeated 2 fold
    fdata <- fdata[sample(1:nrow(fdata)),] ## random sort
    foldsize <- floor(nrow(fdata)/2)
    test.fdata <- fdata[1:foldsize, ]
    train.fdata <- fdata[-(1:foldsize),]
    Sigma <- pd(gam.cov(train.fdata,gamma=gamma,method=method))
    test.Sigma <- gam.cov(test.fdata,gamma=gamma,method=method)
    
    if(precision=="glasso"){
      path <- lapply(lam,function(x)round(QUIC(Sigma,rho=x,msg=0)$X,5))
      for(k in 1:nlambda){
        cvs[j,k] <- sum(diag(test.Sigma %*% path[[k]])) - log(det(path[[k]]))
      }
    }
    
    if(precision=="clime"){
      path <- lapply(lam,function(x)round(clime(Sigma,lambda=x,sigma=T,perturb=T)$Omega[[1]],5))
      for(k in 1:nlambda){
        cvs[j,k] <- sum(diag(test.Sigma %*% path[[k]])) - log(det(path[[k]]))
      }
    }

    if(precision=="nodewise"){
      path <- lapply(huge(x=Sigma,lambda=lam,method="mb",verbose=F)$path,function(x)as.matrix(x))
      for(k in 1:nlambda){
        star.set[k,] <- star.set[k,] + abs(sign(path[[k]][upper.tri(path[[k]])]))
      }
    }
  }
  
  if(precision=="glasso"){
    cv.mean <- apply(cvs,2,mean)
    opt.tune <- which.min(cv.mean)
  }
  
  if(precision=="clime"){
    cv.mean <- apply(cvs,2,mean)
    opt.tune <- which.min(cv.mean)
  }
  
  if(precision=="nodewise"){
    zeta <- 2*(star.set/repeats)*(1-(star.set/repeats))
    instab <- apply(zeta,1,sum)/length(zeta[1,])
    tmp <- instab <= node.beta
    select.instab <- which(tmp==T)
    if(length(select.instab)==0){
      opt.tune <- which.min(instab)
    } else{
        opt.tune <- min(select.instab)
      }
  }
  
  return(as.matrix(all.path[[opt.tune]]))
}


