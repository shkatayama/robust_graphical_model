library(mvtnorm)
library(robustbase)
library(abind)
library(magic)
library(ICSNP)
library(QUIC)
library(igraph)
library(huge)
library(clime)

##covariance matrices
ar1 <- function(p,rho){
	si <- matrix(1,p,p)
	for(j in 2:p){
		for(i in 1:(j-1)){
			si[i,j] <- rho^(abs(i-j))
			si[j,i] <- si[i,j]
		}
	}
	return(si)
}

