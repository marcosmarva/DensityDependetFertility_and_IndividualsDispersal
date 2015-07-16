## We plot nu, the global reproductive number, as function of 
## dispersal rates muJ y muA, in two cases
##    symmetric: muJ = muA
##    asymmetric: muA = 1-muJ
## for large values of the intrinsic reproductive number phi_i(0,0)

rm(list=ls())
library(Matrix)
library(expm)

transientIterations = 500
plottedIterations = 200
iterations_per_value = transientIterations + plottedIterations



sigmaJ1 = .4
sigmaJ2 = .6 
sigmaA1 = .4
sigmaA2 = .3



transient = 500
ttail=200
iterations_per_value = transient + ttail


mu=seq(from = 0, to = 1, by = 0.01)
a1=.6
a2=.65
StrongDensityDependence=TRUE
if(StrongDensityDependence){
  phi1 = function(J, A, muJ, muA){
    psi1 * exp(- a1 * ( J * muJ +  A * muA)) 
  }
  
  phi2 = function(J, A, muJ, muA){
    psi2 * exp(-a2 * ( J * (1 - muJ) +  A * (1 - muA))) 
  }  
} else {
  phi1 = function(J, A, muJ, muA){
    psi1 / (1 + J * muJ +  A * muA) 
  }
  phi2 = function(J, A, muJ, muA){
    psi2 / (1 + J * (1 - muJ1) +  A * (1 - muA1)) 
  }  
}
barPhi = function(J, A, muJ1, muA1){
  muA1 * phi1(J, A, muJ1, muA1) + 
    (1 - muA1) * phi2(J, A, muJ1, muA1)    
}

barSigmaJ = function(muJ1)
  sigmaJ1 * muJ1 + sigmaJ2 * (1 - muJ1)

barSigmaA = function(muA1)
  muA1 * sigmaA1 + (1 - muA1) * sigmaA2

RedSys = function(J, A, muJ1, muA1){
  return(c(A * barPhi(J, A, muJ1, muA1),  J * barSigmaJ(muJ1) + A * barSigmaA(muA1)))
}


psi1 = 5
alpha = 1.1
psi2 = alpha * psi1


## First
mJ1 = matrix(rep(0, length(mu) * iterations_per_value), ncol= iterations_per_value)
mA1 = matrix(rep(0, length(mu) * iterations_per_value), ncol= iterations_per_value)
J0 = 5
A0 = 5
mA1[ , 1] = J0
mJ1[ , 1] = A0
loopIter = 2:(iterations_per_value)
for(r in 1:length(mu)){
  for(c in loopIter){  
    nextIter = RedSys(mJ1[r, c-1], mA1[r, c-1], mu[r], mu[r])  
    mJ1[r, c] = nextIter[1]
    mA1[r, c] = nextIter[2]
  }  
}


mJ2 = matrix(rep(0, length(mu) * iterations_per_value), ncol= iterations_per_value)
mA2 = matrix(rep(0, length(mu) * iterations_per_value), ncol= iterations_per_value)
mA2[ , 1] = J0
mJ2[ , 1] = A0

loopIter = 2:(iterations_per_value)
for(r in 1:length(mu)){
  for(c in loopIter){  
    nextIter = RedSys(mJ2[r, c-1], mA2[r, c-1], mu[r], 1-mu[r])  
    mJ2[r, c] = nextIter[1]
    mA2[r, c] = nextIter[2]
  }  
}

iterations_per_value=200

mJaux1 = matrix(rep(0, iterations_per_value), ncol= iterations_per_value)
mAaux1 = matrix(rep(0, iterations_per_value), ncol= iterations_per_value)
mJaux1[1 , 1] = J0
mAaux1[1 , 1] = A0

loopIter = 2:(iterations_per_value)
for(c in loopIter){  
  nextIter = RedSys(mJaux1[1, c-1], mAaux1[1, c-1], 1, 1)  
  mJaux1[1, c] = nextIter[1]
  mAaux1[1, c] = nextIter[2]
}  

mJaux2 = matrix(rep(0, iterations_per_value), ncol= iterations_per_value)
mAaux2 = matrix(rep(0, iterations_per_value), ncol= iterations_per_value)
mJaux2[1 , 1] = J0
mAaux2[1 , 1] = A0
loopIter = 2:(iterations_per_value)
for(c in loopIter){  
  nextIter = RedSys(mJaux2[1, c-1], mAaux2[1, c-1], 0, 0)  
  mJaux2[1, c] = nextIter[1]
  mAaux2[1, c] = nextIter[2]
}


## Total population size
# Symmetrical dispersal
mN2aux1 = mJaux1 + mAaux1
# Asymmetrical dispersal
mN2aux2 = mJaux2 + mAaux2

## PLOT FIGURES

par(mfrow = c(1,2))


# First bifurcation diagram
mN1= mJ1[, -(1:transient)] + mA1[, -(1:transient)]

matplot(mu, mN1, type="l", col="black", xlab="mu", ylab = "Population size", 
        xlim=c(0,1), ylim=c(0,5.5) )
# Include also the asymptotic local population at patch 1 size if isolated patches
par(new=TRUE)
matplot(mu, mN1*mu, type="l", col="black", xlab="mu", ylab = "Population size", 
        xlim=c(0,1), ylim=c(0,5.5), lty=2 )
par(new=TRUE)
matplot(mu, mN1*(1-mu), type="l", col="black", xlab="mu", ylab = "Population size", 
        xlim=c(0,1), ylim=c(0,5.5), lty=3 )
par(new=TRUE)
matplot(mu, rep(x = mN2aux1[1,iterations_per_value], length(mu)), type="l",  col="darkgray", 
        xlab="mu", ylab = "Population size",  lwd = 2, 
        xlim=c(0,1), ylim=c(0,5.5),lty=2 )
# Include also the asymptotic local population at patch 2 it isolated patches
par(new=TRUE)
matplot(mu, rep(x = mN2aux2[1,iterations_per_value], length(mu)), type="l",  col="darkgray", 
        xlab="mu", ylab = "Population size",  lwd = 2, 
        xlim=c(0,1), ylim=c(0,5.5),lty=3 )
# Include also the asymptotic total population size if isolated patch
par(new=TRUE)
matplot(mu, rep(x = mN2aux1[1,iterations_per_value] + mN2aux2[1,iterations_per_value], length(mu)), type="l",  col="darkgray", 
        xlab="mu", ylab = "Population size",  lwd = 2, 
        xlim=c(0,1), ylim=c(0,5.5),lty=4 )





mN2 = mJ2[, -(1:transient)] + mA2[, -(1:transient)]
matplot(mu, mN2, type="l", col="black", xlab="mu", ylab = "Population size", 
        xlim=c(0,1), ylim=c(0,5.5) )
par(new=TRUE)
matplot(mu, mN2*mu, type="l", col="black", xlab="mu", ylab = "Population size", 
        xlim=c(0,1), ylim=c(0,5.5), lty=2 )
par(new=TRUE)
matplot(mu, mN2*(1-mu), type="l", col="black", xlab="mu", ylab = "Population size", 
        xlim=c(0,1), ylim=c(0,5.5), lty=3 )
# Include also the asymptotic local population at patch 1 size if isolated patches
par(new=TRUE)
matplot(mu, rep(x = mN2aux1[1,iterations_per_value], length(mu)), type="l",  col="darkgray", 
        xlab="mu", ylab = "Population size",  lwd = 2, 
        xlim=c(0,1), ylim=c(0,5.5),lty=2 )
# Include also the asymptotic local population at patch 2 it isolated patches
par(new=TRUE)
matplot(mu, rep(x = mN2aux2[1,iterations_per_value], length(mu)), type="l",  col="darkgray", 
        xlab="mu", ylab = "Population size",  lwd = 2, 
        xlim=c(0,1), ylim=c(0,5.5),lty=3 )
# Include also the asymptotic total population size if isolated patch
par(new=TRUE)
matplot(mu, rep(x = mN2aux1[1,iterations_per_value] + mN2aux2[1,iterations_per_value], length(mu)), type="l",  col="darkgray", 
        xlab="mu", ylab = "Population size",  lwd = 2, 
        xlim=c(0,1), ylim=c(0,5.5),lty=4 )

# 
# matplot(1:iterations_per_value, mN2aux1[1,], type="l",  
#         col="darkgray", lwd = 2, xlab="Time", ylab = "Population size", ylim=c(0,10))
# 
#   
# matplot(1:iterations_per_value, mN2aux2[1,], type="l",  
#         col="darkgray", lwd = 2, xlab="Time", ylab = "Population size", ylim=c(0,10))
