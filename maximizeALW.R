################################################################################
# Author:         Sophie J. Kersting
# Last update:    Nov 2022
# R-version:      4.2.2
################################################################################
#-------------------------------------------------------------------------------
# Maximizing edge imbalance integrals
#-------------------------------------------------------------------------------
# Supplementary material of the manuscript
# "Measuring 3D tree imbalance of graph-theoretical plant models"
# by Sophie Kersting, Luise Kuehn and Mareike Fischer
#-------------------------------------------------------------------------------
################################################################################
# Please first read the corresponding explanations in the manuscript. In this 
# script we use A to refer to A(v). We maximize the edge imbalance integral of
# the inner edge of a two edge path graph for the three parameters A, L and W.
################################################################################

#-------------------------------------------------------------------------------
# 1.) Implementations of all simplified formulas (examples for usage can be 
# found below):
#--- Centroid of pendant tree of edge subdividing node "xtilde"
centr_xtilde <- function(A,L,W,x){
  cent <- (1/(x+W*L)) * (W*L* round(matrix(c(cos(A), -sin(A),
                                          sin(A),  cos(A) ),
                                        byrow = T, nrow = 2), 10) %*% 
                        c(0,-L/2) + x * c(0,x/2))
  return(as.vector(cent))
}
#--- Imbalance values of "xtilde"
A_xtilde <- function(A,L,W,x){
  cent <- round(centr_xtilde(A,L,W,x),10)
  if(cent[1]==0 && cent[2]==x){
    return(0)
  }else{
    A_xtilde <- treeDbalance::angle3dVec(a=c(0,-1,0), b=c(cent,0)-c(0,x,0))
    return(A_xtilde)
  }
}
alpha_xtilde <- function(A,L,W,x){
  cent <- round(centr_xtilde(A,L,W,x),10)
  if(cent[1]==0 && cent[2]==x){
    return(0)
  }else{
    A_xtilde <- treeDbalance::angle3dVec(a=c(0,-1,0), b=c(cent,0)-c(0,x,0))
    if(A_xtilde<=pi/2){
      return(A_xtilde)
    }else{
      return(pi-A_xtilde)
    }
  }
}
#--- Edge imbalance integrals of e_v
A_edge_imbal <- function(A,L,W){
  eimbal <- stats::integrate(f = function(x){
    sapply(x, function(y){A_xtilde(A,L,W,y)})},
    lower = 0, upper = 1)
  return(eimbal$value)
}
alpha_edge_imbal <- function(A,L,W){
  eimbal <- stats::integrate(f = function(x){
    sapply(x, function(y){alpha_xtilde(A,L,W,y)})},
    lower = 0, upper = 1)
  return(eimbal$value)
}
#--- Imbalance indices for the complete two edge path graph T
A_index_2EPG <- function(A,L,W, weight = c("w","l")){
  e_imbal <- A_edge_imbal(A,L,W)
  if(weight=="l"){
    return(e_imbal/(1+L))
  }else if(weight=="w"){
    return(e_imbal/(1+W*L))
  }else{
    stop("Unsuitable weighting method.\n")
  }
}
alpha_index_2EPG <- function(A,L,W, weight = c("w","l")){
  e_imbal <- alpha_edge_imbal(A,L,W)
  if(weight=="l"){
    return(e_imbal/(1+L))
  }else if(weight=="w"){
    return(e_imbal/(1+W*L))
  }else{
    stop("Unsuitable weighting method.\n")
  }
}
if(FALSE){ # Examples for the usage of these functions:
  #A=pi/2; L=1; W=1; x=0.5
  centr_xtilde(A=pi,L=2,W=2,x=0.2) # [1] 0.0000000 0.9571429
  A_xtilde(A=pi/1.2,L=2,W=2,x=0.2) # [1] 2.494008
  alpha_xtilde(A=pi/1.2,L=2,W=2,x=0.2) # [1] 0.6475849
  A_edge_imbal(A=pi,L=2.4142,W=1) # [1] 3.141593
  A_edge_imbal(A=pi,L=sqrt(1/1+1)+1,W=100) # [1] 3.141593
  alpha_edge_imbal(A=1.570926,L=10000,W=1000) # [1] 1.570742
  A_index_2EPG(A=1.5*pi/2,L=2,W=2, weight = "w") # [1] 0.3536337
  alpha_index_2EPG(A=1.5*pi/2,L=2,W=2, weight = "w") # [1] 0.2362989
}

#-------------------------------------------------------------------------------
# 2.) Maximize the edge imbalance functions over all three parameters:
func2optA <- function(z){return(A_edge_imbal(A=z[1], L=z[2], W=z[3]))}
r_A <- stats::optim(c(pi/2, 1, 1), # starting values
                    f = func2optA, # function to optimize
                    control = list(fnscale=-1), # maximize (not minimize)
                    method = "L-BFGS-B", # method (>1 parameters with bounds)
                    lower = c(0,0.001,0.001), # lower bounds
                    upper = c(pi,10000,10000)) # upper bounds
func2optalpha <- function(z){return(alpha_edge_imbal(A=z[1], L=z[2], W=z[3]))}
r_alpha <- stats::optim(c(pi/2, 1, 1), 
                        f = func2optalpha, 
                        control = list(fnscale=-1),
                        method = "L-BFGS-B", 
                        lower = c(0,0.001,0.001), 
                        upper = c(pi,10000,10000))
# Observation: Both edge imbalance integrals get very close to their supremum:
pi-r_A$value # for A=3.141593, L=3.652100, W=1.458276 (r_A$par)
pi/2-r_alpha$value # for A=1.570945, L=4526.482396, W=917.099081 (r_alpha$par)
# For 3D imbalance measure A, in contrast to alpha, we need an angle of 180°
# with comparably small L and W. However, for alpha the angle maximizing angle
# is slightly larger than a right angle (pi/2=1.570796) and L and W are 
# significantly larger.

#-------------------------------------------------------------------------------
# 3.) Maximize the edge imbalance functions over two parameters (length fixed):
func2optA <- function(z){return(A_edge_imbal(A=z[1], L=0.5, W=z[2]))}
r_A <- stats::optim(c(pi/2, 1), # starting values
                    f = func2optA, # function to optimize
                    control = list(fnscale=-1), # maximize (not minimize)
                    method = "L-BFGS-B", # method (>1 parameters with bounds)
                    lower = c(0,0.001), # lower bounds
                    upper = c(pi,10000)) # upper bounds
func2optalpha <- function(z){return(alpha_edge_imbal(A=z[1], L=0.5, W=z[2]))}
r_alpha <- stats::optim(c(pi/2, 1), 
                        f = func2optalpha, 
                        control = list(fnscale=-1),
                        method = "L-BFGS-B", 
                        lower = c(0,0.001), 
                        upper = c(pi,10000))

#-------------------------------------------------------------------------------
# 4.) Maximize the edge imbalance functions over one parameter:
# set W=1, optimize the angle A for given L in [0.01,3]
ws <- 2 # try e.g. W=0.5 and W=2 here
ls <- seq(0.01,3, len=100) # lengths
Avs_A <- rep(NA, length(ls)) # maximizing angle A(v) for A-measurement
Avs_alpha <- rep(NA, length(ls)) # maximizing angle A(v) for alpha-measurement
A_edges <- rep(NA, length(ls)) # respective edge imbalance A_e
alpha_edges <- rep(NA, length(ls)) # respective edge imbalance alpha_e
for(i in 1:length(ls)){
  # search maximum for A-measurement
  res <- stats::optimize(f=function(y){return(A_edge_imbal(A=y, L=ls[i], 
                                                           W=ws))}, 
                         maximum = T, interval = c(0,pi))
  Avs_A[i] <- res$maximum
  A_edges[i] <- res$objective
  # search maximum for alpha-measurement
  res <- stats::optimize(f=function(y){return(alpha_edge_imbal(A=y, L=ls[i], 
                                                           W=ws))}, 
                         maximum = T, interval = c(0,pi))
  Avs_alpha[i] <- res$maximum
  alpha_edges[i] <- res$objective
}
# Gather all data in a data frame
mydf <- data.frame(W=rep(1,4*length(ls)), L=c(ls,ls,ls,ls), 
                   angles=c(Avs_A,A_edges,Avs_alpha,alpha_edges),
                   AorI=c(rep("Av",length(ls)),rep("edge_imbal",length(ls)),
                             rep("Av",length(ls)),rep("edge_imbal",length(ls))),
                   measurement = c(rep("A",2*length(ls)),
                                   rep("alpha",2*length(ls))))
mydf$measurement <- as.factor(mydf$measurement)
mydf$AorI <- as.factor(mydf$AorI)
mydf$size <- mydf$AorI
levels(mydf$size) <- c(0.6,0.4)
mydf$size <- as.numeric(as.vector(mydf$size))
summary(mydf)

#-------------------------------------------------------------------------------
# 5.) Plot the results of 4.):
plot(ls,Avs_A, ylim = c(0,3.2), pch= 16, col = "black",
     xlab = "Length L", ylab = "A(v) and edge imbalance",
     main = "Maximizing A(v) and maximal edge imbalance for W=2 and given L")
points(ls,A_edges, pch=17, col = "black")
points(ls,Avs_alpha, pch=16, col = "gray")
points(ls,alpha_edges, pch=17, col = "gray")
legend("bottomright",
       legend = c("A(v) for A","A(v) for alpha",
                  "edge imbalance for A","edge imbalance for alpha"), 
       pch = c(16,16,17,17), col = c("black", "gray", "black", "gray"))


# local maximum of A(v) for alpha
ls[which.max(Avs_alpha)]

