library(dplyr)
library(geepack)
library(readr)
library(igraph)

working <- c('independence','exchangeable','ar1','unstructured','userdefined')
std.errs <- c('san.se',"fij",'jack', 'j1s')

betas0 <- c(0.1, 0.2, 0.2, 0.1)
medias1 <- betas0[1] + betas0[3]*seq(1,3)
medias2 <- betas0[1] + betas0[2] + (betas0[3] + betas0[4])*seq(1,3)

numvueltas <- 200
numind <- 25
am <- 0.1
ros_un <- matrix(c(1,0.5,0.25,0.5,1,0.5,0.25,0.5,1), ncol=3)
ros_ex <- matrix(c(1,am,am,am,1,am,am,am,1), ncol=3)
ros_ind <- matrix(c(1,0,0,0,1,0,0,0,1), ncol=3)
ros_AR1 <- matrix(c(1,am,am^2,am,1,am,am^2,am,1), ncol=3)
rosl <- list(ros_ind, ros_ex, ros_AR1 , ros_un)
indros <- 3

logmedias1 <- exp(medias1)
logmedias2 <- exp(medias2)
logmedias <- list(logmedias1, logmedias2)

beta0s <- c()
beta1s <- c()
beta2s <- c()
beta3s <- c()

for (vueltas in seq(1,numvueltas)) {

  data1 <- data.frame()
  for (casos in 1:2) {
  l <- 1
  Ss <- list()
  betas <- c()
  alphas <- matrix(runif(9),ncol=3)
  conds <- upper.tri(alphas, diag = TRUE)
  for (i in seq(1,3)) {
    for (j in seq(1,3)) {
      if (conds[i,j]) {
        alphas[i,j] <- rosl[[indros]][i,j]*sqrt(logmedias[[casos]][i]*logmedias[[casos]][j])
      } else {
        alphas[i,j] <- NA
      }
    }
  }
  alphas[is.na(alphas)] <- 0
  
  while(sum(is.na(alphas)) != 9) {
    alphas[alphas == 0] <- NA
    
    if (sum(is.na(alphas)) == 9) {
      break
    }
    
    betal <- min(alphas, na.rm = TRUE)
    rs <- which.min(alphas)
    betas[l] <- betal
    aux <- t(matrix(seq(1,3^2),nrow=3))
    rs1 <- which(aux == rs, arr.ind = TRUE)
    
    g1 <- graph_from_adjacency_matrix(alphas, weighted=TRUE)
    S <- c(rs1[[1]], rs1[[2]])
    for (iii in 1:3) {
      if (sum(iii != S) == 2) {
        auxs <- subcomponent(g1 ,iii, mode='all')
        if (sum(auxs == S) == 2) {
          S <- append(S, iii)
          break
        }
      }
    }
    # for (k in seq(1,3)) {
    #   for (i in seq(1,3)) {
    #     for (j in seq(1,3)) {
    #       cond1 <- sum(rs1[[1]] == c(k,i))
    #       cond2 <- sum(rs1[[2]] == c(j,k))
    #       cond3 <- sum(rs1[[1]] == c(i,k))
    #       cond4 <- sum(rs1[[2]] == c(k,j))
    #       if(!(is.na(alphas[i,k])) && !(is.na(alphas[k,j])) && !(is.na(alphas[k,i])) && !(is.na(alphas[j,k])) && cond1 == 1 && cond2 == 1 && cond3 == 1 && cond4 == 1)
    #         S <- append(S,k)
    #     }
    #   }
    # }
    
    # for (k in seq(1,3)) {
    #   cond1 <- !(is.na(alphas[rs1[[1]],k])) && !(is.na(alphas[k,rs1[[2]]]))
    #   cond2 <- !(is.na(alphas[rs1[[2]],k])) && !(is.na(alphas[k,rs1[[1]]]))
    #   cond3 <- !(is.na(alphas[k,rs1[[1]]])) && !(is.na(alphas[rs1[[2]],k]))
    #   cond4 <- !(is.na(alphas[k,rs1[[2]]])) && !(is.na(alphas[rs1[[1]],k]))
    # 
    #   if((cond1 && cond2) || (cond3 && cond4))
    #     S <- append(S,k)
    # }
    
    S <- S[!duplicated(S)]
    Ss <- append(Ss, list(S))
    for (i in S) {
      for (j in S) {
        if (!(is.na(alphas[i,j]))) {
          alphas[i,j] <- alphas[i,j] - betal
        }
      }
    }
    #print (sum(na.omit(alphas) < 0 ))
    l = l+1
  }
  
  T <- matrix(runif(3*length(Ss)),nrow = 3)
  for (j in seq(1,length(Ss))) {
    for (i in Ss[[j]]) {
      T[i,j] <- 1
    }
  }
  T[T != 1] <- 0
  
  Y <- matrix(runif(6*numind),nrow=6)
  for (l in seq(1,length(betas))) {
    Y[l,] <- rpois(numind,betas[l])
  }
  
  Z <- T%*%Y
  
  data <- data.frame(Y=as.numeric(Z),X=rep((casos-1),numind*3), J=rep(1:3, numind))
  data1 <- rbind(data1, data)
  }

  #fit <- glm(Y ~ X + J + X:J, data = data1, family = poisson)
  fit <- geeglm(Y ~ X + J + X:J, data = data1, family = poisson,
                corstr = working[indros], id=seq(1,numind*6))
  beta0s <- c(beta0s, fit$coefficients[1])
  beta1s <- c(beta1s, fit$coefficients[2])
  beta2s <- c(beta2s, fit$coefficients[3])
  beta3s <- c(beta3s, fit$coefficients[4])
}
mean(beta0s)
mean(beta1s)
mean(beta2s)
mean(beta3s)
