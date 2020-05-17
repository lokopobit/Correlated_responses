library(dplyr)
library(geepack)
library(readr)

working <- c('independence','exchangeable','ar1','unstructured','userdefined')
std.errs <- c('san.se',"fij",'jack', 'j1s')

### INICIALMENTE: VECTOR DE MEDIAS Y MATRIZ CORRELACION
k <- 3
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
indros <- 1

logmedias1 <- 1/(1+exp(-medias1))
logmedias2 <- 1/(1+exp(-medias2))
logmedias <- list(logmedias1, logmedias2)

beta0s <- c()
beta1s <- c()
beta2s <- c()
beta3s <- c()

x0 <- c()
x1 <- c()

for (vueltas in seq(1,numvueltas)) {
  
  data1 <- data.frame()
  for (casos in 1:2) {
  ### PASO 0: INICIALIZAR ALPHAS
    alphas <- matrix(runif(k^2), ncol=k)
    for (i in seq(1,k)) {
      for (j in seq(1,k)) {
        aux <- (1-logmedias[[casos]][i])/logmedias[[casos]][i]
        aux1 <- (1-logmedias[[casos]][j])/logmedias[[casos]][j]
        alphas[i,j] <- log(1 + rosl[[indros]][i,j]*sqrt(aux*aux1))
      }
    }
  
  ### PASO 1: DEFINIR Tl, Betal y Sl
    
    
    l = 0
    aux00 <- TRUE
    betals <- c()
    Sls <- list()
    
    while (aux00) {
      l<-l+1
      # if (l == 15) {
      # break
      # }
      Tl <- list()
      Sl <- list()
      for (i in seq(1,k)) {
        for (j in seq(1,k)) {
          if (alphas[i,j] > 0) {
            Tl <- append(Tl, alphas[i,j])
            Sl <- append(Sl, list(c(i,j))) #segun internet
          } else {
            Tl <- append(Tl, NA)
          }
        }
      }
      
      Tl <- as.numeric(Tl)
      betal <- min(Tl, na.rm =TRUE)
      rs <- which.min(Tl)
      betals[l] <- betal
      auxm <- t(matrix(seq(1,k^2),nrow=k))
      rs1 <- which(auxm == rs, arr.ind = TRUE)
      
      cond1 <- abs(alphas[rs1[1],rs1[1]]) < 1e-10
      cond2 <- abs(alphas[rs1[2],rs1[2]]) < 1e-10
      if (cond1 | cond2) {
        break
      }
      
      #CALCULO DE Sl y ACTUALIZACION DE LOS alphas
      
      Sl <- append(Sl, list(c(rs1[1],rs1[2])))
      Sl <- Sl[!duplicated(Sl)]
      Sls[[l]] <- Sl
      for (k1 in seq(1,length(Sl))) {
        alphas[Sl[[k1]][1],Sl[[k1]][2]] <- alphas[Sl[[k1]][1],Sl[[k1]][2]] - betal
      }
      
      aux00 <- sum(alphas == 0) != (k^2)
      aux_ext <- TRUE
      alphas1 <- 0
      betals1  <- 0
      Sl1 <- 0
      
      if (aux00 == FALSE) {
        aux_ext <- FALSE
        alphas1 <- alphas
        betals1 <- betals
        Sl1 <- Sls
      }
    }
  
    Yes <- c()
    for (ind in seq(1:numind)) {
    for (i in seq(1,k)) {
      auxb <- 0
      aux1b <- 0
      for (l in seq(1,length(betals1))) {
        Sl1new <- c()
        for (l1 in seq(1,length(Sl1[[l]]))) {
          Sl1new <- c(Sl1new, Sl1[[l]][[l1]])
        }
        aux1b <- any(Sl1new == i)
        auxb <- auxb + rpois(1, betals1[l])*aux1b
      }
      Yes <- c(Yes, auxb)
    }
    Zs <- ifelse(Yes == 0,1,0)
    }
    
    data <- data.frame(Y=as.numeric(Zs), X=rep((casos-1),numind*3), J=rep(1:3, numind))
    data1 <- rbind(data1, data)

  }
  
  #fit <- glm(Y ~ X + J + X:J, data = data1, family = binomial)
  fit <- geeglm(Y ~ X + J + X:J, data = data1, family = binomial,
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
