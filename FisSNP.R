## Run like this in terminal:
## Rscript GlobalFstSNR.R SNPsMac3NAfiltered.txt

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "PWFst_SNP3_out.txt"
}

## program...

#installed.packages()
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

m <- as.matrix(read.table(args[1], header=TRUE, colClasses = "character"))
recombPoints <- read.table(args[2])
recombPoints <- as.numeric(unlist(recombPoints))

getFrq <- function(m) {
  
  pcol <- rep(0, times = nrow(m))
  qcol <- rep(0, times = nrow(m))
  Pcol <- rep(0, times = nrow(m))
  Qcol <- rep(0, times = nrow(m))
  Hcol <- rep(0, times = nrow(m))
  Pecol <- rep(0, times = nrow(m))
  Qecol <- rep(0, times = nrow(m))
  Hecol <- rep(0, times = nrow(m))
  missing <- rep(0, times = nrow(m))
  
  for (locus in 1:nrow(m)) {
    q <- 0
    p <- 0
    Q <- 0
    P <- 0
    H <- 0
    tot <- 0
    miss <- 0
    for (ind in 1:ncol(m)) {
      if (is.na(m[locus,ind])== TRUE) {
        miss <- miss + 1
      } else if ((m[locus,ind] == "00") == TRUE) {
        p <- p + 2
        P <- P + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "11") == TRUE) {
        q <- q + 2
        Q <- Q + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "10") == TRUE | (m[locus,ind] == "01") == TRUE) {
        q <- q + 1
        p <- p + 1
        H <- H + 1
        tot <- tot + 2
      } else {
        miss <- miss + 1
      }
    }
    qcol[locus] <- q/tot
    pcol[locus] <- p/tot
    Qcol[locus] <- Q/(0.5*tot)
    Pcol[locus] <- P/(0.5*tot)
    Hcol[locus] <- H/(0.5*tot)
    Pecol[locus] <- (p/tot)^2
    Qecol[locus] <- (q/tot)^2
    Hecol[locus] <- 2*(p/tot)*(q/tot)
    missing[locus] <- miss
  }
  
  m_dat <- cbind(p = pcol, q = qcol,  P = Pcol, Q = Qcol, H = Hcol, Pe = Pecol, Qe = Qecol, He = Hecol, NAs = missing)
  
  return(m_dat)
  
}

m1 <- m[,1:15]
m2 <- m[,16:29]
m3 <- m[,30:45]
m4 <- m[,46:61]
m5 <- m[,62:75]
m6 <- m[,76:89]
m7 <- m[,90:103]
m8 <- m[,104:119]
m9 <- m[,120:135]
m10 <- m[,136:150]

m1_dat <- getFrq(m1)
m2_dat <- getFrq(m2)
m3_dat <- getFrq(m3)
m4_dat <- getFrq(m4)
m5_dat <- getFrq(m5)
m6_dat <- getFrq(m6)
m7_dat <- getFrq(m7)
m8_dat <- getFrq(m8)
m9_dat <- getFrq(m9)
m10_dat <- getFrq(m10)


registerDoParallel(32)
permutSNP <- function(gametes, recombPoints) {
  R <- 64
  boot <- c()
  
  getFrqM <- function(M){
    pcol <- rep(0, times = nrow(M))
    qcol <- rep(0, times = nrow(M))
    Pcol <- rep(0, times = nrow(M))
    Qcol <- rep(0, times = nrow(M))
    Hcol <- rep(0, times = nrow(M))
    Pecol <- rep(0, times = nrow(M))
    Qecol <- rep(0, times = nrow(M))
    Hecol <- rep(0, times = nrow(M))
    missing <- rep(0, times = nrow(M))
    
    for (locus in 1:nrow(M)) {
      q <- 0
      p <- 0
      Q <- 0
      P <- 0
      H <- 0
      tot <- 0
      miss <- 0
      for (ind in 1:ncol(M)) {
        if (is.na(M[locus,ind])){
          miss <- miss + 1
        } else if (M[locus,ind] == 0){
          p <- p + 2
          P <- P + 1
        } else if (M[locus,ind] == 1){
          p <- p + 1
          q <- q + 1
          H <- H + 1
        } else if (M[locus,ind] == 2){
          q <- q + 2
          Q <- Q + 1
        }
      }
      tot <- ncol(M) - miss
      qcol[locus] <- q/(tot*2)
      pcol[locus] <- p/(tot*2)
      Qcol[locus] <- Q/tot
      Pcol[locus] <- P/tot
      Hcol[locus] <- H/tot
      Pecol[locus] <- (p/(2*tot))^2
      Qecol[locus] <- (q/(2*tot))^2
      Hecol[locus] <- 2*(p/(2*tot))*(q/(2*tot))
      missing[locus] <- miss
    }
    
    m_dat <- cbind(p = pcol, q = qcol,  P = Pcol, Q = Qcol, H = Hcol, Pe = Pecol, Qe = Qecol, He = Hecol, NAs = missing)
    
    return(m_dat)
    
  }
  
  shuffleGametes <- function(gametes,recombPoints){
    shuffledGametes <- matrix(nrow = nrow(gametes), ncol = ncol(gametes))
    sectionStart <- 1
    sectionLength <- 0
    for (i in 1 :(length(recombPoints)-1)) {
      section <- as.matrix(gametes[(recombPoints[i]:(recombPoints[i+1]-1)),])
      if (nrow(section) == ncol(gametes)) {
        section <- t(section)
      }
      sectionLength <- nrow(section)
      for (j in 1:(ncol(gametes)/2)){
        randomNumber <- sample(c(0,1),1)
        if (randomNumber) {
          tmp <- section[,2*j]
          section[,2*j] <- section[,2*j-1]
          section[,2*j-1] <- tmp
        }
      }
      shuffledGametes[(sectionStart:(sectionStart+sectionLength-1)),] <-  section
      sectionStart <- sectionStart + sectionLength
    }
    return(shuffledGametes)
  }
  
  shuffleColumns <- function(gametes){
    rand <- sample(c(1:ncol(gametes)), replace = F)
    gametes <- gametes[,rand]
  }
  
  gamsToInds <- function(gametes){
    M <- matrix(nrow=nrow(gametes), ncol=(ncol(gametes)/2))
    for (col in 1:ncol(M)){
      M[,col] <- gametes[,(col*2)]+gametes[,(col*2-1)]
    }
    return(M)
  }
  
  boot <- foreach (icount(R), .combine = c, .inorder = F) %dopar% {
    
    recombined <- shuffleGametes(gametes, recombPoints)
    
    pool <- shuffleColumns(recombined)
    
    M <- gamsToInds(pool)

    M_dat <- getFrqM(M)
    
    N <- mean(na.omit(M_dat[,5]))
    D <- mean(na.omit(M_dat[,8]))
    1-(N/D)
    
  }
  return(boot)
}

gametes <- matrix(nrow=nrow(m), ncol=ncol(m)*2)
for (row in 1:nrow(gametes)) {
  for (col in 1:ncol(m)) {
    if (is.na(m[row,col]) == 0) {
      gametes[row, col*2] <- as.numeric(substr(m[row,col],2,2))
      gametes[row, (col*2)-1] <- as.numeric(substr(m[row,col],1,1))
    }
  }
}

gametes1 <- gametes[,(1:30)]
gametes2 <- gametes[,(31:58)]
gametes3 <- gametes[,(59:90)]
gametes4 <- gametes[,(91:122)]
gametes5 <- gametes[,(123:150)]
gametes6 <- gametes[,(151:178)]
gametes7 <- gametes[,(179:206)]
gametes8 <- gametes[,(207:238)]
gametes9 <- gametes[,(239:270)]
gametes10 <- gametes[,(271:300)]

boot1 <- permutSNP(gametes1, recombPoints)
boot2 <- permutSNP(gametes2, recombPoints)
boot3 <- permutSNP(gametes3, recombPoints)
boot4 <- permutSNP(gametes4, recombPoints)
boot5 <- permutSNP(gametes5, recombPoints)
boot6 <- permutSNP(gametes6, recombPoints)
boot7 <- permutSNP(gametes7, recombPoints)
boot8 <- permutSNP(gametes8, recombPoints)
boot9 <- permutSNP(gametes9, recombPoints)
boot10 <- permutSNP(gametes10, recombPoints)

FisTable <- matrix(nrow = 10, ncol = 4)
FisTable[1,1] <- 1 -(mean(m1_dat[,5])/(mean(m1_dat[,8])))
FisTable[2,1] <- 1 -(mean(m2_dat[,5])/(mean(m2_dat[,8])))
FisTable[3,1] <- 1 -(mean(m3_dat[,5])/(mean(m3_dat[,8])))
FisTable[4,1] <- 1 -(mean(m4_dat[,5])/(mean(m4_dat[,8])))
FisTable[5,1] <- 1 -(mean(m5_dat[,5])/(mean(m5_dat[,8])))
FisTable[6,1] <- 1 -(mean(m6_dat[,5])/(mean(m6_dat[,8])))
FisTable[7,1] <- 1 -(mean(m7_dat[,5])/(mean(m7_dat[,8])))
FisTable[8,1] <- 1 -(mean(m8_dat[,5])/(mean(m8_dat[,8])))
FisTable[9,1] <- 1 -(mean(m9_dat[,5])/(mean(m9_dat[,8])))
FisTable[10,1] <- 1 -(mean(m10_dat[,5])/(mean(m10_dat[,8])))

FisTable[1,2] <- mean(boot1 > FisTable[1,1])
FisTable[2,2] <- mean(boot1 > FisTable[2,1])
FisTable[3,2] <- mean(boot1 > FisTable[3,1])
FisTable[4,2] <- mean(boot1 > FisTable[4,1])
FisTable[5,2] <- mean(boot1 > FisTable[5,1])
FisTable[6,2] <- mean(boot1 > FisTable[6,1])
FisTable[7,2] <- mean(boot1 > FisTable[7,1])
FisTable[8,2] <- mean(boot1 > FisTable[8,1])
FisTable[9,2] <- mean(boot1 > FisTable[9,1])
FisTable[10,2] <- mean(boot1 > FisTable[10,1])

FisTable[1,3:4] <- quantile(boot1, c(0.025, 0.975))
FisTable[2,3:4] <- quantile(boot2, c(0.025, 0.975))
FisTable[3,3:4] <- quantile(boot3, c(0.025, 0.975))
FisTable[4,3:4] <- quantile(boot4, c(0.025, 0.975))
FisTable[5,3:4] <- quantile(boot5, c(0.025, 0.975))
FisTable[6,3:4] <- quantile(boot6, c(0.025, 0.975))
FisTable[7,3:4] <- quantile(boot7, c(0.025, 0.975))
FisTable[8,3:4] <- quantile(boot8, c(0.025, 0.975))
FisTable[9,3:4] <- quantile(boot9, c(0.025, 0.975))
FisTable[10,3:4] <- quantile(boot10, c(0.025, 0.975))

print(FisTable) 

