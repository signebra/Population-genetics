## Run like this in terminal:
## Rscript GlobalFstSNR.R SNPsMac3NAfiltered.txt

## program...
library(doParallel)
library(foreach)
library(iterators)

m <- as.matrix(read.table(file = 'SNPsMac15.txt', header=TRUE, colClasses = "character"))

# convert strings "00", "01", "10", and "11" to numbers 0, 1, 1, and 2, respectively

class(m) <- "numeric"

for (row in 1:nrow(m)){
  for (col in 1:ncol(m)){
    if (is.na(m[row,col])){
    } else if (m[row,col] == 10){
      m[row,col] <- 1
    } else if (m[row,col] == 11){
      m[row,col] <- 2
    }
  }
}

# getFrq makes a table of allele and genotype frequencies + expected heterozygosity and number of missing values

getFrq <- function(m) {
  
  pcol <- rep(0, times = nrow(m))
  qcol <- rep(0, times = nrow(m))
  Pcol <- rep(0, times = nrow(m))
  Qcol <- rep(0, times = nrow(m))
  Hcol <- rep(0, times = nrow(m))
  Hecol <- rep(0, times = nrow(m))
  missing <- rep(0, times = nrow(m))
  
  for (locus in 1:nrow(m)) {
    miss <- ncol(m) - length(na.omit(m[locus,]))
    tot <- ncol(m)-miss
    Q <- sum(na.omit(m[locus,]) == 2)
    P <- sum(na.omit(m[locus,]) == 0)
    H <- tot - P - Q
    p <- 2*P+H
    q <- 2*Q+H
    qcol[locus] <- q/(2*tot)
    pcol[locus] <- p/(2*tot)
    Qcol[locus] <- Q/(tot)
    Pcol[locus] <- P/(tot)
    Hcol[locus] <- H/(tot)
    Hecol[locus] <- 2*qcol[locus]*pcol[locus]
    missing[locus] <- miss
  }
  
  m_dat <- cbind(p = pcol, q = qcol,  P = Pcol, Q = Qcol, H = Hcol, He = Hecol, NAs = missing)
  
  return(m_dat)
  
}

m_dat <- getFrq(m)

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

registerDoParallel(20)
bootGlobalFstSNP <- function(m) {
  R <- 1000
  boot <- c()
  
  getFrq <- function(m) {
    
    pcol <- rep(0, times = nrow(m))
    qcol <- rep(0, times = nrow(m))
    Pcol <- rep(0, times = nrow(m))
    Qcol <- rep(0, times = nrow(m))
    Hcol <- rep(0, times = nrow(m))
    Hecol <- rep(0, times = nrow(m))
    missing <- rep(0, times = nrow(m))
    
    for (locus in 1:nrow(m)) {
      miss <- ncol(m) - length(na.omit(m[locus,]))
      tot <- ncol(m)-miss
      Q <- sum(na.omit(m[locus,]) == 2)
      P <- sum(na.omit(m[locus,]) == 0)
      H <- tot - P - Q
      p <- 2*P+H
      q <- 2*Q+H
      qcol[locus] <- q/(2*tot)
      pcol[locus] <- p/(2*tot)
      Qcol[locus] <- Q/(tot)
      Pcol[locus] <- P/(tot)
      Hcol[locus] <- H/(tot)
      Hecol[locus] <- 2*qcol[locus]*pcol[locus]
      missing[locus] <- miss
    }
    
    m_dat <- cbind(p = pcol, q = qcol,  P = Pcol, Q = Qcol, H = Hcol, He = Hecol, NAs = missing)
    
    return(m_dat)
    
  }
  
  m_dat <- getFrq(m)
  
  boot <- foreach (icount(R), .combine = c, .inorder = F) %dopar% {
    
    ind <- c(1:150)
    ind1 <- sample(ind, 15, replace = F)
    M1 <- m[,ind1]
    ind <- ind[! ind %in% ind1]
    ind2 <- sample(ind, 14, replace = F)
    M2 <- m[,ind2]
    ind <- ind[! ind %in% ind2]
    ind3 <- sample(ind, 16, replace = F)
    M3 <- m[,ind3]
    ind <- ind[! ind %in% ind3]
    ind4 <- sample(ind, 16, replace = F)
    M4 <- m[,ind4]
    ind <- ind[! ind %in% ind4]
    ind5 <- sample(ind, 14, replace = F)
    M5 <- m[,ind5]
    ind <- ind[! ind %in% ind5]
    ind6 <- sample(ind, 14, replace = F)
    M6 <- m[,ind6]
    ind <- ind[! ind %in% ind6]
    ind7 <- sample(ind, 14, replace = F)
    M7 <- m[,ind7]
    ind <- ind[! ind %in% ind7]
    ind8 <- sample(ind, 16, replace = F)
    M8 <- m[,ind8]
    ind <- ind[! ind %in% ind8]
    ind9 <- sample(ind, 16, replace = F)
    M9 <- m[,ind9]
    ind <- ind[! ind %in% ind9]
    ind10 <- sample(ind, 15, replace = F)
    M10 <- m[,ind10]
    
    M1_dat <- getFrq(M1)
    M2_dat <- getFrq(M2)
    M3_dat <- getFrq(M3)
    M4_dat <- getFrq(M4)
    M5_dat <- getFrq(M5)
    M6_dat <- getFrq(M6)
    M7_dat <- getFrq(M7)
    M8_dat <- getFrq(M8)
    M9_dat <- getFrq(M9)
    M10_dat <- getFrq(M10)
    
    k <- 10
    Ht <- m_dat[,6]
    Hs <- (M1_dat[,6] + M2_dat[,6] + M3_dat[,6] + M4_dat[,6] +
             M5_dat[,6] + M6_dat[,6] + M7_dat[,6] + M8_dat[,6] +
             M9_dat[,6] + M10_dat[,6])/10
    N <- k * mean(na.omit(Ht - Hs))
    D <- mean(na.omit(((k*Ht) - Hs) * (1- Hs)))
    N/D
    
  }
  return(boot)
}


system.time(boot <- bootGlobalFstSNP(m))

k <- 10
Ht <- m_dat[,6]
Hs <- (m1_dat[,6] + m2_dat[,6] + m3_dat[,6] + m4_dat[,6] +
       m5_dat[,6] + m6_dat[,6] + m7_dat[,6] + m8_dat[,6] +
       m9_dat[,6] + m10_dat[,6])/10
N <- k * mean(Ht - Hs)
D <- mean(((k*Ht) - Hs) * (1- Hs))
Gst <- N/D

p <-  mean(boot > Gst)
hist(boot)

res <-  list(Gst, p, quantile(boot, c(0.025, 0.975)), boot)

print(res) 

