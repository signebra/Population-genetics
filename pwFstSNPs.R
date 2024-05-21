
## program...
library(doParallel)
library(foreach)
library(iterators)

df <- as.matrix(read.table(file = 'SNPsMac15.txt', header=TRUE, colClasses = "character"))

loci <- df[,1]
m <- df[,-1]
class(m) <- "numeric"

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
bootPWFstSNP <- function(m1,m2) {
  R <- 80
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
  
  mc <- cbind(m1,m2)
  
  m_dat <- getFrq(mc)
  
  boot <- foreach (icount(R), .combine = c, .inorder = F) %dopar% {
    
    ind <- c(1:ncol(mc))
    ind1 <- sample(ind, ncol(m1), replace = F)
    M1 <- mc[,ind1]
    ind <- ind[! ind %in% ind1]
    ind2 <- sample(ind, ncol(m2), replace = F)
    M2 <- mc[,ind2]
    
    M1_dat <- getFrq(M1)
    M2_dat <- getFrq(M2)
    
    k <- 2
    Ht <- m_dat[,6]
    Hs <- (M1_dat[,6] + M2_dat[,6])/2
    N <- k * mean(na.omit(Ht - Hs))
    D <- mean(na.omit(((k*Ht) - Hs) * (1- Hs)))
    N/D
    
  }
  return(boot)
  
}

GetPWFstSNP <- function(m1, m2, boot) {
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
  mc <- cbind(m1, m2)
  k <- ncol(mc)
  Gst_vector <- matrix(nrow = nrow(mc), ncol = 1)
  m_dat <- getFrq(mc)
  m1_dat <- getFrq(m1)
  m2_dat <- getFrq(m2)
  
  k <- 2
  Ht <- m_dat[,6]
  Hs <- (m1_dat[,6] + m2_dat[,6])/2
  N <- k * mean(na.omit(Ht - Hs))
  D <- mean(na.omit(((k*Ht) - Hs) * (1- Hs)))
  Gst <- N/D
  
  p <-  mean(boot > Gst)
  res <- c(Gst, p)
  return(res)
}

boot1.2 <- bootPWFstSNP(m1,m2) 
boot1.3 <- bootPWFstSNP(m1,m3) 
boot1.4 <- bootPWFstSNP(m1,m4) 
boot1.5 <- bootPWFstSNP(m1,m5) 
boot1.6 <- bootPWFstSNP(m1,m6) 
boot1.7 <- bootPWFstSNP(m1,m7) 
boot1.8 <- bootPWFstSNP(m1,m8) 
boot1.9 <- bootPWFstSNP(m1,m9) 
boot1.10 <- bootPWFstSNP(m1,m10) 
boot2.3 <- bootPWFstSNP(m2,m3)
boot2.4 <- bootPWFstSNP(m2,m4)
boot2.5 <- bootPWFstSNP(m2,m5)
boot2.6 <- bootPWFstSNP(m2,m6)
boot2.7 <- bootPWFstSNP(m2,m7)
boot2.8 <- bootPWFstSNP(m2,m8)
boot2.9 <- bootPWFstSNP(m2,m9)
boot2.10 <- bootPWFstSNP(m2,m10)
boot3.4 <- bootPWFstSNP(m3,m4)
boot3.5 <- bootPWFstSNP(m3,m5)
boot3.6 <- bootPWFstSNP(m3,m6)
boot3.7 <- bootPWFstSNP(m3,m7)
boot3.8 <- bootPWFstSNP(m3,m8)
boot3.9 <- bootPWFstSNP(m3,m9)
boot3.10 <- bootPWFstSNP(m3,m10)
boot4.5 <- bootPWFstSNP(m4,m5)
boot4.6 <- bootPWFstSNP(m4,m6)
boot4.7 <- bootPWFstSNP(m4,m7)
boot4.8 <- bootPWFstSNP(m4,m8)
boot4.9 <- bootPWFstSNP(m4,m9)
boot4.10 <- bootPWFstSNP(m4,m10)
boot5.6 <- bootPWFstSNP(m5,m6)
boot5.7 <- bootPWFstSNP(m5,m7)
boot5.8 <- bootPWFstSNP(m5,m8)
boot5.9 <- bootPWFstSNP(m5,m9)
boot5.10 <- bootPWFstSNP(m5,m10)
boot6.7 <- bootPWFstSNP(m6,m7)
boot6.8 <- bootPWFstSNP(m6,m8)
boot6.9 <- bootPWFstSNP(m6,m9)
boot6.10 <- bootPWFstSNP(m6,m10)
boot7.8 <- bootPWFstSNP(m7,m8)
boot7.9 <- bootPWFstSNP(m7,m9)
boot7.10 <- bootPWFstSNP(m7,m10)
boot8.9 <- bootPWFstSNP(m8,m9)
boot8.10 <- bootPWFstSNP(m8,m10)
boot9.10 <- bootPWFstSNP(m9,m10)

PWFstMatrix <- matrix(nrow = 13, ncol = 10)
res1.2 <- GetPWFstSNP(m1,m2,boot1.2)
PWFstMatrix[1,2] <- res1.2[1]
PWFstMatrix[2,1] <- res1.2[2]
res1.3 <- GetPWFstSNP(m1,m3,boot1.3)
PWFstMatrix[1,3] <- res1.3[1]
PWFstMatrix[3,1] <- res1.3[2]
res1.4 <- GetPWFstSNP(m1,m4,boot1.4)
PWFstMatrix[1,4] <- res1.4[1]
PWFstMatrix[4,1] <- res1.4[2]
res1.5 <- GetPWFstSNP(m1,m5,boot1.5)
PWFstMatrix[1,5] <- res1.5[1]
PWFstMatrix[5,1] <- res1.5[2]
res1.6 <- GetPWFstSNP(m1,m6,boot1.6)
PWFstMatrix[1,6] <- res1.6[1]
PWFstMatrix[6,1] <- res1.6[2]
res1.7 <- GetPWFstSNP(m1,m7,boot1.7)
PWFstMatrix[1,7] <- res1.7[1]
PWFstMatrix[7,1] <- res1.7[2]
res1.8 <- GetPWFstSNP(m1,m8,boot1.8)
PWFstMatrix[1,8] <- res1.8[1]
PWFstMatrix[8,1] <- res1.8[2]
res1.9 <- GetPWFstSNP(m1,m9,boot1.9)
PWFstMatrix[1,9] <- res1.9[1]
PWFstMatrix[9,1] <- res1.9[2]
res1.10 <- GetPWFstSNP(m1,m10,boot1.10)
PWFstMatrix[1,10] <- res1.10[1]
PWFstMatrix[10,1] <- res1.10[2]
res2.3 <- GetPWFstSNP(m2,m3,boot2.3)
PWFstMatrix[2,3] <- res2.3[1]
PWFstMatrix[3,2] <- res2.3[2]
res2.4 <- GetPWFstSNP(m2,m4,boot2.4)
PWFstMatrix[2,4] <- res2.4[1]
PWFstMatrix[4,2] <- res2.4[2]
res2.5 <- GetPWFstSNP(m2,m5,boot2.5)
PWFstMatrix[2,5] <- res2.5[1]
PWFstMatrix[5,2] <- res2.5[2]
res2.6 <- GetPWFstSNP(m2,m6,boot2.6)
PWFstMatrix[2,6] <- res2.6[1]
PWFstMatrix[6,2] <- res2.6[2]
res2.7 <- GetPWFstSNP(m2,m7,boot2.7)
PWFstMatrix[2,7] <- res2.7[1]
PWFstMatrix[7,2] <- res2.7[2]
res2.8 <- GetPWFstSNP(m2,m8,boot2.8)
PWFstMatrix[2,8] <- res2.8[1]
PWFstMatrix[8,2] <- res2.8[2]
res2.9 <- GetPWFstSNP(m2,m9,boot2.9)
PWFstMatrix[2,9] <- res2.9[1]
PWFstMatrix[9,2] <- res2.9[2]
res2.10 <- GetPWFstSNP(m2,m10,boot2.10)
PWFstMatrix[2,10] <- res2.10[1]
PWFstMatrix[10,2] <- res2.10[2]
res3.4 <- GetPWFstSNP(m3,m4,boot3.4)
PWFstMatrix[3,4] <- res3.4[1]
PWFstMatrix[4,3] <- res3.4[2]
res3.5 <- GetPWFstSNP(m3,m5,boot3.5)
PWFstMatrix[3,5] <- res3.5[1]
PWFstMatrix[5,3] <- res3.5[2]
res3.6 <- GetPWFstSNP(m3,m6,boot3.6)
PWFstMatrix[3,6] <- res3.6[1]
PWFstMatrix[6,3] <- res3.6[2]
res3.7 <- GetPWFstSNP(m3,m7,boot3.7)
PWFstMatrix[3,7] <- res3.7[1]
PWFstMatrix[7,3] <- res3.7[2]
res3.8 <- GetPWFstSNP(m3,m8,boot3.8)
PWFstMatrix[3,8] <- res3.8[1]
PWFstMatrix[8,3] <- res3.8[2]
res3.9 <- GetPWFstSNP(m3,m9,boot3.9)
PWFstMatrix[3,9] <- res3.9[1]
PWFstMatrix[9,3] <- res3.9[2]
res3.10 <- GetPWFstSNP(m3,m10,boot3.10)
PWFstMatrix[3,10] <- res3.10[1]
PWFstMatrix[10,3] <- res3.10[2]
res4.5 <- GetPWFstSNP(m4,m5,boot4.5)
PWFstMatrix[4,5] <- res4.5[1]
PWFstMatrix[5,4] <- res4.5[2]
res4.6 <- GetPWFstSNP(m4,m6,boot4.6)
PWFstMatrix[4,6] <- res4.6[1]
PWFstMatrix[6,4] <- res4.6[2]
res4.7 <- GetPWFstSNP(m4,m7,boot4.7)
PWFstMatrix[4,7] <- res4.7[1]
PWFstMatrix[7,4] <- res4.7[2]
res4.8 <- GetPWFstSNP(m4,m8,boot4.8)
PWFstMatrix[4,8] <- res4.8[1]
PWFstMatrix[8,4] <- res4.8[2]
res4.9 <- GetPWFstSNP(m4,m9,boot4.9)
PWFstMatrix[4,9] <- res4.9[1]
PWFstMatrix[9,4] <- res4.9[2]
res4.10 <- GetPWFstSNP(m4,m10,boot4.10)
PWFstMatrix[4,10] <- res4.10[1]
PWFstMatrix[10,4] <- res4.10[2]
res5.6 <- GetPWFstSNP(m5,m6,boot5.6)
PWFstMatrix[5,6] <- res5.6[1]
PWFstMatrix[6,5] <- res5.6[2]
res5.7 <- GetPWFstSNP(m5,m7,boot5.7)
PWFstMatrix[5,7] <- res5.7[1]
PWFstMatrix[7,5] <- res5.7[2]
res5.8 <- GetPWFstSNP(m5,m8,boot5.8)
PWFstMatrix[5,8] <- res5.8[1]
PWFstMatrix[8,5] <- res5.8[2]
res5.9 <- GetPWFstSNP(m5,m9,boot5.9)
PWFstMatrix[5,9] <- res5.9[1]
PWFstMatrix[9,5] <- res5.9[2]
res5.10 <- GetPWFstSNP(m5,m10,boot5.10)
PWFstMatrix[5,10] <- res5.10[1]
PWFstMatrix[10,5] <- res5.10[2]
res6.7 <- GetPWFstSNP(m6,m7,boot6.7)
PWFstMatrix[6,7] <- res6.7[1]
PWFstMatrix[7,6] <- res6.7[2]
res6.8 <- GetPWFstSNP(m6,m8,boot6.8)
PWFstMatrix[6,8] <- res6.8[1]
PWFstMatrix[8,6] <- res6.8[2]
res6.9 <- GetPWFstSNP(m6,m9,boot6.9)
PWFstMatrix[6,9] <- res6.9[1]
PWFstMatrix[9,6] <- res6.9[2]
res6.10 <- GetPWFstSNP(m6,m10,boot6.10)
PWFstMatrix[6,10] <- res6.10[1]
PWFstMatrix[10,6] <- res6.10[2]
res7.8 <- GetPWFstSNP(m7,m8,boot7.8)
PWFstMatrix[7,8] <- res7.8[1]
PWFstMatrix[8,7] <- res7.8[2]
res7.9 <- GetPWFstSNP(m7,m9,boot7.9)
PWFstMatrix[7,9] <- res7.9[1]
PWFstMatrix[9,7] <- res7.9[2]
res7.10 <- GetPWFstSNP(m7,m10,boot7.10)
PWFstMatrix[7,10] <- res7.10[1]
PWFstMatrix[10,7] <- res7.10[2]
res8.9 <- GetPWFstSNP(m8,m9,boot8.9)
PWFstMatrix[8,9] <- res8.9[1]
PWFstMatrix[9,8] <- res8.9[2]
res8.10 <- GetPWFstSNP(m8,m10,boot8.10)
PWFstMatrix[8,10] <- res8.10[1]
PWFstMatrix[10,8] <- res8.10[2]
res9.10 <- GetPWFstSNP(m9,m10,boot9.10)
PWFstMatrix[9,10] <- res9.10[1]
PWFstMatrix[10,9] <- res9.10[2]
99
PWFstMatrix[12,1] <- mean(c(res1.2[1],res1.3[1],res1.4[1],res1.5[1],res1.6[1],res1.7[1],res1.8[1],res1.9[1],res1.10[1]))
PWFstMatrix[12,2] <- mean(c(res1.2[1],res2.3[1],res2.4[1],res2.5[1],res2.6[1],res2.7[1],res2.8[1],res2.9[1],res2.10[1]))
PWFstMatrix[12,3] <- mean(c(res1.3[1],res2.3[1],res3.4[1],res3.5[1],res3.6[1],res3.7[1],res3.8[1],res3.9[1],res3.10[1]))
PWFstMatrix[12,4] <- mean(c(res1.4[1],res2.4[1],res3.4[1],res4.5[1],res4.6[1],res4.7[1],res4.8[1],res4.9[1],res4.10[1]))
PWFstMatrix[12,5] <- mean(c(res1.5[1],res2.5[1],res3.5[1],res4.5[1],res5.6[1],res5.7[1],res5.8[1],res5.9[1],res5.10[1]))
PWFstMatrix[12,6] <- mean(c(res1.6[1],res2.6[1],res3.6[1],res4.6[1],res5.6[1],res6.7[1],res6.8[1],res6.9[1],res6.10[1]))
PWFstMatrix[12,7] <- mean(c(res1.7[1],res2.7[1],res3.7[1],res4.7[1],res5.7[1],res6.7[1],res7.8[1],res7.9[1],res7.10[1]))
PWFstMatrix[12,8] <- mean(c(res1.8[1],res2.8[1],res3.8[1],res4.8[1],res5.8[1],res6.8[1],res7.8[1],res8.9[1],res8.10[1]))
PWFstMatrix[12,9] <- mean(c(res1.9[1],res2.9[1],res3.9[1],res4.9[1],res5.9[1],res6.9[1],res7.9[1],res8.9[1],res9.10[1]))
PWFstMatrix[12,10] <- mean(c(res1.10[1],res2.10[1],res3.10[1],res4.10[1],res5.10[1],res6.10[1],res7.10[1],res8.10[1],res9.10[1]))

boot1 <- rowMeans(cbind(boot1.2,boot1.3,boot1.4,boot1.5,boot1.6,boot1.7,boot1.8,boot1.9,boot1.10))
boot2 <- rowMeans(cbind(boot1.2,boot2.3,boot2.4,boot2.5,boot2.6,boot2.7,boot2.8,boot2.9,boot2.10))
boot3 <- rowMeans(cbind(boot1.3,boot2.3,boot3.4,boot3.5,boot3.6,boot3.7,boot3.8,boot3.9,boot3.10))
boot4 <- rowMeans(cbind(boot1.4,boot2.4,boot3.4,boot4.5,boot4.6,boot4.7,boot4.8,boot4.9,boot4.10))
boot5 <- rowMeans(cbind(boot1.5,boot2.5,boot3.5,boot4.5,boot5.6,boot5.7,boot5.8,boot5.9,boot5.10))
boot6 <- rowMeans(cbind(boot1.6,boot2.6,boot3.6,boot4.6,boot5.6,boot6.7,boot6.8,boot6.9,boot6.10))
boot7 <- rowMeans(cbind(boot1.7,boot2.7,boot3.7,boot4.7,boot5.7,boot6.7,boot7.8,boot7.9,boot7.10))
boot8 <- rowMeans(cbind(boot1.8,boot2.8,boot3.8,boot4.8,boot5.8,boot6.8,boot7.8,boot8.9,boot8.10))
boot9 <- rowMeans(cbind(boot1.9,boot2.9,boot3.9,boot4.9,boot5.9,boot6.9,boot7.9,boot8.9,boot9.10))
boot10 <- rowMeans(cbind(boot1.10,boot2.10,boot3.10,boot4.10,boot5.10,boot6.10,boot7.10,boot8.10,boot9.10))

PWFstMatrix[13,1] <- mean(boot1 > PWFstMatrix[12,1])
PWFstMatrix[13,2] <- mean(boot2 > PWFstMatrix[12,2])
PWFstMatrix[13,3] <- mean(boot3 > PWFstMatrix[12,3])
PWFstMatrix[13,4] <- mean(boot4 > PWFstMatrix[12,4])
PWFstMatrix[13,5] <- mean(boot5 > PWFstMatrix[12,5])
PWFstMatrix[13,6] <- mean(boot6 > PWFstMatrix[12,6])
PWFstMatrix[13,7] <- mean(boot7 > PWFstMatrix[12,7])
PWFstMatrix[13,8] <- mean(boot8 > PWFstMatrix[12,8])
PWFstMatrix[13,9] <- mean(boot9 > PWFstMatrix[12,9])
PWFstMatrix[13,10] <- mean(boot10 > PWFstMatrix[12,10])

print(PWFstMatrix)

