library(doParallel)
library(foreach)
library(iterators) 

df <- (read.table(file ="haplotypesMac3.txt"))

# Make table of allele counts

table <- c(1,0)
for (row in 2:nrow(df)){
  if (as.character(df[row,2] == as.character(df[row-1,2]))){
  }
  else {
    table <- rbind(table, c(row,0))
  }
}
for (row in 1:nrow(table)-1){
  table[row,2] <- table[row+1,1]-table[row,1]
}
table[nrow(table),2] <- nrow(df)-table[nrow(table),1]+1

m <- as.matrix(df[,-c(1:2)])
loci <- (df[,c(1:2)])
class(m) <- "numeric"

for (col in 1:ncol(m)){
  for (row in 1:nrow(m)){
    if (is.na(m[row,col])){
      m[row,col] <- 0
    }
  }
}


getFrq <- function(m, table){
  m_dat <- matrix(0,nrow(table), 10)
  colnames(m_dat) <- c("a","b","c","d","e","f","HE","HO","alleles","missing")
  for (row in 1:nrow(table)){
    m_dat[row,10] <- ncol(m)-(0.5*sum(m[c(table[row,1]:(table[row,1]+table[row,2]-1)),]))
    for (allele in 0:(table[row,2]-1)){
      m_dat[row,allele+1] <- sum(m[(table[row,1]+allele),])/(2*(ncol(m)-m_dat[row,10]))
    }
    m_dat[row,7] <- 1-m_dat[row,1]^2-m_dat[row,2]^2-m_dat[row,3]^2-m_dat[row,4]^2-
      m_dat[row,5]^2-m_dat[row,6]^2
    m_dat[row,8] <- sum((m[c(table[row,1]:(table[row,1]+table[row,2]-1)),])==1)/
      (sum(m[c(table[row,1]:(table[row,1]+table[row,2]-1)),]))
    m_dat[row,9] <- 6 - sum(m_dat[row,c(1:6)] == 0)
    
  }
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

m1_dat <- getFrq(m1, table)
m2_dat <- getFrq(m2, table)
m3_dat <- getFrq(m3, table)
m4_dat <- getFrq(m4, table)
m5_dat <- getFrq(m5, table)
m6_dat <- getFrq(m6, table)
m7_dat <- getFrq(m7, table)
m8_dat <- getFrq(m8, table)
m9_dat <- getFrq(m9, table)
m10_dat <- getFrq(m10, table)
m_dat <- getFrq(m, table)

registerDoParallel(64)
bootGlobalFstHaplo <- function(m) {
  R <- 10000
  boot <- c()
  
  getFrq <- function(m, table){
    m_dat <- matrix(0,nrow(table), 10)
    colnames(m_dat) <- c("a","b","c","d","e","f","HE","HO","alleles","missing")
    for (row in 1:nrow(table)){
      m_dat[row,10] <- ncol(m)-(0.5*sum(m[c(table[row,1]:(table[row,1]+table[row,2]-1)),]))
      for (allele in 0:(table[row,2]-1)){
        m_dat[row,allele+1] <- sum(m[(table[row,1]+allele),])/(2*(ncol(m)-m_dat[row,10]))
      }
      m_dat[row,7] <- 1-m_dat[row,1]^2-m_dat[row,2]^2-m_dat[row,3]^2-m_dat[row,4]^2-
        m_dat[row,5]^2-m_dat[row,6]^2
      m_dat[row,8] <- sum((m[c(table[row,1]:(table[row,1]+table[row,2]-1)),])==1)/
        (sum(m[c(table[row,1]:(table[row,1]+table[row,2]-1)),]))
      m_dat[row,9] <- 6 - sum(m_dat[row,c(1:6)] == 0)
      
    }
    return(m_dat)
  }
  
  m_dat <- getFrq(m,table)
  
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
    
    M1_dat <- getFrq(M1, table)
    M2_dat <- getFrq(M2, table)
    M3_dat <- getFrq(M3, table)
    M4_dat <- getFrq(M4, table)
    M5_dat <- getFrq(M5, table)
    M6_dat <- getFrq(M6, table)
    M7_dat <- getFrq(M7, table)
    M8_dat <- getFrq(M8, table)
    M9_dat <- getFrq(M9, table)
    M10_dat <- getFrq(M10, table)
    
    k <- 10
    Ht <- m_dat[,7]
    Hs <- (M1_dat[,7] + M2_dat[,7] + M3_dat[,7] + M4_dat[,7] +
             M5_dat[,7] + M6_dat[,7] + M7_dat[,7] + M8_dat[,7] +
             M9_dat[,7] + M10_dat[,7])/10
    N <- k * mean(Ht - Hs)
    D <- mean(((k*Ht) - Hs) * (1- Hs))
    N/D
    
  }
  return(boot)
}

system.time(boot <- bootGlobalFstHaplo(m))

k <- 10
Ht <- m_dat[,7]
Hs <- (m1_dat[,7] + m2_dat[,7] + m3_dat[,7] + m4_dat[,7] +
         m5_dat[,7] + m6_dat[,7] + m7_dat[,7] + m8_dat[,7] +
         m9_dat[,7] + m10_dat[,7])/10
N <- k * mean(Ht - Hs)
D <- mean(((k*Ht) - Hs) * (1- Hs))
Gst <- N/D

p <-  mean(boot > Gst)
hist(boot)

res <-  list(Gst, p, quantile(boot, c(0.025, 0.975)), boot)

print(res)  
