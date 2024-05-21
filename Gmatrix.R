
#### SNPs ####

m <- as.matrix(read.table("SNPsMac3.txt",header=TRUE, colClasses = "character"))

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

m_dat <- getFrq(m)

M <- matrix(nrow = nrow(m), ncol = ncol(m))
for (row in 1:nrow(m)) {
  for (col in 1:ncol(m)) {
    if (is.na(m[row,col])) {
    }
    else if (m[row, col] == as.character("00")) {
      M[row, col] <- 2
    }
    else if (m[row, col] == "01" | m[row, col] == "10") {
      M[row, col] <- 1
    }
    else if (m[row, col] == as.character("11")) {
      M[row, col] <- 0
    }
  } 
}

#### new method ####

m <- M

#### Make S ####
C <- !is.na(m)
S <- matrix(nrow=ncol(C), ncol=ncol(C))
q1 <- c()
q2 <- c()
for (loc in 1:nrow(m)){
  p <- sum(na.omit(m[loc,]))/(2*sum(!is.na(m[loc,])))
  q1 <- append(q1,p)
  q2 <- append(q2,(1-p))
}
for (i in 1:nrow(S)){
  print(i)
  for (j in 1:ncol(S)){
    S[i,j] <- 2*q1%*%((C[,i]*C[,j])*q2)
  }
}

#### Make GRM ####
X <- m
for (row in 1:nrow(m)){
  p <- sum(na.omit(m[row,]))/(2*sum(!is.na(m[loc,])))
  X[row,] <- m[row,] - 2*p
}
X <- t(X)

X[is.na(X)] <- 0

G2 <- (X %*% t(X)) / S

rho <- 0
for (row in 1:nrow(m)){
  p <- sum(na.omit(m[row,]))/(2*sum(!is.na(m[loc,])))
  rho <- rho + p*(1-p)
}
rho

G1 <- (X %*% t(X)) / (2*rho)


write.table(G1,file="GSNP3.txt",row.names=FALSE)

fgrm <- 1-diag(G1)

write.table(fgrm,file="FgrmSNP3.txt",row.names=FALSE) 


#### Haplotypes ####

m <- as.matrix(read.table("Haplo15.txt",header=TRUE, colClasses = "character"))

getFrqHaplo <- function(m) {
  
  acol <- matrix(nrow = nrow(m), ncol =1)
  bcol <- matrix(nrow = nrow(m), ncol =1)
  ccol <- matrix(nrow = nrow(m), ncol =1)
  dcol <- matrix(nrow = nrow(m), ncol =1)
  ecol <- matrix(nrow = nrow(m), ncol =1)
  fcol <- matrix(nrow = nrow(m), ncol =1)
  gcol <- matrix(nrow = nrow(m), ncol =1)
  hcol <- matrix(nrow = nrow(m), ncol =1)
  Acol <- matrix(nrow = nrow(m), ncol =1)
  Bcol <- matrix(nrow = nrow(m), ncol =1)
  Ccol <- matrix(nrow = nrow(m), ncol =1)
  Dcol <- matrix(nrow = nrow(m), ncol =1)
  Ecol <- matrix(nrow = nrow(m), ncol =1)
  Fcol <- matrix(nrow = nrow(m), ncol =1)
  Gcol <- matrix(nrow = nrow(m), ncol =1)
  Hcol <- matrix(nrow = nrow(m), ncol =1)
  Hobscol <- matrix(nrow = nrow(m), ncol =1)
  Hexpcol <- matrix(nrow = nrow(m), ncol =1)
  missing <- matrix(nrow = nrow(m), ncol =1)
  
  for (locus in 1:nrow(m)) {
    a <- 0
    b <- 0
    c <- 0
    d <- 0
    e <- 0
    f <- 0
    g <- 0
    h <- 0
    A <- 0
    B <- 0
    C <- 0
    D <- 0
    E <- 0
    FF <- 0
    G <- 0
    H <- 0
    Hobs <- 0
    tot <- 0
    miss <- 0
    for (ind in 1:ncol(m)) {
      if (is.na(m[locus,ind])== TRUE) {
        miss <- miss + 1
      } else if ((m[locus,ind] == "00")) {
        a <- a + 2
        A <- A + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "11")) {
        b <- b + 2
        B <- B + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "01")) {
        a <- a + 1
        b <- b + 1
        Hobs <- Hobs + 1
        tot <- tot + 2   
      } else if ((m[locus,ind] == "02")) {
        a <- a + 1
        c <- c + 1
        Hobs <- Hobs + 1
        tot <- tot + 2  
      } else if ((m[locus,ind] == "12")) {
        b <- b + 1
        c <- c + 1
        Hobs <- Hobs + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "22")) {
        c <- c + 2
        C <- C + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "03")) {
        a <- a + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2     
      } else if ((m[locus,ind] == "13")) {
        b <- b + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2    
      } else if ((m[locus,ind] == "23")) {
        c <- c + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2   
      } else if ((m[locus,ind] == "33")) {
        d <- d + 2
        D <- D + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "04")) {
        a <- a + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2   
      } else if ((m[locus,ind] == "14")) {
        b <- b + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "24")) {
        c <- c + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "34")) {
        d <- d + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "44")) {
        e <- e + 2
        E <- E + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "05")) {
        a <- a + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2      
      } else if ((m[locus,ind] == "15")) {
        b <- b + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2      
      } else if ((m[locus,ind] == "25")) {
        c <- c + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2 
      } else if ((m[locus,ind] == "35")) {
        d <- d + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2     
      } else if ((m[locus,ind] == "45")) {
        e <- e + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2   
      } else if ((m[locus,ind] == "55")) {
        f <- f + 2
        FF <- FF + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "06")) {
        a <- a + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2       
      } else if ((m[locus,ind] == "16")) {
        b <- b + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2       
      } else if ((m[locus,ind] == "26")) {
        c <- c + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2     
      } else if ((m[locus,ind] == "36")) {
        d <- d + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2       
      } else if ((m[locus,ind] == "46")) {
        e <- e + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2      
      } else if ((m[locus,ind] == "56")) {
        f <- f + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "66")) {
        g <- g + 2
        G <- G + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "07")) {
        a <- a + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "17")) {
        b <- b + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "27")) {
        c <- c + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "37")) {
        d <- d + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "47")) {
        e <- e + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "57")) {
        f <- f + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "67")) {
        g <- g + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "77")) {
        h <- h + 2
        H <- H + 1
        tot <- tot + 2
      } else {
        miss <- miss + 1
      }
    }
    acol[locus] <- a/tot
    bcol[locus] <- b/tot
    ccol[locus] <- c/tot
    dcol[locus] <- d/tot
    ecol[locus] <- e/tot
    fcol[locus] <- f/tot
    gcol[locus] <- g/tot
    hcol[locus] <- h/tot
    Acol[locus] <- A/(0.5*tot)
    Bcol[locus] <- B/(0.5*tot)
    Ccol[locus] <- C/(0.5*tot)
    Dcol[locus] <- D/(0.5*tot)
    Ecol[locus] <- E/(0.5*tot)
    Fcol[locus] <- FF/(0.5*tot)
    Gcol[locus] <- G/(0.5*tot)
    Hcol[locus] <- H/(0.5*tot)
    Hobscol[locus] <- Hobs/(0.5*tot)
    Hexpcol[locus] <- 1 - (a/tot)^2 - (b/tot)^2 - (c/tot)^2 - (d/tot)^2 - (e/tot)^2 - (f/tot)^2 - (g/tot)^2 - (h/tot)^2
    missing[locus] <- miss
  }
  
  m_dat <- cbind(acol,bcol,ccol,dcol,ecol,fcol,gcol,hcol, Acol,Bcol,Ccol,Dcol,Ecol,Fcol,Gcol,Hcol, Hobscol, Hexpcol, missing)
  
  return(m_dat)
  
}

m_dat <- getFrqHaplo(m)

maxvar <- 8

H <- array(dim = c(nrow(m), maxvar, ncol(m)))
for (loc in 1:nrow(m)) {
  for (var in 1:maxvar) {
    for (ind in 1:ncol(m)) {
      H[loc,var,ind] <- 0
    }
  }
} # fills H with 0s

for (loc in 1:nrow(m)) {
  for (var in 1:maxvar) {
    for (ind in 1:ncol(m)) {
      if (is.na(m[loc,ind])== TRUE) {
      } else if ((m[loc,ind] == "00")) {
        H[loc, 1, ind] <- 2
      } else if ((m[loc,ind] == "11")) {
        H[loc, 2, ind] <- 2
      } else if (m[loc,ind] == "01" | (m[loc,ind] == "10")) {
        H[loc, 1, ind] <- 1
        H[loc, 2, ind] <- 1
      } else if ((m[loc,ind] == "02" | (m[loc,ind] == "20"))) {
        H[loc, 1, ind] <- 1
        H[loc, 3, ind] <- 1  
      } else if ((m[loc,ind] == "12") | (m[loc,ind] == "21")) {
        H[loc, 2, ind] <- 1
        H[loc, 3, ind] <- 1        
      } else if ((m[loc,ind] == "22")) {
        H[loc, 3, ind] <- 2
      } else if ((m[loc,ind] == "03") | (m[loc,ind] == "30")) {
        H[loc, 1, ind] <- 1
        H[loc, 4, ind] <- 1    
      } else if ((m[loc,ind] == "13") | (m[loc,ind] == "31")) {
        H[loc, 2, ind] <- 1
        H[loc, 4, ind] <- 1   
      } else if ((m[loc,ind] == "23") | (m[loc,ind] == "32")) {
        H[loc, 3, ind] <- 1
        H[loc, 4, ind] <- 1   
      } else if ((m[loc,ind] == "33")) {
        H[loc, 4, ind] <- 2
      } else if ((m[loc,ind] == "04") | (m[loc,ind] == "40")) {
        H[loc, 1, ind] <- 1
        H[loc, 4, ind] <- 1  
      } else if ((m[loc,ind] == "14") | (m[loc,ind] == "41")) {
        H[loc, 2, ind] <- 1
        H[loc, 5, ind] <- 1       
      } else if ((m[loc,ind] == "24") | (m[loc,ind] == "42")) {
        H[loc, 3, ind] <- 1
        H[loc, 5, ind] <- 1        
      } else if ((m[loc,ind] == "34") | (m[loc,ind] == "43")) {
        H[loc, 4, ind] <- 1
        H[loc, 5, ind] <- 1
      } else if ((m[loc,ind] == "44")) {
        H[loc, 5, ind] <- 2
      } else if ((m[loc,ind] == "05") | (m[loc,ind] == "50")) {
        H[loc, 1, ind] <- 1
        H[loc, 6, ind] <- 1     
      } else if ((m[loc,ind] == "15") | (m[loc,ind] == "51")) {
        H[loc, 2, ind] <- 1
        H[loc, 6, ind] <- 1       
      } else if ((m[loc,ind] == "25") | (m[loc,ind] == "52")) {
        H[loc, 3, ind] <- 1
        H[loc, 6, ind] <- 1 
      } else if ((m[loc,ind] == "35") | (m[loc,ind] == "53")) {
        H[loc, 4, ind] <- 1
        H[loc, 6, ind] <- 1     
      } else if ((m[loc,ind] == "45") | (m[loc,ind] == "54")) {
        H[loc, 5, ind] <- 1
        H[loc, 6, ind] <- 1   
      } else if ((m[loc,ind] == "55")) {
        H[loc, 6, ind] <- 2
      } else if ((m[loc,ind] == "06") | (m[loc,ind] == "60")) {
        H[loc, 1, ind] <- 1
        H[loc, 7, ind] <- 1     
      } else if ((m[loc,ind] == "16") | (m[loc,ind] == "61")) {
        H[loc, 2, ind] <- 1
        H[loc, 7, ind] <- 1       
      } else if ((m[loc,ind] == "26") | (m[loc,ind] == "62")) {
        H[loc, 3, ind] <- 1
        H[loc, 7, ind] <- 1 
      } else if ((m[loc,ind] == "36") | (m[loc,ind] == "63")) {
        H[loc, 4, ind] <- 1
        H[loc, 7, ind] <- 1     
      } else if ((m[loc,ind] == "46") | (m[loc,ind] == "64")) {
        H[loc, 5, ind] <- 1
        H[loc, 7, ind] <- 1   
      } else if ((m[loc,ind] == "56") | (m[loc,ind] == "65")) {
        H[loc, 6, ind] <- 1
        H[loc, 7, ind] <- 1
      } else if ((m[loc,ind] == "66")) {
        H[loc, 7, ind] <- 2
      } else if ((m[loc,ind] == "07") | (m[loc,ind] == "70")) {
        H[loc, 1, ind] <- 1
        H[loc, 8, ind] <- 1     
      } else if ((m[loc,ind] == "17") | (m[loc,ind] == "71")) {
        H[loc, 2, ind] <- 1
        H[loc, 8, ind] <- 1       
      } else if ((m[loc,ind] == "27") | (m[loc,ind] == "72")) {
        H[loc, 3, ind] <- 1
        H[loc, 8, ind] <- 1 
      } else if ((m[loc,ind] == "37") | (m[loc,ind] == "73")) {
        H[loc, 4, ind] <- 1
        H[loc, 8, ind] <- 1     
      } else if ((m[loc,ind] == "47") | (m[loc,ind] == "74")) {
        H[loc, 5, ind] <- 1
        H[loc, 8, ind] <- 1   
      } else if ((m[loc,ind] == "57") | (m[loc,ind] == "75")) {
        H[loc, 6, ind] <- 1
        H[loc, 8, ind] <- 1
      } else if ((m[loc,ind] == "67") | (m[loc,ind] == "76")) {
        H[loc, 7, ind] <- 1
        H[loc, 8, ind] <- 1
      } else if ((m[loc,ind] == "77")) {
        H[loc, 8, ind] <- 2
      } else {
        print("else")
      }
    }
  }
  print(loc)
} # fills in H

#### new method ####

M <- matrix(nrow = (nrow(m)*maxvar), ncol= ncol(m))
for (a in 1:nrow(M)) {
  for (ind in 1:ncol(M)) {
    M[a,ind] <- 0
  }
}
# fills H with 0s

for (loc in 1:nrow(m)) {
  for (ind in 1:ncol(M)) {
    if (is.na(m[loc,ind])== TRUE) {
      M[loc*8-7,ind] <- NA
      M[loc*8-6,ind] <- NA
      M[loc*8-5,ind] <- NA
      M[loc*8-4,ind] <- NA
      M[loc*8-3,ind] <- NA
      M[loc*8-2,ind] <- NA
      M[loc*8-1,ind] <- NA
      M[loc*8-0,ind] <- NA
    } else if ((m[loc,ind] == "00")) {
      M[loc*8-7,ind] <- 2
    } else if ((m[loc,ind] == "11")) {
      M[loc*8-6,ind] <- 2
    } else if (m[loc,ind] == "01" | (m[loc,ind] == "10")) {
      M[loc*8-7,ind] <- 1
      M[loc*8-6,ind] <- 1
    } else if ((m[loc,ind] == "02" | (m[loc,ind] == "20"))) {
      M[loc*8-7,ind] <- 1
      M[loc*8-5,ind] <- 1  
    } else if ((m[loc,ind] == "12") | (m[loc,ind] == "21")) {
      M[loc*8-6,ind] <- 1
      M[loc*8-5,ind] <- 1        
    } else if ((m[loc,ind] == "22")) {
      M[loc*8-5,ind] <- 2
    } else if ((m[loc,ind] == "03") | (m[loc,ind] == "30")) {
      M[loc*8-7,ind] <- 1
      M[loc*8-4,ind] <- 1    
    } else if ((m[loc,ind] == "13") | (m[loc,ind] == "31")) {
      M[loc*8-6,ind] <- 1
      M[loc*8-4,ind] <- 1   
    } else if ((m[loc,ind] == "23") | (m[loc,ind] == "32")) {
      M[loc*8-5,ind] <- 1
      M[loc*8-4,ind] <- 1   
    } else if ((m[loc,ind] == "33")) {
      M[loc*8-4,ind] <- 2
    } else if ((m[loc,ind] == "04") | (m[loc,ind] == "40")) {
      M[loc*8-7,ind] <- 1
      M[loc*8-3,ind] <- 1  
    } else if ((m[loc,ind] == "14") | (m[loc,ind] == "41")) {
      M[loc*8-6,ind] <- 1
      M[loc*8-3,ind] <- 1       
    } else if ((m[loc,ind] == "24") | (m[loc,ind] == "42")) {
      M[loc*8-5,ind] <- 1
      M[loc*8-3,ind] <- 1        
    } else if ((m[loc,ind] == "34") | (m[loc,ind] == "43")) {
      M[loc*8-4,ind] <- 1
      M[loc*8-3,ind] <- 1
    } else if ((m[loc,ind] == "44")) {
      M[loc*8-3,ind] <- 2
    } else if ((m[loc,ind] == "05") | (m[loc,ind] == "50")) {
      M[loc*8-7,ind] <- 1
      M[loc*8-2,ind] <- 1     
    } else if ((m[loc,ind] == "15") | (m[loc,ind] == "51")) {
      M[loc*8-6,ind] <- 1
      M[loc*8-2,ind] <- 1       
    } else if ((m[loc,ind] == "25") | (m[loc,ind] == "52")) {
      M[loc*8-5,ind] <- 1
      M[loc*8-2,ind] <- 1 
    } else if ((m[loc,ind] == "35") | (m[loc,ind] == "53")) {
      M[loc*8-4,ind] <- 1
      M[loc*8-2,ind] <- 1     
    } else if ((m[loc,ind] == "45") | (m[loc,ind] == "54")) {
      M[loc*8-3,ind] <- 1
      M[loc*8-2,ind] <- 1   
    } else if ((m[loc,ind] == "55")) {
      M[loc*8-2,ind] <- 2
    } else if ((m[loc,ind] == "06") | (m[loc,ind] == "60")) {
      M[loc*8-7,ind] <- 1
      M[loc*8-1,ind] <- 1     
    } else if ((m[loc,ind] == "16") | (m[loc,ind] == "61")) {
      M[loc*8-6,ind] <- 1
      M[loc*8-1,ind] <- 1       
    } else if ((m[loc,ind] == "26") | (m[loc,ind] == "62")) {
      M[loc*8-5,ind] <- 1
      M[loc*8-1,ind] <- 1 
    } else if ((m[loc,ind] == "36") | (m[loc,ind] == "63")) {
      M[loc*8-4,ind] <- 1
      M[loc*8-1,ind] <- 1     
    } else if ((m[loc,ind] == "46") | (m[loc,ind] == "64")) {
      M[loc*8-3,ind] <- 1
      M[loc*8-1,ind] <- 1   
    } else if ((m[loc,ind] == "56") | (m[loc,ind] == "65")) {
      M[loc*8-2,ind] <- 1
      M[loc*8-1,ind] <- 1
    } else if ((m[loc,ind] == "66")) {
      M[loc*8-1,ind] <- 2
    } else if ((m[loc,ind] == "07") | (m[loc,ind] == "70")) {
      M[loc*8-7,ind] <- 1
      M[loc*8-0,ind] <- 1     
    } else if ((m[loc,ind] == "17") | (m[loc,ind] == "71")) {
      M[loc*8-6,ind] <- 1
      M[loc*8-0,ind] <- 1       
    } else if ((m[loc,ind] == "27") | (m[loc,ind] == "72")) {
      M[loc*8-5,ind] <- 1
      M[loc*8-0,ind] <- 1 
    } else if ((m[loc,ind] == "37") | (m[loc,ind] == "73")) {
      M[loc*8-4,ind] <- 1
      M[loc*8-0,ind] <- 1     
    } else if ((m[loc,ind] == "47") | (m[loc,ind] == "74")) {
      M[loc*8-3,ind] <- 1
      M[loc*8-0,ind] <- 1   
    } else if ((m[loc,ind] == "57") | (m[loc,ind] == "75")) {
      M[loc*8-2,ind] <- 1
      M[loc*8-0,ind] <- 1
    } else if ((m[loc,ind] == "67") | (m[loc,ind] == "76")) {
      M[loc*8-1,ind] <- 1
      M[loc*8-0,ind] <- 1
    } else if ((m[loc,ind] == "77")) {
      M[loc*8-0,ind] <- 2
    } else {
      print("else")
    }
  }
}

 # fills in M

m <- M

#### remove empty ####

remove <- c()
for (row in 1: nrow(m)){
  if (sum(na.omit(m[row,])) == 0) {
    remove <- append(remove, row)
  }
}
length(remove)
head(remove)
m <- m[-remove,]

#### MAC filtering ####

P <- m[,c(1:ncol(m))]
keep <- c()
for (row in 1:nrow(m)){
  if (sum(na.omit(P[row,]))>=14){
    keep <- append(keep,row)
  }
}
length(keep)
m <- m[keep,]

#### Remove non-polymorphic ####
remove <- c()
for (row in 1:nrow(m)){
  if (sum(na.omit(m[row,]))> 2*ncol(m)-2*sum(is.na(m[row,]))-14){
    remove <- append(remove,row)
  }
}
remove
m <- m[-remove,]

dim(m)

#### Make S ####
C <- !is.na(m)
S <- matrix(nrow=ncol(C), ncol=ncol(C))
q1 <- c()
q2 <- c()
for (loc in 1:nrow(m)){
  p <- sum(na.omit(m[loc,]))/(2*sum(!is.na(m[loc,])))
  q1 <- append(q1,p)
  q2 <- append(q2,(1-p))
}
for (i in 1:nrow(S)){
  print(i)
  for (j in 1:ncol(S)){
    S[i,j] <- 2*q1%*%((C[,i]*C[,j])*q2)
  }
}

#### Make GRM ####
X <- m
for (row in 1:nrow(m)){
  p <- sum(na.omit(m[row,]))/(2*sum(!is.na(m[loc,])))
  X[row,] <- m[row,] - 2*p
}
X <- t(X)

X[is.na(X)] <- 0

G2 <- (X %*% t(X)) / S

rho <- 0
for (row in 1:nrow(m)){
  p <- sum(na.omit(m[row,]))/(2*sum(!is.na(m[loc,])))
  rho <- rho + p*(1-p)
}
rho

G1 <- (X %*% t(X)) / (2*rho)

write.table(G1,file="Ghaplo3.txt",row.names=FALSE) 

fgrm <- diag(G1)

write.table(fgrm,file="FgrmHaplo3.txt",row.names=FALSE) 
