## Run like this in terminal:
## Rscript Fst_values.R SNPsMac3NAfiltered.txt

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "PWFst_SNP3_out.txt"
}

## program...
library(doParallel)
library(foreach)
library(iterators)
m <- as.matrix(read.table(args[1], header=TRUE, colClasses = "character"))

getFrqHaplo <- function(m) {
  
  acol <- rep(0, times = nrow(m))
  bcol <- rep(0, times = nrow(m))
  ccol <- rep(0, times = nrow(m))
  dcol <- rep(0, times = nrow(m))
  ecol <- rep(0, times = nrow(m))
  fcol <- rep(0, times = nrow(m))
  gcol <- rep(0, times = nrow(m))
  hcol <- rep(0, times = nrow(m))
  Acol <- rep(0, times = nrow(m))
  Bcol <- rep(0, times = nrow(m))
  Ccol <- rep(0, times = nrow(m))
  Dcol <- rep(0, times = nrow(m))
  Ecol <- rep(0, times = nrow(m))
  Fcol <- rep(0, times = nrow(m))
  Gcol <- rep(0, times = nrow(m))
  Hcol <- rep(0, times = nrow(m))
  Hobscol <- rep(0, times = nrow(m))
  Hexpcol <- rep(0, times = nrow(m))
  missing <- rep(0, times = nrow(m))
  
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
      } else if ((m[locus,ind] == "22")) {
        c <- c + 2
        C <- C + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "33")) {
        d <- d + 2
        D <- D + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "44")) {
        e <- e + 2
        E <- E + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "55")) {
        f <- f + 2
        FF <- FF + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "66")) {
        g <- g + 2
        G <- G + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "77")) {
        h <- h + 2
        H <- H + 1
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
      } else if ((m[locus,ind] == "03")) {
        a <- a + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "04")) {
        a <- a + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "05")) {
        a <- a + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "06")) {
        a <- a + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "07")) {
        a <- a + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "12")) {
        b <- b + 1
        c <- c + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "13")) {
        b <- b + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "14")) {
        b <- b + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "15")) {
        b <- b + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "16")) {
        b <- b + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "17")) {
        b <- b + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "23")) {
        c <- c + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "24")) {
        c <- c + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "25")) {
        c <- c + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "26")) {
        c <- c + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "27")) {
        c <- c + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "34")) {
        d <- d + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "35")) {
        d <- d + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "36")) {
        d <- d + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "37")) {
        d <- d + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "45")) {
        e <- e + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "46")) {
        e <- e + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "47")) {
        e <- e + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "56")) {
        f <- f + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "57")) {
        f <- f + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "66")) {
        g <- g + 1
        h <- h + 1
        Hobs <- Hobs + 1
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


registerDoParallel(32)
bootPWFstHaplo <- function(m1, m2) {
  R <- 10000
  boot <- c()
  getFrqHaplo <- function(m) {
    
    acol <- rep(0, times = nrow(m))
    bcol <- rep(0, times = nrow(m))
    ccol <- rep(0, times = nrow(m))
    dcol <- rep(0, times = nrow(m))
    ecol <- rep(0, times = nrow(m))
    fcol <- rep(0, times = nrow(m))
    gcol <- rep(0, times = nrow(m))
    hcol <- rep(0, times = nrow(m))
    Acol <- rep(0, times = nrow(m))
    Bcol <- rep(0, times = nrow(m))
    Ccol <- rep(0, times = nrow(m))
    Dcol <- rep(0, times = nrow(m))
    Ecol <- rep(0, times = nrow(m))
    Fcol <- rep(0, times = nrow(m))
    Gcol <- rep(0, times = nrow(m))
    Hcol <- rep(0, times = nrow(m))
    Hobscol <- rep(0, times = nrow(m))
    Hexpcol <- rep(0, times = nrow(m))
    missing <- rep(0, times = nrow(m))
    
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
        } else if ((m[locus,ind] == "22")) {
          c <- c + 2
          C <- C + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "33")) {
          d <- d + 2
          D <- D + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "44")) {
          e <- e + 2
          E <- E + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "55")) {
          f <- f + 2
          FF <- FF + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "66")) {
          g <- g + 2
          G <- G + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "77")) {
          h <- h + 2
          H <- H + 1
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
        } else if ((m[locus,ind] == "03")) {
          a <- a + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "04")) {
          a <- a + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "05")) {
          a <- a + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "06")) {
          a <- a + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "07")) {
          a <- a + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "12")) {
          b <- b + 1
          c <- c + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "13")) {
          b <- b + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "14")) {
          b <- b + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "15")) {
          b <- b + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "16")) {
          b <- b + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "17")) {
          b <- b + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "23")) {
          c <- c + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "24")) {
          c <- c + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "25")) {
          c <- c + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "26")) {
          c <- c + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "27")) {
          c <- c + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "34")) {
          d <- d + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "35")) {
          d <- d + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "36")) {
          d <- d + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "37")) {
          d <- d + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "45")) {
          e <- e + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "46")) {
          e <- e + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "47")) {
          e <- e + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "56")) {
          f <- f + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "57")) {
          f <- f + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "66")) {
          g <- g + 1
          h <- h + 1
          Hobs <- Hobs + 1
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
  mc <- cbind(m1,m2)
  m_dat <- getFrqHaplo(mc)
  
  boot <- foreach (icount(R), .combine = c, .inorder = F) %dopar% {
    
    ind <- c(1:ncol(mc))
    ind1 <- sample(ind, ncol(m1), replace = F)
    M1 <- mc[,ind1]
    ind <- ind[! ind %in% ind1]
    ind2 <- sample(ind, ncol(m2), replace = F)
    M2 <- mc[,ind2]
    
    M1_dat <- getFrqHaplo(M1)
    M2_dat <- getFrqHaplo(M2)
        
    k <- 2
    Ht <- m_dat[,18]
    Hs <- (M1_dat[,18] + M2_dat[,18])/2
    N <- k * mean(na.omit((Ht - Hs)))
    D <- mean(na.omit(((k*Ht) - Hs) * (1- Hs)))
    N/D

  }
  
  return(boot)
}

GetPWFstHaplo <- function(m1, m2, boot) {
  getFrqHaplo <- function(m) {
    
    acol <- rep(0, times = nrow(m))
    bcol <- rep(0, times = nrow(m))
    ccol <- rep(0, times = nrow(m))
    dcol <- rep(0, times = nrow(m))
    ecol <- rep(0, times = nrow(m))
    fcol <- rep(0, times = nrow(m))
    gcol <- rep(0, times = nrow(m))
    hcol <- rep(0, times = nrow(m))
    Acol <- rep(0, times = nrow(m))
    Bcol <- rep(0, times = nrow(m))
    Ccol <- rep(0, times = nrow(m))
    Dcol <- rep(0, times = nrow(m))
    Ecol <- rep(0, times = nrow(m))
    Fcol <- rep(0, times = nrow(m))
    Gcol <- rep(0, times = nrow(m))
    Hcol <- rep(0, times = nrow(m))
    Hobscol <- rep(0, times = nrow(m))
    Hexpcol <- rep(0, times = nrow(m))
    missing <- rep(0, times = nrow(m))
    
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
        } else if ((m[locus,ind] == "22")) {
          c <- c + 2
          C <- C + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "33")) {
          d <- d + 2
          D <- D + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "44")) {
          e <- e + 2
          E <- E + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "55")) {
          f <- f + 2
          FF <- FF + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "66")) {
          g <- g + 2
          G <- G + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "77")) {
          h <- h + 2
          H <- H + 1
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
        } else if ((m[locus,ind] == "03")) {
          a <- a + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "04")) {
          a <- a + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "05")) {
          a <- a + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "06")) {
          a <- a + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "07")) {
          a <- a + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "12")) {
          b <- b + 1
          c <- c + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "13")) {
          b <- b + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "14")) {
          b <- b + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "15")) {
          b <- b + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "16")) {
          b <- b + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "17")) {
          b <- b + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "23")) {
          c <- c + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "24")) {
          c <- c + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "25")) {
          c <- c + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "26")) {
          c <- c + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "27")) {
          c <- c + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "34")) {
          d <- d + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "35")) {
          d <- d + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "36")) {
          d <- d + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "37")) {
          d <- d + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "45")) {
          e <- e + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "46")) {
          e <- e + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "47")) {
          e <- e + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "56")) {
          f <- f + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "57")) {
          f <- f + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "66")) {
          g <- g + 1
          h <- h + 1
          Hobs <- Hobs + 1
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
  mc <- cbind(m1, m2)
  k <- ncol(mc)
  Gst_vector <- matrix(nrow = nrow(mc), ncol = 1)
  m_dat <- getFrqHaplo(mc)
  m1_dat <- getFrqHaplo(m1)
  m2_dat <- getFrqHaplo(m2)
  
  k <- 2
  Ht <- m_dat[,18]
  Hs <- (m1_dat[,18] + m2_dat[,18])/2
  N <- k * mean(na.omit((Ht - Hs)))
  D <- mean(na.omit(((k*Ht) - Hs) * (1- Hs)))
  Gst <- N/D
  
  p <-  mean(boot > Gst)
  res <- c(Gst, p)
  return(res)
}

boot1.2 <- bootPWFstHaplo(m1,m2) 
boot1.3 <- bootPWFstHaplo(m1,m3) 
boot1.4 <- bootPWFstHaplo(m1,m4) 
boot1.5 <- bootPWFstHaplo(m1,m5) 
boot1.6 <- bootPWFstHaplo(m1,m6) 
boot1.7 <- bootPWFstHaplo(m1,m7) 
boot1.8 <- bootPWFstHaplo(m1,m8) 
boot1.9 <- bootPWFstHaplo(m1,m9) 
boot1.10 <- bootPWFstHaplo(m1,m10) 
boot2.3 <- bootPWFstHaplo(m2,m3)
boot2.4 <- bootPWFstHaplo(m2,m4)
boot2.5 <- bootPWFstHaplo(m2,m5)
boot2.6 <- bootPWFstHaplo(m2,m6)
boot2.7 <- bootPWFstHaplo(m2,m7)
boot2.8 <- bootPWFstHaplo(m2,m8)
boot2.9 <- bootPWFstHaplo(m2,m9)
boot2.10 <- bootPWFstHaplo(m2,m10)
boot3.4 <- bootPWFstHaplo(m3,m4)
boot3.5 <- bootPWFstHaplo(m3,m5)
boot3.6 <- bootPWFstHaplo(m3,m6)
boot3.7 <- bootPWFstHaplo(m3,m7)
boot3.8 <- bootPWFstHaplo(m3,m8)
boot3.9 <- bootPWFstHaplo(m3,m9)
boot3.10 <- bootPWFstHaplo(m3,m10)
boot4.5 <- bootPWFstHaplo(m4,m5)
boot4.6 <- bootPWFstHaplo(m4,m6)
boot4.7 <- bootPWFstHaplo(m4,m7)
boot4.8 <- bootPWFstHaplo(m4,m8)
boot4.9 <- bootPWFstHaplo(m4,m9)
boot4.10 <- bootPWFstHaplo(m4,m10)
boot5.6 <- bootPWFstHaplo(m5,m6)
boot5.7 <- bootPWFstHaplo(m5,m7)
boot5.8 <- bootPWFstHaplo(m5,m8)
boot5.9 <- bootPWFstHaplo(m5,m9)
boot5.10 <- bootPWFstHaplo(m5,m10)
boot6.7 <- bootPWFstHaplo(m6,m7)
boot6.8 <- bootPWFstHaplo(m6,m8)
boot6.9 <- bootPWFstHaplo(m6,m9)
boot6.10 <- bootPWFstHaplo(m6,m10)
boot7.8 <- bootPWFstHaplo(m7,m8)
boot7.9 <- bootPWFstHaplo(m7,m9)
boot7.10 <- bootPWFstHaplo(m7,m10)
boot8.9 <- bootPWFstHaplo(m8,m9)
boot8.10 <- bootPWFstHaplo(m8,m10)
boot9.10 <- bootPWFstHaplo(m9,m10)

PWFstMatrix <- matrix(nrow = 13, ncol = 10)
res1.2 <- GetPWFstHaplo(m1,m2,boot1.2)
PWFstMatrix[1,2] <- res1.2[1]
PWFstMatrix[2,1] <- res1.2[2]
res1.3 <- GetPWFstHaplo(m1,m3,boot1.3)
PWFstMatrix[1,3] <- res1.3[1]
PWFstMatrix[3,1] <- res1.3[2]
res1.4 <- GetPWFstHaplo(m1,m4,boot1.4)
PWFstMatrix[1,4] <- res1.4[1]
PWFstMatrix[4,1] <- res1.4[2]
res1.5 <- GetPWFstHaplo(m1,m5,boot1.5)
PWFstMatrix[1,5] <- res1.5[1]
PWFstMatrix[5,1] <- res1.5[2]
res1.6 <- GetPWFstHaplo(m1,m6,boot1.6)
PWFstMatrix[1,6] <- res1.6[1]
PWFstMatrix[6,1] <- res1.6[2]
res1.7 <- GetPWFstHaplo(m1,m7,boot1.7)
PWFstMatrix[1,7] <- res1.7[1]
PWFstMatrix[7,1] <- res1.7[2]
res1.8 <- GetPWFstHaplo(m1,m8,boot1.8)
PWFstMatrix[1,8] <- res1.8[1]
PWFstMatrix[8,1] <- res1.8[2]
res1.9 <- GetPWFstHaplo(m1,m9,boot1.9)
PWFstMatrix[1,9] <- res1.9[1]
PWFstMatrix[9,1] <- res1.9[2]
res1.10 <- GetPWFstHaplo(m1,m10,boot1.10)
PWFstMatrix[1,10] <- res1.10[1]
PWFstMatrix[10,1] <- res1.10[2]
res2.3 <- GetPWFstHaplo(m2,m3,boot2.3)
PWFstMatrix[2,3] <- res2.3[1]
PWFstMatrix[3,2] <- res2.3[2]
res2.4 <- GetPWFstHaplo(m2,m4,boot2.4)
PWFstMatrix[2,4] <- res2.4[1]
PWFstMatrix[4,2] <- res2.4[2]
res2.5 <- GetPWFstHaplo(m2,m5,boot2.5)
PWFstMatrix[2,5] <- res2.5[1]
PWFstMatrix[5,2] <- res2.5[2]
res2.6 <- GetPWFstHaplo(m2,m6,boot2.6)
PWFstMatrix[2,6] <- res2.6[1]
PWFstMatrix[6,2] <- res2.6[2]
res2.7 <- GetPWFstHaplo(m2,m7,boot2.7)
PWFstMatrix[2,7] <- res2.7[1]
PWFstMatrix[7,2] <- res2.7[2]
res2.8 <- GetPWFstHaplo(m2,m8,boot2.8)
PWFstMatrix[2,8] <- res2.8[1]
PWFstMatrix[8,2] <- res2.8[2]
res2.9 <- GetPWFstHaplo(m2,m9,boot2.9)
PWFstMatrix[2,9] <- res2.9[1]
PWFstMatrix[9,2] <- res2.9[2]
res2.10 <- GetPWFstHaplo(m2,m10,boot2.10)
PWFstMatrix[2,10] <- res2.10[1]
PWFstMatrix[10,2] <- res2.10[2]
res3.4 <- GetPWFstHaplo(m3,m4,boot3.4)
PWFstMatrix[3,4] <- res3.4[1]
PWFstMatrix[4,3] <- res3.4[2]
res3.5 <- GetPWFstHaplo(m3,m5,boot3.5)
PWFstMatrix[3,5] <- res3.5[1]
PWFstMatrix[5,3] <- res3.5[2]
res3.6 <- GetPWFstHaplo(m3,m6,boot3.6)
PWFstMatrix[3,6] <- res3.6[1]
PWFstMatrix[6,3] <- res3.6[2]
res3.7 <- GetPWFstHaplo(m3,m7,boot3.7)
PWFstMatrix[3,7] <- res3.7[1]
PWFstMatrix[7,3] <- res3.7[2]
res3.8 <- GetPWFstHaplo(m3,m8,boot3.8)
PWFstMatrix[3,8] <- res3.8[1]
PWFstMatrix[8,3] <- res3.8[2]
res3.9 <- GetPWFstHaplo(m3,m9,boot3.9)
PWFstMatrix[3,9] <- res3.9[1]
PWFstMatrix[9,3] <- res3.9[2]
res3.10 <- GetPWFstHaplo(m3,m10,boot3.10)
PWFstMatrix[3,10] <- res3.10[1]
PWFstMatrix[10,3] <- res3.10[2]
res4.5 <- GetPWFstHaplo(m4,m5,boot4.5)
PWFstMatrix[4,5] <- res4.5[1]
PWFstMatrix[5,4] <- res4.5[2]
res4.6 <- GetPWFstHaplo(m4,m6,boot4.6)
PWFstMatrix[4,6] <- res4.6[1]
PWFstMatrix[6,4] <- res4.6[2]
res4.7 <- GetPWFstHaplo(m4,m7,boot4.7)
PWFstMatrix[4,7] <- res4.7[1]
PWFstMatrix[7,4] <- res4.7[2]
res4.8 <- GetPWFstHaplo(m4,m8,boot4.8)
PWFstMatrix[4,8] <- res4.8[1]
PWFstMatrix[8,4] <- res4.8[2]
res4.9 <- GetPWFstHaplo(m4,m9,boot4.9)
PWFstMatrix[4,9] <- res4.9[1]
PWFstMatrix[9,4] <- res4.9[2]
res4.10 <- GetPWFstHaplo(m4,m10,boot4.10)
PWFstMatrix[4,10] <- res4.10[1]
PWFstMatrix[10,4] <- res4.10[2]
res5.6 <- GetPWFstHaplo(m5,m6,boot5.6)
PWFstMatrix[5,6] <- res5.6[1]
PWFstMatrix[6,5] <- res5.6[2]
res5.7 <- GetPWFstHaplo(m5,m7,boot5.7)
PWFstMatrix[5,7] <- res5.7[1]
PWFstMatrix[7,5] <- res5.7[2]
res5.8 <- GetPWFstHaplo(m5,m8,boot5.8)
PWFstMatrix[5,8] <- res5.8[1]
PWFstMatrix[8,5] <- res5.8[2]
res5.9 <- GetPWFstHaplo(m5,m9,boot5.9)
PWFstMatrix[5,9] <- res5.9[1]
PWFstMatrix[9,5] <- res5.9[2]
res5.10 <- GetPWFstHaplo(m5,m10,boot5.10)
PWFstMatrix[5,10] <- res5.10[1]
PWFstMatrix[10,5] <- res5.10[2]
res6.7 <- GetPWFstHaplo(m6,m7,boot6.7)
PWFstMatrix[6,7] <- res6.7[1]
PWFstMatrix[7,6] <- res6.7[2]
res6.8 <- GetPWFstHaplo(m6,m8,boot6.8)
PWFstMatrix[6,8] <- res6.8[1]
PWFstMatrix[8,6] <- res6.8[2]
res6.9 <- GetPWFstHaplo(m6,m9,boot6.9)
PWFstMatrix[6,9] <- res6.9[1]
PWFstMatrix[9,6] <- res6.9[2]
res6.10 <- GetPWFstHaplo(m6,m10,boot6.10)
PWFstMatrix[6,10] <- res6.10[1]
PWFstMatrix[10,6] <- res6.10[2]
res7.8 <- GetPWFstHaplo(m7,m8,boot7.8)
PWFstMatrix[7,8] <- res7.8[1]
PWFstMatrix[8,7] <- res7.8[2]
res7.9 <- GetPWFstHaplo(m7,m9,boot7.9)
PWFstMatrix[7,9] <- res7.9[1]
PWFstMatrix[9,7] <- res7.9[2]
res7.10 <- GetPWFstHaplo(m7,m10,boot7.10)
PWFstMatrix[7,10] <- res7.10[1]
PWFstMatrix[10,7] <- res7.10[2]
res8.9 <- GetPWFstHaplo(m8,m9,boot8.9)
PWFstMatrix[8,9] <- res8.9[1]
PWFstMatrix[9,8] <- res8.9[2]
res8.10 <- GetPWFstHaplo(m8,m10,boot8.10)
PWFstMatrix[8,10] <- res8.10[1]
PWFstMatrix[10,8] <- res8.10[2]
res9.10 <- GetPWFstHaplo(m9,m10,boot9.10)
PWFstMatrix[9,10] <- res9.10[1]
PWFstMatrix[10,9] <- res9.10[2]

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
