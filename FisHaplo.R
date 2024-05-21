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
library(doParallel)
library(foreach)
library(iterators)
m <- as.matrix(read.table(args[1], header=TRUE, colClasses = "character"))

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
      } else if ((m[locus,ind] == "01") || (m[locus,ind] == "10")) {
        a <- a + 1
        b <- b + 1
        Hobs <- Hobs + 1
        tot <- tot + 2   
      } else if ((m[locus,ind] == "02") || (m[locus,ind] == "20")) {
        a <- a + 1
        c <- c + 1
        Hobs <- Hobs + 1
        tot <- tot + 2  
      } else if ((m[locus,ind] == "12") || (m[locus,ind] == "21")) {
        b <- b + 1
        c <- c + 1
        Hobs <- Hobs + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "22")) {
        c <- c + 2
        C <- C + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "03") || (m[locus,ind] == "30")) {
        a <- a + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2     
      } else if ((m[locus,ind] == "13") || (m[locus,ind] == "31")) {
        b <- b + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2    
      } else if ((m[locus,ind] == "23") || (m[locus,ind] == "32")) {
        c <- c + 1
        d <- d + 1
        Hobs <- Hobs + 1
        tot <- tot + 2   
      } else if ((m[locus,ind] == "33")) {
        d <- d + 2
        D <- D + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "04") || (m[locus,ind] == "40")) {
        a <- a + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2   
      } else if ((m[locus,ind] == "14") || (m[locus,ind] == "41")) {
        b <- b + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "24") || (m[locus,ind] == "42")) {
        c <- c + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "34") || (m[locus,ind] == "43")) {
        d <- d + 1
        e <- e + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "44")) {
        e <- e + 2
        E <- E + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "05") || (m[locus,ind] == "50")) {
        a <- a + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2      
      } else if ((m[locus,ind] == "15") || (m[locus,ind] == "51")) {
        b <- b + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2      
      } else if ((m[locus,ind] == "25") || (m[locus,ind] == "52")) {
        c <- c + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2 
      } else if ((m[locus,ind] == "35") || (m[locus,ind] == "53")) {
        d <- d + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2     
      } else if ((m[locus,ind] == "45") || (m[locus,ind] == "54")) {
        e <- e + 1
        f <- f + 1
        Hobs <- Hobs + 1
        tot <- tot + 2   
      } else if ((m[locus,ind] == "55")) {
        f <- f + 2
        FF <- FF + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "06") || (m[locus,ind] == "60")) {
        a <- a + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2       
      } else if ((m[locus,ind] == "16") || (m[locus,ind] == "61")) {
        b <- b + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2       
      } else if ((m[locus,ind] == "26") || (m[locus,ind] == "62")) {
        c <- c + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2     
      } else if ((m[locus,ind] == "36") || (m[locus,ind] == "63")) {
        d <- d + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2       
      } else if ((m[locus,ind] == "46") || (m[locus,ind] == "64")) {
        e <- e + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2      
      } else if ((m[locus,ind] == "56") || (m[locus,ind] == "65")) {
        f <- f + 1
        g <- g + 1
        Hobs <- Hobs + 1
        tot <- tot + 2        
      } else if ((m[locus,ind] == "66")) {
        g <- g + 2
        G <- G + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "07") || (m[locus,ind] == "70")) {
        a <- a + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "17") || (m[locus,ind] == "71")) {
        b <- b + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "27") || (m[locus,ind] == "72")) {
        c <- c + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "37") || (m[locus,ind] == "73")) {
        d <- d + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "47") || (m[locus,ind] == "74")) {
        e <- e + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "57") || (m[locus,ind] == "75")) {
        f <- f + 1
        h <- h + 1
        Hobs <- Hobs + 1
        tot <- tot + 2
      } else if ((m[locus,ind] == "67") || (m[locus,ind] == "76")) {
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

m1_dat <- getFrqHaplo(m1)
m2_dat <- getFrqHaplo(m2)
m3_dat <- getFrqHaplo(m3)
m4_dat <- getFrqHaplo(m4)
m5_dat <- getFrqHaplo(m5)
m6_dat <- getFrqHaplo(m6)
m7_dat <- getFrqHaplo(m7)
m8_dat <- getFrqHaplo(m8)
m9_dat <- getFrqHaplo(m9)
m10_dat <- getFrqHaplo(m10)


registerDoParallel(32)
permutHaplo <- function(gametes) {
  R <- 10000
  boot <- c()
  
  shuffleGametes <- function(gametes){
    shuffledGametes <- matrix(nrow = nrow(gametes), ncol = ncol(gametes))
    sectionEnd <- 0
    for (i in 1 :nrow(gametes)) {
      section <- gametes[i,]
      for (j in 1:(ncol(gametes)/2)){
        randomNumber <- sample(c(0,1),1)
        if (randomNumber) {
          tmp <- section[2*j]
          section[2*j] <- section[2*j-1]
          section[2*j-1] <- tmp
        }
      }
      shuffledGametes[i,] <-  section
    }
    return(shuffledGametes)
  }
  
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
        } else if ((m[locus,ind] == "01") || (m[locus,ind] == "10")) {
          a <- a + 1
          b <- b + 1
          Hobs <- Hobs + 1
          tot <- tot + 2   
        } else if ((m[locus,ind] == "02") || (m[locus,ind] == "20")) {
          a <- a + 1
          c <- c + 1
          Hobs <- Hobs + 1
          tot <- tot + 2  
        } else if ((m[locus,ind] == "12") || (m[locus,ind] == "21")) {
          b <- b + 1
          c <- c + 1
          Hobs <- Hobs + 1
          tot <- tot + 2        
        } else if ((m[locus,ind] == "22")) {
          c <- c + 2
          C <- C + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "03") || (m[locus,ind] == "30")) {
          a <- a + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2     
        } else if ((m[locus,ind] == "13") || (m[locus,ind] == "31")) {
          b <- b + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2    
        } else if ((m[locus,ind] == "23") || (m[locus,ind] == "32")) {
          c <- c + 1
          d <- d + 1
          Hobs <- Hobs + 1
          tot <- tot + 2   
        } else if ((m[locus,ind] == "33")) {
          d <- d + 2
          D <- D + 1
          tot <- tot + 2        
        } else if ((m[locus,ind] == "04") || (m[locus,ind] == "40")) {
          a <- a + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2   
        } else if ((m[locus,ind] == "14") || (m[locus,ind] == "41")) {
          b <- b + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2        
        } else if ((m[locus,ind] == "24") || (m[locus,ind] == "42")) {
          c <- c + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2        
        } else if ((m[locus,ind] == "34") || (m[locus,ind] == "43")) {
          d <- d + 1
          e <- e + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "44")) {
          e <- e + 2
          E <- E + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "05") || (m[locus,ind] == "50")) {
          a <- a + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2      
        } else if ((m[locus,ind] == "15") || (m[locus,ind] == "51")) {
          b <- b + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2      
        } else if ((m[locus,ind] == "25") || (m[locus,ind] == "52")) {
          c <- c + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2 
        } else if ((m[locus,ind] == "35") || (m[locus,ind] == "53")) {
          d <- d + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2     
        } else if ((m[locus,ind] == "45") || (m[locus,ind] == "54")) {
          e <- e + 1
          f <- f + 1
          Hobs <- Hobs + 1
          tot <- tot + 2   
        } else if ((m[locus,ind] == "55")) {
          f <- f + 2
          FF <- FF + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "06") || (m[locus,ind] == "60")) {
          a <- a + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2       
        } else if ((m[locus,ind] == "16") || (m[locus,ind] == "61")) {
          b <- b + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2       
        } else if ((m[locus,ind] == "26") || (m[locus,ind] == "62")) {
          c <- c + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2     
        } else if ((m[locus,ind] == "36") || (m[locus,ind] == "63")) {
          d <- d + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2       
        } else if ((m[locus,ind] == "46") || (m[locus,ind] == "64")) {
          e <- e + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2      
        } else if ((m[locus,ind] == "56") || (m[locus,ind] == "65")) {
          f <- f + 1
          g <- g + 1
          Hobs <- Hobs + 1
          tot <- tot + 2        
        } else if ((m[locus,ind] == "66")) {
          g <- g + 2
          G <- G + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "07") || (m[locus,ind] == "70")) {
          a <- a + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "17") || (m[locus,ind] == "71")) {
          b <- b + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "27") || (m[locus,ind] == "72")) {
          c <- c + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "37") || (m[locus,ind] == "73")) {
          d <- d + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "47") || (m[locus,ind] == "74")) {
          e <- e + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "57") || (m[locus,ind] == "75")) {
          f <- f + 1
          h <- h + 1
          Hobs <- Hobs + 1
          tot <- tot + 2
        } else if ((m[locus,ind] == "67") || (m[locus,ind] == "76")) {
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
  
  shuffleColumns <- function(gametes){
    rand <- sample(c(1:ncol(gametes)), replace = F)
    gametes <- gametes[,rand]
  }
  
  gamsToInds <- function(gametes){
    M <- matrix(nrow=nrow(gametes), ncol=(ncol(gametes)/2))
    for (col in 1:ncol(M)){
      M[,col] <-  paste(gametes[,(col*2)], gametes[,(col*2-1)], sep ="")
    }
    for (col in 1:ncol(M)){
      for (row in 1:nrow(M)){
        if (nchar(M[row,col]) > 2) {
          M[row,col] <- NA
        }
      }
    }
    
    return(M)
  }

  
  boot <- foreach(icount(R), .combine = c, .inorder = F) %dopar% {
    
    recombined <- shuffleGametes(gametes)
    
    pool <- shuffleColumns(recombined)
    
    M <- gamsToInds(pool)
    
    M_dat <- getFrqHaplo(M)
    
    N <- mean(na.omit(M_dat[,17]))
    D <- mean(na.omit(M_dat[,18]))
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

boot1 <- permutHaplo(gametes1)
boot2 <- permutHaplo(gametes2)
boot3 <- permutHaplo(gametes3)
boot4 <- permutHaplo(gametes4)
boot5 <- permutHaplo(gametes5)
boot6 <- permutHaplo(gametes6)
boot7 <- permutHaplo(gametes7)
boot8 <- permutHaplo(gametes8)
boot9 <- permutHaplo(gametes9)
boot10 <- permutHaplo(gametes10)

FisTable <- matrix(nrow = 10, ncol = 4)
FisTable[1,1] <- 1 -(mean(m1_dat[,17])/(mean(m1_dat[,18])))
FisTable[2,1] <- 1 -(mean(m2_dat[,17])/(mean(m2_dat[,18])))
FisTable[3,1] <- 1 -(mean(m3_dat[,17])/(mean(m3_dat[,18])))
FisTable[4,1] <- 1 -(mean(m4_dat[,17])/(mean(m4_dat[,18])))
FisTable[5,1] <- 1 -(mean(m5_dat[,17])/(mean(m5_dat[,18])))
FisTable[6,1] <- 1 -(mean(m6_dat[,17])/(mean(m6_dat[,18])))
FisTable[7,1] <- 1 -(mean(m7_dat[,17])/(mean(m7_dat[,18])))
FisTable[8,1] <- 1 -(mean(m8_dat[,17])/(mean(m8_dat[,18])))
FisTable[9,1] <- 1 -(mean(m9_dat[,17])/(mean(m9_dat[,18])))
FisTable[10,1] <- 1 -(mean(m10_dat[,17])/(mean(m10_dat[,18])))

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
