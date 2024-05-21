library(vcfR)
library(adegenet)

vcf <- read.vcfR("sacharina_latissima_151_mac3_haplotype_filtered.vcf")

genind <- vcfR2genind(vcf)
df <- genind2df(genind)
df <- t(df)
loci <- c(rownames(df))
m <- as.matrix(df)
class(m) <- "numeric"
m <- m[,-55]


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

# MAC filtering ####

mac <- 3
remove <- c()

for (row in 1: nrow(m)){
  if (sum(na.omit(m[row,])) < mac | (2*length(na.omit(m[row,])) - sum(na.omit(m[row,]))) < mac){
    remove <- append(remove, row)
  }
}

print(gettextf("Identified %i SNP loci with mac under %i.", length(remove), mac))

m <- m[-remove,]
loci <- loci[-remove]

#### Remove low completeness loci####

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

m_dat <- getFrq(m)
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


compfilter <- function(m, m_dat, min) {
  remove <- c()
  for (row in 1:nrow(m)) {
    # if (is.na(m_dat[row,6])){
    #   remove <- append(remove, row)
    # }
    if (ncol(m)-(m_dat[row,7]) < min) {
      remove <- append(remove, row)
    }
  }
  return(remove)
}

NAtot <- 120

NApop <- 10

totalList <- compfilter(m, m_dat, NAtot)
acrossPop <- length(totalList)
totalList <- append(totalList, compfilter(m1, m1_dat, NApop))
totalList <- append(totalList, compfilter(m2, m2_dat, NApop))
totalList <- append(totalList, compfilter(m3, m3_dat, NApop))
totalList <- append(totalList, compfilter(m4, m4_dat, NApop))
totalList <- append(totalList, compfilter(m5, m5_dat, NApop))
totalList <- append(totalList, compfilter(m6, m6_dat, NApop))
totalList <- append(totalList, compfilter(m7, m7_dat, NApop))
totalList <- append(totalList, compfilter(m8, m8_dat, NApop))
totalList <- append(totalList, compfilter(m9, m9_dat, NApop))
totalList <- append(totalList, compfilter(m10, m10_dat, NApop))

totalList <- unique(totalList)

print(gettextf("Identified %i loci that are present in less than %i individuals of the total population.", acrossPop, NAtot))
print(gettextf("Identified %i loci that are present in less than than %i individuals in any of the subpopulations.", (length(totalList)-acrossPop), NApop))
print(gettextf("A total of %i loci has been removed.", length(totalList)))

m <- m[-totalList,]
loci <- loci[-totalList]
dim(m)

table <- cbind(loci,m)
#write.table(table,file="SNPsMac3.txt",row.names=FALSE) 