library (readr)

df <- read_tsv(file = 'haplotypes.tsv')
loci <- df[,c(1:2)]
m <- as.matrix(df[,-c(1:2)])
class(m) <- "numeric"
# m <-m[,-55] # remove duplicate sample


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

# Completeness filtering

compFilter <- function(min, m, remove){
  for (row in 1:nrow(m)){
    if (length(na.omit(m[row,])) < min){
      remove <- append(remove, row)
    }
  }
  remove <- unique(remove)
  print(gettextf("Identified %i alleles from low completeness loci.", length(remove)))
  return(remove)
}

min <- 10
remove <- compFilter(0.8*150, m, c())

remove <- append(remove,compFilter(min, m1, remove))
remove <- append(remove,compFilter(min, m2, remove))
remove <- append(remove,compFilter(min, m3, remove))
remove <- append(remove,compFilter(min, m4, remove))
remove <- append(remove,compFilter(min, m5, remove))
remove <- append(remove,compFilter(min, m6, remove))
remove <- append(remove,compFilter(min, m7, remove))
remove <- append(remove,compFilter(min, m8, remove))
remove <- append(remove,compFilter(min, m9, remove))
remove <- append(remove,compFilter(min, m10, remove))

m <- m[-remove,]
loci <- loci[-remove,]

# MAC filtering

mac <- 3
remove <- c()
for (row in 1:nrow(m)){
  count <- sum(na.omit(m[row,]))
  if (count < mac){
    remove <- append(remove, row)
  }
}
print(gettextf("Identified %i low frequency alleles.", length(remove)))
m <- m[-remove,]
loci <- loci[-remove,]

# Remove non-polymorphic loci

remove <- c()
for (row in 2:(nrow(loci)-1)){
  if (loci[row,1] != loci[row-1,1] && loci[row,1] != loci[row+1,1]){
    remove <- append(remove,row)
  }
}
if (loci[nrow(loci),1] != loci[nrow(loci)-1,1]){
  remove <- append(remove,nrow(loci))
}
print(gettextf("Identified %i non-polymorphic loci.", length(remove)))
#remove
m <- m[-remove,]
loci <- loci[-remove,]
dim(m)
length(unlist(unique(loci[,1])))

table <- cbind(loci,m)

#write.table(table,file="haplotypesMac3.txt",row.names=FALSE) 

