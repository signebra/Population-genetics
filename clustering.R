library(dendextend)
G <- read.table(file = "GSNP3.txt")
G <- G[-1,]
D <- dist(G)
hclust.complete <- hclust(D, method = "complete")
plot(hclust.complete)

#### Determine K ####
  
  hclust.complete <- hclust(D, method = "complete")
  
  clusters <- cutree(hclust(D, method = "complete"), k = K)
  
  dist <- as.matrix(D)
  for (i in 1:nrow(dist)){
    dist[i,i] <- NA
  }
  
  getMeanS <- function(clusters,K,dist){
    
    getA <- function(clusters,ind,dist){
      c <- clusters[ind]
      C <- (clusters == c)
      A <- dist[ind,C]
      a <- mean(na.omit(A))
      return(a)
    }
  
    getB <- function(clusters,ind,K,dist){
      b <- Inf
      for (c in 1:K){
        if (c != clusters[ind]){
          C <- (clusters == c)
          B <- dist[ind,C]
          temp <- mean(B)
          if (temp < b){
            b <- temp
          }
        }
      }
      return(b)
    }
    
    tot <- 0
    for (ind in 1:150){
      if (sum(clusters == clusters[ind]) > 1){
        a <- getA(clusters,ind,dist)
        b <- getB(clusters,ind,K,dist)
        S <- (b-a)/max(b,a)
        tot <- tot + S
      }
    }
    meanS <- tot/150
    return(meanS)
  }

  silhouettes <- c()
  for (K in 2:15){
    clusters <- cutree(hclust(D, method = "complete"), k = K)
    silhouettes <- append(silhouettes, getMeanS(clusters,K,dist))
  }
  silhouettes
  plot(silhouettes, type="l")

#### Make cluster tables ####

# SNPS:
# mac3 -> K = 3
# mac15 -> K = 7 (inkl. 1 outlier)
# Haplotypes:
# mac3 -> K = 6
# mac15 -> K = 8

G <- read.table(file = "G.txt")
G <- G[-1,]
D <- dist(G)
K = 5
clusters <- cutree(hclust(D, method = "complete"), k = K)

table <- matrix(ncol = 3, nrow = 150)
table[,1] <- clusters
table[,2] <- c(1:150)
table[,3] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
               4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
               5,5,5,5,5,5,5,5,5,5,5,5,5,5,
               6,6,6,6,6,6,6,6,6,6,6,6,6,6,
               7,7,7,7,7,7,7,7,7,7,7,7,7,7,
               8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
               9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
               10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
table <- table[order(table[,1],decreasing=FALSE),]
table1 <- table

table2 <- matrix(nrow = K, ncol = 10)
all <- c(1:150)
for (cluster in 1:K){
  p1 <- 0
  p2 <- 0
  p3 <- 0
  p4 <- 0
  p5 <- 0
  p6 <- 0
  p7 <- 0
  p8 <- 0
  p9 <- 0
  p10 <- 0
  for (row in as.vector(all[(table1[,1] == cluster)])) {
    if (table1[row,3] == 1) {
      p1 <- p1 + 1
    }
    if (table1[row,3] == 2) {
      p2 <- p2 + 1
    }
    if (table1[row,3] == 3) {
      p3 <- p3 + 1
    }
    if (table1[row,3] == 4) {
      p4 <- p4 + 1
    }
    if (table1[row,3] == 5) {
      p5 <- p5 + 1
    }
    if (table1[row,3] == 6) {
      p6 <- p6 + 1
    }
    if (table1[row,3] == 7) {
      p7 <- p7 + 1
    }
    if (table1[row,3] == 8) {
      p8 <- p8 + 1
    }
    if (table1[row,3] == 9) {
      p9 <- p9 + 1
    }
    if (table1[row,3] == 10) {
      p10 <- p10 + 1
    }
  }
  table2[cluster, 1] <- p1
  table2[cluster, 2] <- p2
  table2[cluster, 3] <- p3
  table2[cluster, 4] <- p4
  table2[cluster, 5] <- p5
  table2[cluster, 6] <- p6
  table2[cluster, 7] <- p7
  table2[cluster, 8] <- p8
  table2[cluster, 9] <- p9
  table2[cluster, 10] <- p10

}
table2
#write.table(table2, file = "ClustersSNP3.txt")

G <- read.table(file = "GSNP15.txt")
D <- dist(G)
K = 8
clusters <- cutree(hclust(D, method = "complete"), k = K)

table1 <- matrix(ncol = 3, nrow = 150)
table1[,1] <- clusters
table1[,2] <- c(1:150)
table1[,3] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
               4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
               5,5,5,5,5,5,5,5,5,5,5,5,5,5,
               6,6,6,6,6,6,6,6,6,6,6,6,6,6,
               7,7,7,7,7,7,7,7,7,7,7,7,7,7,
               8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
               9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
               10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
table1 <- table1[order(table1[,1],decreasing=FALSE),]
table <- cbind (table, table1)

table2 <- matrix(nrow = K, ncol = 10)
all <- c(1:150)
for (cluster in 1:K){
  p1 <- 0
  p2 <- 0
  p3 <- 0
  p4 <- 0
  p5 <- 0
  p6 <- 0
  p7 <- 0
  p8 <- 0
  p9 <- 0
  p10 <- 0
  for (row in as.vector(all[(table1[,1] == cluster)])) {
    if (table1[row,3] == 1) {
      p1 <- p1 + 1
    }
    if (table1[row,3] == 2) {
      p2 <- p2 + 1
    }
    if (table1[row,3] == 3) {
      p3 <- p3 + 1
    }
    if (table1[row,3] == 4) {
      p4 <- p4 + 1
    }
    if (table1[row,3] == 5) {
      p5 <- p5 + 1
    }
    if (table1[row,3] == 6) {
      p6 <- p6 + 1
    }
    if (table1[row,3] == 7) {
      p7 <- p7 + 1
    }
    if (table1[row,3] == 8) {
      p8 <- p8 + 1
    }
    if (table1[row,3] == 9) {
      p9 <- p9 + 1
    }
    if (table1[row,3] == 10) {
      p10 <- p10 + 1
    }
  }
  table2[cluster, 1] <- p1
  table2[cluster, 2] <- p2
  table2[cluster, 3] <- p3
  table2[cluster, 4] <- p4
  table2[cluster, 5] <- p5
  table2[cluster, 6] <- p6
  table2[cluster, 7] <- p7
  table2[cluster, 8] <- p8
  table2[cluster, 9] <- p9
  table2[cluster, 10] <- p10
  
}
write.table(table2, file = "ClustersSNP15.txt")

G <- read.table(file = "GHaplo3.txt")
D <- dist(G)
K = 6
clusters <- cutree(hclust(D, method = "complete"), k = K)

table1 <- matrix(ncol = 3, nrow = 150)
table1[,1] <- clusters
table1[,2] <- c(1:150)
table1[,3] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
                9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
                10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
table1 <- table1[order(table1[,1],decreasing=FALSE),]
table <- cbind (table, table1)

table2 <- matrix(nrow = K, ncol = 10)
all <- c(1:150)
for (cluster in 1:K){
  p1 <- 0
  p2 <- 0
  p3 <- 0
  p4 <- 0
  p5 <- 0
  p6 <- 0
  p7 <- 0
  p8 <- 0
  p9 <- 0
  p10 <- 0
  for (row in as.vector(all[(table1[,1] == cluster)])) {
    if (table1[row,3] == 1) {
      p1 <- p1 + 1
    }
    if (table1[row,3] == 2) {
      p2 <- p2 + 1
    }
    if (table1[row,3] == 3) {
      p3 <- p3 + 1
    }
    if (table1[row,3] == 4) {
      p4 <- p4 + 1
    }
    if (table1[row,3] == 5) {
      p5 <- p5 + 1
    }
    if (table1[row,3] == 6) {
      p6 <- p6 + 1
    }
    if (table1[row,3] == 7) {
      p7 <- p7 + 1
    }
    if (table1[row,3] == 8) {
      p8 <- p8 + 1
    }
    if (table1[row,3] == 9) {
      p9 <- p9 + 1
    }
    if (table1[row,3] == 10) {
      p10 <- p10 + 1
    }
  }
  table2[cluster, 1] <- p1
  table2[cluster, 2] <- p2
  table2[cluster, 3] <- p3
  table2[cluster, 4] <- p4
  table2[cluster, 5] <- p5
  table2[cluster, 6] <- p6
  table2[cluster, 7] <- p7
  table2[cluster, 8] <- p8
  table2[cluster, 9] <- p9
  table2[cluster, 10] <- p10
  
}
write.table(table2, file = "ClustersHaplo3.txt")

G <- read.table(file = "GHaplo15.txt")
D <- dist(G)
K = 8
clusters <- cutree(hclust(D, method = "complete"), k = K)

table1 <- matrix(ncol = 3, nrow = 150)
table1[,1] <- clusters
table1[,2] <- c(1:150)
table1[,3] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
                5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                7,7,7,7,7,7,7,7,7,7,7,7,7,7,
                8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
                9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,
                10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
table1 <- table1[order(table1[,1],decreasing=FALSE),]
table <- cbind (table, table1)

table2 <- matrix(nrow = K, ncol = 10)
all <- c(1:150)
for (cluster in 1:K){
  p1 <- 0
  p2 <- 0
  p3 <- 0
  p4 <- 0
  p5 <- 0
  p6 <- 0
  p7 <- 0
  p8 <- 0
  p9 <- 0
  p10 <- 0
  for (row in as.vector(all[(table1[,1] == cluster)])) {
    if (table1[row,3] == 1) {
      p1 <- p1 + 1
    }
    if (table1[row,3] == 2) {
      p2 <- p2 + 1
    }
    if (table1[row,3] == 3) {
      p3 <- p3 + 1
    }
    if (table1[row,3] == 4) {
      p4 <- p4 + 1
    }
    if (table1[row,3] == 5) {
      p5 <- p5 + 1
    }
    if (table1[row,3] == 6) {
      p6 <- p6 + 1
    }
    if (table1[row,3] == 7) {
      p7 <- p7 + 1
    }
    if (table1[row,3] == 8) {
      p8 <- p8 + 1
    }
    if (table1[row,3] == 9) {
      p9 <- p9 + 1
    }
    if (table1[row,3] == 10) {
      p10 <- p10 + 1
    }
  }
  table2[cluster, 1] <- p1
  table2[cluster, 2] <- p2
  table2[cluster, 3] <- p3
  table2[cluster, 4] <- p4
  table2[cluster, 5] <- p5
  table2[cluster, 6] <- p6
  table2[cluster, 7] <- p7
  table2[cluster, 8] <- p8
  table2[cluster, 9] <- p9
  table2[cluster, 10] <- p10
  
}
write.table(table2, file = "ClustersHaplo15.txt")




