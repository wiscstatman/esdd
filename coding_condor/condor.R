#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
#library(partitions)
## idea

## for a putative gsk.set A = { j }, 

## 1. use available data on A^c to form clusters of targets
##    find a possible coding   
## 2. measure the quality of the gsk-set by homogeneity within
##    clusters on the activities of non-gsk compounds

## y_ij = activity indicator u[i,j] > threshold
## sum_{block in cluster} sum_{i,i' in block}  d_{i,i'}
##
## where d_{i,i'} = 1 -  sum_[ j in A^c ] [ y_ij * y_i'j ]/sum [ y_ij max y_i'j]

# note that in summing over pairs within blocks there's more things to
# sum (given number of clusters) in unbalanced compared to balanced blockings
# and the distances will be small if 

#y <- 1*(u > .5 )

target <- as.numeric(args[1])
nclust <- as.numeric(args[2])

load("pkis1.RData")

y <- y[-target,] ## remove the first row for cross validation
m <- nrow(y) ## number of targets
n <- ncol(y) ## number of compounds

source("imeasure3.R")

specs <- list( nclust=nclust, isize=16, phi=5, nbig=1000000 )

## let's do a randomized search

nbig <- specs$nbig 

nc <- specs$nc  ## number of clusters
isize <- specs$isize  ## gsk set size
phi <- specs$phi  ## penalty for extrapolation

Asum <- rep(0,n)
Anum <- rep(0,n)

a <- sample(1:n, size=isize )
A <- rep(FALSE, n)
A[a] <- TRUE
im <- imeasure3(A, y, nclust=nc )

Anum[a] <- 1

best <- vector( mode="list" )

best$A <- A
best$score <- im$sumd - phi*length(im$mapping) 
best$clust <- im$cluster
best$map <- im$mapping

Asum[a] <- best$score

for( ii in 2:nbig )
{
  a <- sample(1:n, size=isize )
  A <- rep(FALSE, n)
  A[a] <- TRUE
  
  im <- imeasure3(A, y, nclust=nc )
  score <- im$sumd - phi*length( im$mapping )  ## longer mapping's better
  Asum[a] <- Asum[a] + score ## override Michael's code
  Anum[a] <- Anum[a] +1
  
  if( score < best$score )
  {
    best$score <- score
    best$A <- A
    best$clust <- im$cluster
    best$map <- im$mapping
  }
  print(ii)
}


svname <- paste(target, nclust, sep="_")
svname <- paste("tar", svname, sep="")
svname <- paste(svname, "RData", sep=".")
## svname <- paste(target, svname, sep="/")

save( Asum, Anum, best, file=svname)

