#m= # snps
#n= # inds
## assuming all SNPs have allele frequency of 50%
sim <- function(m, n, rep = 1000) {
  nuniq <- 0
  nuniq2 <- 0
  first.match <- NULL;
  for(i in 1:rep) {
    mat1 <- matrix(rbinom(m*n, 1, 0.5), n, m);
    mat2 <- matrix(rbinom(m*n, 1, 0.5), n, m);

    reads <- NULL;

    for(j in 1:m) {
      if(sample(1,1:2)==1) {
        reads <- c(reads, mat1[1,j]);
      } else {
        reads <- c(reads, mat2[1,j]);
      }
    }

    match <- t(apply(mat1,1,"==",reads)) | t(apply(mat2,1,"==",reads));

    if ( nrow(unique(matrix(rbinom(m*n, 2, 0.5), n, m))) == n ) { nuniq <- nuniq + 1}

    if ( length(which(rowSums(match)==m)) == 1 ) {nuniq2 <- nuniq2+1;}

    first.match <- c(first.match, which(rowSums(match)==m)[[1]]);

    ##browser();
  }
  return (list(nuniq/rep, nuniq2/rep, first.match))
}

M <- seq(2,100,2);
N <- seq(2,64,2);

out <- matrix(NA, nrow=length(M), ncol=length(N))

for(m.i in 1:length(M)) {
  print(m.i);
  m <- M[m.i];
  for(n.i in 1:length(N)) {
    n <- N[n.i];
    out[m.i,n.i] <- sim(m,n)[[2]];
  }
}
