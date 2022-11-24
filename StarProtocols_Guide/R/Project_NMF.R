# Create Project.NMF function
project.NMF <- function(input.array, nmf.result){
  require(MASS)
  m <- as.matrix(input.array)
  gs.names <- rownames(input.array)
  gs.descs <- rownames(input.array)
  sample.names <- colnames(input.array)
  
  W <- as.data.frame(basis(nmf.result))
  W.row.names <- row.names(W)
  W.row.descs <- row.names(W)
  W.names <- names(W)
  
  overlap <- intersect(gs.names, W.row.names)
  
  cat("Size of Input dataset:", length(gs.names), "genes\n")
  cat("Size of W matrix (rows):", length(W.row.names), "genes\n")
  cat("Size of overlap:", length(overlap), "genes\n")
  
  locations.m <- match(overlap, gs.names, nomatch=0)
  m2 <- m[locations.m, ]
  locations.W <- match(overlap, W.row.names, nomatch=0)
  W2 <- W[locations.W, ]
  W2 <- as.matrix(W2)
  
  H <- ginv(W2) %*% m2
  
  n.col <- length(H[1,])
  for (i in 1:n.col) {
    S.2 <- sqrt(sum(H[,i]*H[,i]))
    #        S.2 <- sum(H[,i])
    H[,i] <- H[,i]/S.2
  }
  
  V <- data.frame(H)
  names(V) <- sample.names
  row.names(V) <- W.names
  return(V)
}
