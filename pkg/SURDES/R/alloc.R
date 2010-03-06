alloc <- function(mdist, vini=rep(0, dim(mdist)[2]), vtarget=rep(1, dim(mdist)[2]), criteria, sdint=rep(1,length(criteria)), conditions, iter){

  mdistmatrix <- as.matrix(mdist)
  if (any(diag(mdistmatrix) != 0) | any(mdistmatrix[lower.tri(t(mdistmatrix))] != mdistmatrix[lower.tri(mdistmatrix)])){
    stop('mdist must be a distance matrix')
  }

  if (length(vini) != (dim(mdist)[2]) | (length(vtarget)) != (dim(mdist)[2]) | (dim(conditions)[2]) != (dim(mdist)[2])){
   stop('Dimensions of input matrices and vectors do not match')
  }

  if ((length(vtarget[vtarget != 0]) - (length(vini[vini != 0]))) < iter){
    stop('Number of iteractions greater than available sites')
    }

  pmmatrix <- matrix(NA, ncol=dim(mdist)[1], nrow=iter)
  selmatrix <- matrix(NA, ncol=dim(mdist)[1], nrow=iter)
  selNames <- vector(length=iter)

  fmin <- function(x,n){
  n=0
    ifelse(x == min(x, na.rm=TRUE), 1, 0)
  }

  fmax <- function(x,n){
    n=0
    ifelse(x==max(x, na.rm=TRUE), 1, 0)
  }

  rangemax <- function(x,n){
    sd <- n*sd(x, na.rm=TRUE)
    if (is.na(sd)) sd <- 0
    max <- max(x, na.rm=TRUE)
    if (is.na(sd)) x[is.na(x)] <- 0
    if (is.na(sd)) x[x > 0] <- 1 else
      ifelse(x >= (max - sd), 1, 0)
  }

  rangemin <- function(x,n){
    sd <- n*sd(x, na.rm=TRUE)
    if (is.na(sd)) sd <- 0
    min <- min(x, na.rm=TRUE)
    if (is.na(sd)) x[is.na(x)] <- 0
    if(is.na(sd)) x[x > 0] <- 1 else
      ifelse(x <= (min + sd), 1, 0)
  }

  cond <- function(x, ...) UseMethod("cond")
  cond.max <- function(x,n) fmax(x,n)
  cond.min <- function(x,n) fmin(x,n)
  cond.rangemax <- function(x,n) rangemax(x,n)
  cond.rangemin <- function(x,n) rangemin(x,n)
  cond.default <- function(x) print("Criterion not specified")

  conditionsSplit <- split(conditions, row(conditions))
  criteriaList <- as.list(criteria)

  for(j in 1:length(criteria)){
    if (criteriaList[j] == "max") attr(conditionsSplit[[j]], "class") <- "max" else
    if (criteriaList[j] == "min") attr(conditionsSplit[[j]], "class") <- "min" else
    if (criteriaList[j] == "rangemax") attr(conditionsSplit[[j]], "class") <- "rangemax" else
    if (criteriaList[j] == "rangemin") attr(conditionsSplit[[j]], "class") <- "rangemin" else
    print("Criterion not specified")
  }

  for (i in 1:iter){

    prueba2 <- apply(mdist, 1, function(x){
      ifelse(vini == 1, x - 0, 1)})
    f <- apply(prueba2, 2, function(x) min(x))
    prueba2 <- apply(mdist, 2, function(x){
      ifelse(x < f, x - 0, f)})

    pmedian <- colSums(as.matrix(prueba2))
    pmedian <- ifelse(rowSums(prueba2) == 0, 0, pmedian)
    pmmatrix[i, ] <- pmedian
    mini <- pmedian * ifelse(vtarget < 1, 0, 1)
    min.sd <- min(mini[mini != 0]) + sd(mini[mini != 0])
    vselectedE <- ifelse(mini == 0 | mini > min.sd, 0, 1)


    vselectedP <- matrix(NA, ncol=ncol(conditions), nrow=(nrow(conditions)) + 1)
    names1 <- c(rep("crit", length(criteria) + 1))
    names2 <- c("pmedian", 1:length(criteria))
    names <- paste(names1, names2, sep="")
    rownames(vselectedP) <- names
    vselectedP[1, ] <- vselectedE

    for(l in 1:length(criteria)){
      a <- conditionsSplit[[l]]
      b <- vselectedP[l, ]
      a[which(b == 0)] <- NA
      c <- cond(a,sdint[l])
      c[is.na(c)] <- 0
      vselectedP[l + 1, ] <- c
      selection <- vselectedP[l + 1, ]
    }


    if (any(rowSums(vselectedP, na.rm=TRUE) == 1)){
      print(paste('Variable selected is: ', min(names(which(rowSums(vselectedP, na.rm=TRUE) == 1))), sep=''))
    } else{
      print('No variable found. Random variable selected')
    }

    f5a <- runif(length(selection), 1, 10) * selection
    selection <- ifelse(f5a == max(f5a), 1, 0)
    selmatrix[i, ] <- selection
    vini <- vini + selection
  }

  result <- list(selmatrix, pmmatrix)
  names(result) <- c("selmatrix", "pmmatrix")
  return(result)
}

