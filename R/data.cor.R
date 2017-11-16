#' Differential correlations between all variables within a dataset
#'
#' @data the dataset to be used
#'
#' @group the variable to distinguish between correlation groups
#'
#' @odered which column should the table be ordered by?
#'
#' @limit how many combinations should be displayed?
#'
#' @author Emily Mears, \email{mears.emilyrose@gmail.com}
#'
#' @export
#'


data.cor <-  function(data, group, ordered = "None", limit = NA){
  
  cor.prob <- function(X, dfr = nrow(X) -2){
    R <- cor(X, use = "pairwise.complete.obs")
    above <- row(R) < col(R)
    r2 <- R[above]^2
    Fstat <- r2 * dfr/(1 - r2)
    R[above] <- 1 - pf(Fstat, 1, dfr)
    R[row(R) == col(R)] <-NA
    R
  }
  
  flattenSquareMatrix <- function(m) {
    if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.")
    if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
    ut <- upper.tri(m)
    data.frame(i = rownames(m)[row(m)[ut]],
               j = rownames(m)[col(m)[ut]],
               cor = t(m)[ut],
               p = m[ut])
  }
  
  compcorr <- function(n1, r1, n2, r2){
    # Fisher's Z-transformation
    # ad hoc process
    num1a <- which(r1 >= 0.99)
    num2a <- which(r2 >= 0.99)
    r1[num1a] <- 0.99
    r2[num2a] <- 0.99
    num1b <- which(r1 <= -0.99)
    num2b <- which(r2 <= -0.99)
    r1[num1b] <- -0.99
    r2[num2b] <- -0.99
    # z
    z1 <- atanh(r1)
    z2 <- atanh(r2)
    # difference
    dz <- (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
    # p-value
    pv <- 2*( 1 - pnorm(abs(dz)) )
    return(list(diff=dz, pval=pv))
  }
  
  my_compcorr <- function(x){
    compcorr(n1 = 6, r1 = x[,3], n2 = 6, r2 = x[,5])$pval
  }
  
  data[,which(colnames(data) == group)] = as.character(data[,which(colnames(data) == group)])
  
  group_1_name = as.character(data[,group][1])
  group_2_name = as.character(data[data[,group] != as.character(data[,group][1]),group][1])
  
  data[which(as.character(data[,which(colnames(data) == group)]) == as.character(data[,which(colnames(data) == group)])[1]), which(colnames(data) == group)] = "1"
  data[which(as.character(data[,which(colnames(data) == group)]) != as.character(data[,which(colnames(data) == group)])[1]), which(colnames(data) == group)] = "2"
  
  group1Cors = flattenSquareMatrix(cor.prob(data[data[,group] == 1, lapply(data, class) == "numeric"]))
  group2Cors = flattenSquareMatrix(cor.prob(data[data[,group] == 2, lapply(data, class) == "numeric"]))
  group1and2Cors = cbind(group1Cors, group2Cors[3:4])
  
  colnames(group1and2Cors) = c("Var1", "Var2", paste(group_1_name, "Correlation"), paste(group_1_name, "P-value"), paste(group_2_name, "Correlation"), paste(group_2_name, "P-value"))
  Fisher_Value = as.array(my_compcorr(group1and2Cors))
  
  group1and2compcor = cbind(group1and2Cors, Fisher_Value)
  BH_p_adjust <- p.adjust(Fisher_Value, method = "BH", n = length(Fisher_Value))
  group1and2compcor <- cbind(group1and2compcor, BH_p_adjust)
  
  #Ordering options
  if(ordered == "None"){
    group1and2compcor = group1and2compcor
  }
  
  else if(ordered == "Group_1_correlation"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,3],decreasing=FALSE),]
  }
  
  else if(ordered == "Group_1_Pvalue"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,4],decreasing=FALSE),]
  }
  
  else if(ordered == "Group_2_correlation"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,5],decreasing=FALSE),]
  }
  
  else if(ordered == "Group_2_Pvalue"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,6],decreasing=FALSE),]
  }
  
  else if(ordered == "Fisher_value"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,7],decreasing=FALSE),]
  }
  
  else if(ordered == "BH_p_adjust"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,8],decreasing=FALSE),]
  }
  
  if(is.na(limit)){
    return(group1and2compcor)
  }
  
  else if(!is.na(limit)){
    return(group1and2compcor[1:limit,])
  } 
}