#' Calculates differential correlation statistics for all the possible variable-variable combinations within a single dataset.
#'
#' Creates a table of pairwise correlations statistics for all variables within a dataframe, with separation of two groups for comparison. This can be used as an exploratory tool to investigate differential correlations within a dataset.
#' Requires grouping variable to compare two groups.
#'
#' @param data dataframe for analysis
#' @param group character string defining the grouping variable for comparative differential correlations
#' @param ordered character string defining the variable to order output by. Choose from \code{g1cor}, \code{g1p}, \code{g2cor}, \code{g2r}, \code{fisher} (default) and \code{BH}.
#' @param limit numeric input to limit number of output correlation pairings
#'
#' @return a table (or dataframe) with correlation coefficients, p-values, fisher r-to-z p-value and adjusted p-value using Benjamini-Hochberg corretion for each correlation pair, with separation of two groups for comparison.
#'
#' @author Emily Mears, \email{mears.emilyrose@gmail.com}, Matthew Grant, \email{mgra576@aucklanduni.ac.nz}
#' @author Ben Day, \email{benjamindayengineer@gmail.com}
#'
#' @examples
#'## Load example dataframes
#'df <- read.csv("example_data/excorr_df2.csv")
#'
#'## Run function
#'data.cor(df, group = "sex")
#'data.cor(df, group = "sex", ordered = "fisher", limit = 20)
#'
#' @export
#'



data.cor <-  function(data, group, ordered = "fisher", limit = NA){

  cor.prob <- function(X, dfr = nrow(X) -2){
    R <- cor(X, use = "pairwise.complete.obs")
    above <- row(R) < col(R)
    R2 <- R[above]^2
    Fstat <- R2 * dfr/(1 - R2)
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

  # Check for all-zero variables/cols
  if (TRUE %in% lapply(colnames(data), function(x) all(data[, x] == 0)) == TRUE) {
    data <- data[, colSums(data != 0) > 0]
    warning(paste("Some columns contain only zeros and were removed before calculation."))
  }

  # Check if group column contains more than 2 groups
  if(length(unique(as.character(data[data[,group] != as.character(data[,group][1]),group]))) > 1) {
    stop("More than 2 groups in grouping variable")}

  data[,which(colnames(data) == group)] = as.character(data[,which(colnames(data) == group)])

  group_1_name = as.character(data[,group][1])
  group_2_name = as.character(data[data[,group] != as.character(data[,group][1]),group][1])

  group1Cors = flattenSquareMatrix(cor.prob(data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], lapply(data, class) == "numeric"]))
  group2Cors = flattenSquareMatrix(cor.prob(data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], lapply(data, class) == "numeric"]))
  group1and2Cors = cbind(group1Cors, group2Cors[3:4])



  fisher_r_to_z = as.array(compcorr(n1 = as.numeric(table(data[,which(colnames(data) == group)])[1]),
                                    r1 = group1and2Cors[,3],
                                    n2 = as.numeric(table(data[,which(colnames(data) == group)])[2]),
                                    r2 = group1and2Cors[,5]))$pval

  group1and2compcor = cbind(group1and2Cors, fisher_r_to_z)
  BH_p_adjust <- p.adjust(fisher_r_to_z, method = "BH", n = length(fisher_r_to_z))
  group1and2compcor <- cbind(group1and2compcor, BH_p_adjust)

  colnames(group1and2compcor) = c("Var1",
                               "Var2",
                               paste(group_1_name, "correlation coefficent"),
                               paste(group_1_name, "correlation p-value"),
                               paste(group_2_name, "correlation coefficent"),
                               paste(group_2_name, "correlation p-value"),
                               "Fisher R to Z p-value",
                               "BH p-value adjustment")



  # # Rounding correlation coefficients
  group1and2compcor[,3] <- round(group1and2compcor[,3], digits = 2)
  group1and2compcor[,5] <- round(group1and2compcor[,5], digits = 2)
  # # Significant figures p-values
  group1and2compcor[,4] <- signif(group1and2compcor[,4], digits = 4)
  group1and2compcor[,6] <- signif(group1and2compcor[,6], digits = 4)
  group1and2compcor[,7] <- signif(group1and2compcor[,7], digits = 4)
  group1and2compcor[,8] <- signif(group1and2compcor[,8], digits = 4)

  #group1and2compcor[,3:8] <- as.numeric(group1and2compcor[,3:8])


  #Ordering options
  if(ordered == "g1cor"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,3], decreasing=FALSE),]
  }

  else if(ordered == "g1p"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,4], decreasing=FALSE),]
  }

  else if(ordered == "g2cor"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,5], decreasing=FALSE),]
  }

  else if(ordered == "g2p"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,6], decreasing=FALSE),]
  }

  else if(ordered == "fisher"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,7], decreasing=FALSE),]
  }

  else if(ordered == "BH"){
    group1and2compcor = group1and2compcor[order(group1and2compcor[,8], decreasing=FALSE),]
  }

  else {
    group1and2compcor = group1and2compcor[order(group1and2compcor[,7], decreasing=FALSE),]
    warning("Undefined ordering column.")
  }

  # Last step is to apply limit
  if(is.na(limit)){
    return(group1and2compcor)
  }

  else if(!is.na(limit)){
    return(group1and2compcor[1:limit,])
  }
}
