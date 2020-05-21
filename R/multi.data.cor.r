#' Calculates differential correlation statistics for all the possible variable-variable combinations between two datasets.
#'
#' Requires grouping variable to contain 2 groups. Requires same samples across datasets (arranged in same order).
#' Declare columns (descriptors) column to both datasets.
#'
#' @param data1 first dataframe
#' @param data2 second dataframe
#' @param sample_col character string defining the sample or identifying column in both dataframes
#' @param group character string defining the grouping variable for comparative differential correlations
#' @param ordered character string defining the variable to order output by. Choose from \code{g1cor}, \code{g1p}, \code{g2cor}, \code{g2r}, \code{fisher} (default) and \code{BH}.
#' @param limit numeric input to limit number of output correlation pairings
#'
#' @return a table (or dataframe) with the correlation coefficients, p-values, fisher r-to-z statistic and BH p-value correlation for each correlation pair
#'
#' @author Emily Mears, \email{mears.emilyrose@gmail.com}, Matthew Grant, \email{mgra576@aucklanduni.ac.nz}
#' @author Ben Day, \email{benjamindayengineer@gmail.com}
#'
#' @examples
#'## Load example dataframes
#'df1 <- read.csv("example_data/excorr_df1.csv")
#'df2 <- read.csv("example_data/excorr_df2.csv")
#'
#'## Run function
#'multi.data.cor(df1, df2, sample_col = "sample", common_cols = c("sex", "sample"), group = "sex")
#'multi.data.cor(df1, df2, sample_col = "sample", common_cols = c("sex", "sample"), group = "sex", ordered = "Group_2_Pvalue", limit = 100)
#'
#' @export
#'



multi.data.cor <-  function(data1, data2, sample_col, common_cols, group, ordered = "fisher", limit = NA){

  ##--- First check data inputs
  # Check for all-zero variables/cols
  if (TRUE %in% lapply(colnames(data1), function(x) all(data1[, x] == 0)) == TRUE) {
    data1 <- data1[, colSums(data1 != 0) > 0]
    #a <- sapply(colnames(data1), function(x) all(data1[, x] == 0), simplify = FALSE)
    #zvars <- names(a[a == TRUE])
    warning(paste("Some columns contain only zeros and were removed before calculation."))
    #warning(paste(zvars, "removed.", sep = " "))
  }
  if (TRUE %in% lapply(colnames(data2), function(x) all(data2[, x] == 0)) == TRUE) {
    data2 <- data2[, colSums(data2 != 0) > 0]
    #a <- sapply(colnames(data2), function(x) all(data2[, x] == 0), simplify = FALSE)
    #zvars <- names(a[a == TRUE])
    warning(paste("Some columns contain only zeros and were removed before calculation."))
    #cat("Removed variables: ", zvars, "\n", sep="\t")
  }

  # Stop function if sample columns do not match
  if (any(data1[, sample_col] != data2[, sample_col]) == TRUE)
    stop("Sample column must be identical across datasets (and in the same order).")

  # Check if length of datasets is same
  if(max(lengths(data1)) != max(lengths(data2))) {
    stop("Datasets must contain same number of rows")
  }

  # Convert group columns to char
  data1[, which(colnames(data1) == group)] = as.character(data1[, which(colnames(data1) == group)])
  data2[, which(colnames(data2) == group)] = as.character(data2[, which(colnames(data2) == group)])

  # Find both group names
  group_1_name = as.character(data2[,group][1])
  group_2_name = as.character(data2[data2[,group] != as.character(data2[,group][1]),group][1])

  # Check if group column contains more than 2 groups
  if(length(unique(as.character(data1[data1[,group] != as.character(data1[,group][1]),group]))) > 1) {
    stop("More than 2 groups in data1 grouping variable")}
  if(length(unique(as.character(data2[data2[,group] != as.character(data2[,group][1]),group]))) > 1) {
    stop("More than 2 groups in data2 grouping variable")}



  ##--- Define existing functions to use
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
    return(list(diff = dz, pval = pv))
  }



  ##--- Tidy datasets
  # Ensure sample cols is treated as char
  data1[, which(colnames(data1) == sample_col)] <- as.character(data1[, which(colnames(data1) == sample_col)])
  data2[, which(colnames(data2) == sample_col)] <- as.character(data2[, which(colnames(data2) == sample_col)])
  # Remove non-numeric vars from datasets 1 & 2
  nums1 <- base::Filter(is.numeric, data1)
  non_nums1 <- data1[, which(colnames(data1) %in% common_cols)]
  data1 <- cbind(non_nums1, nums1)

  nums2 <- base::Filter(is.numeric, data2)
  non_nums2 <- data2[, which(colnames(data2) %in% common_cols)]
  data2 <- cbind(non_nums2, nums2)



  ### METHOD
  # Run multi.cor.var.compare coding (one variable at a time) to cover full dataset

  # Identify all variables in dataset 1
  data1_nums <- colnames(base::Filter(is.numeric, data1))

  # Initialise output dataframe
  output <- data.frame(matrix(NA, nrow = 0, ncol = 8))


  # Iterate through all variables in data1
  for (j in 1:length(data1_nums)) {

    variable <- data1_nums[j]

    # Begin pairwise correlations
    cond1_cor_pvals = list()
    for(i in 1:ncol(data2[, lapply(data2, class) == "numeric"])){
      cond1_cor_pvals[i] = cor.test(x = data1[data1[,group] == as.character(data1[,which(colnames(data1) == group)])[1], variable],
                                    y = data2[data2[,group] == as.character(data2[,which(colnames(data2) == group)])[1], lapply(data2, class) == "numeric"][,i])$p.value
    }

    cond2_cor_pvals = list()
    for(i in 1:ncol(data2[, lapply(data2, class) == "numeric"])){
      cond2_cor_pvals[i] = cor.test(x = data1[data1[,group] != as.character(data1[,which(colnames(data1) == group)])[1], variable],
                                    y = data2[data2[,group] != as.character(data2[,which(colnames(data2) == group)])[1], lapply(data2, class) == "numeric"][,i])$p.value
    }

    # Vectorise cor.test pvalue outputs
    cond1_cor_pvals = unlist(cond1_cor_pvals)
    cond2_cor_pvals = unlist(cond2_cor_pvals)

    cond1_cor <- cor(x = data2[data2[,group] == as.character(data2[,which(colnames(data2) == group)])[1], lapply(data2, class) == "numeric"],
                     y = data1[data1[,group] == as.character(data1[,which(colnames(data1) == group)])[1], variable])
    cond2_cor <- cor(x = data2[data2[,group] != as.character(data2[,which(colnames(data2) == group)])[1], lapply(data2, class) == "numeric"],
                     y = data1[data1[,group] != as.character(data1[,which(colnames(data1) == group)])[1], variable])

    # Combine variable correlations and pvlaues for each group
    var_cor <- cbind(cond1_cor, cond1_cor_pvals, cond2_cor, cond2_cor_pvals)

    # Then append Fisher r to z as a 3rd column and order
    fisher <- compcorr(n1 = length(data1[data1[,group] == as.character(data1[,which(colnames(data1) == group)])[1], variable]),
                      r1 = var_cor[,1],
                      n2 = length(data1[data1[,group] != as.character(data1[,which(colnames(data1) == group)])[1], variable]),
                      r2 = var_cor[,3])$pval

    # Append fisher transformation and pval adjustments for multiple testing
    p_adjust <- p.adjust(fisher, method = "BH", n = length(fisher))
    var_cor_fish <- as.data.frame(cbind(var_cor, fisher, p_adjust, Var1 = variable, Var2 = rownames(var_cor)))

    # Deal with rownames
    rownames(var_cor_fish) <- 1:nrow(var_cor_fish)

    # Bind to output dataframe
    output <- rbind(output, var_cor_fish)

  }

  # Rename columns
  colnames(output) <- c(paste(group_1_name, "correlation coefficent"),
                        paste(group_1_name, "correlation p-value"),
                        paste(group_2_name, "correlation coefficent"),
                        paste(group_2_name, "correlation p-value"),
                        "Fisher R to Z p-value",
                        "BH p-value adjustment",
                        "Variable from data1",
                        "Variable from data2")

  # Reorder Var1 at front
  output <- output[, c(7, 8, 1, 2, 3, 4, 5, 6)]


  # Round correlation coefficients
  output[,3] <- round(as.numeric(as.character(output[,3])), digits = 2)
  output[,5] <- round(as.numeric(as.character(output[,5])), digits = 2)

  # Significant figures p-values
  output[,4] <- signif(as.numeric(as.character(output[,4])), digits = 4)
  output[,6] <- signif(as.numeric(as.character(output[,6])), digits = 4)
  output[,7] <- signif(as.numeric(as.character(output[,7])), digits = 4)
  output[,8] <- signif(as.numeric(as.character(output[,8])), digits = 4)


  # Ordering options
  if (grepl("fisher", ordered) == TRUE | grepl("Fisher", ordered) == TRUE) {
    output <- output[order(output[,"Fisher R to Z p-value"]),]
  }

  else if (ordered == "g1cor"){
    output <- output[order(output[,3]),]
  }

  else if (ordered == "g1p"){
    output <- output[order(output[,"Group 1 cor p-value"]),]
  }

  else if (ordered == "g2cor"){
    output <- output[order(output[,5]),]
  }

  else if (ordered == "g2p"){
    output = output[order(output[,"Group 2 cor p-value"]),]
  }

  else if (grepl("BH", ordered) == TRUE | grepl("bh", ordered) == TRUE) {
    output = output[order(output[,"BH p-value adjustment"]),]
  }

  else {
    output <- output[order(output[,"Fisher R to Z p-value"]),]
    warning("Undefined ordering column.")
  }

  # Last step is to apply limit
  if(is.na(limit)){
    return(output)
  }

  else if(!is.na(limit)){
    return(output[1:limit, ])
  }

}
