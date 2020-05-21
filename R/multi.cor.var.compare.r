#' Calculates differential correlation statistics between a specified variable and all other variables within a different dataset.
#'
#' Creates a table of pairwise correlations statistics, with separation of two groups for comparison.
#' This can be used as an exploratory tool to investigate correlations between a specific variable of interest and all other variables within a separate dataset.
#'
#' @param variable character string indicating the variable of interest for differential correlation (within \code{data1})
#' @param data1 first dataframe containing variable of interest
#' @param data2 second dataframe containing variable to compare
#' @param sample_col character string defining the sample or identifying column in both datasets (common to both dataframes)
#' @param group character string defining the grouping variable for comparative differential correlations
#' @param ordered character string defining which column the table should be ordered by. Choose from \code{g1cor}, \code{g1p}, \code{g2cor}, \code{g2r}, \code{fisher} (default) and \code{BH}.
#'
#' @return a table (or dataframe) with Pearson correlation coefficients (r), associated p-values, fisher r-to-z statistic and BH p-value correlation for each correlation pair
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
#'multi.cor.var.compare(variable = "NTRpFI", df1, df2, sample_col = "sample", group = "sex")
#'multi.cor.var.compare(variable = "NTRpFI", df1, df2, sample_col = "sample", group = "sex", ordered = "fisher")
#'
#' @export
#'

multi.cor.var.compare <- function(variable, data1, data2, sample_col, group, ordered = "fisher"){

  ##--- First check data inputs

  # Convert group columns to char
  data1[, which(colnames(data1) == group)] = as.character(data1[, which(colnames(data1) == group)])
  data2[, which(colnames(data2) == group)] = as.character(data2[, which(colnames(data2) == group)])

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

  # Check variable exists in data
  if (!(variable %in% colnames(data1)))
    stop("Variable not in dataset 1")

  # Check if length of datasets is same
  if(max(lengths(data1)) != max(lengths(data2))) {
    stop("Datasets must contain same number of rows")
  }

  # Check if group column contains more than 2 groups
  if(length(unique(as.character(data1[data1[,group] != as.character(data1[,group][1]),group]))) > 1) {
    stop("More than 2 groups in data1 grouping variable")}
  if(length(unique(as.character(data2[data2[,group] != as.character(data2[,group][1]),group]))) > 1) {
    stop("More than 2 groups in data2 grouping variable")}






  ##----- Define functions to be used

  # Fisher's Z-transformation
  compcor <- function(n1, r1, n2, r2){

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





  # Find both group names
  group_1_name = as.character(data2[,group][1])
  group_2_name = as.character(data2[data2[,group] != as.character(data2[,group][1]),group][1])

  # Check if group column contains more than 2 groups
  if(length(unique(as.character(data2[data2[,group] != as.character(data2[,group][1]),group]))) > 1)
    stop("More than 2 groups in grouping variable")

  #-----------------

  # Begin pairwise correlations
  cond1_cor_pvals = list()
  for(i in 1:ncol(data2[, lapply(data2, class) == "numeric"])){
    cond1_cor_pvals[i] = cor.test(method = "pearson", x = data1[data1[,group] == as.character(data1[,which(colnames(data1) == group)])[1], variable],
                                  y = data2[data2[,group] == as.character(data2[,which(colnames(data2) == group)])[1], lapply(data2, class) == "numeric"][,i])$p.value
  }

  cond2_cor_pvals = list()
  for(i in 1:ncol(data2[, lapply(data2, class) == "numeric"])){
    cond2_cor_pvals[i] = cor.test(method = "pearson", x = data1[data1[,group] != as.character(data1[,which(colnames(data1) == group)])[1], variable],
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
  fisher <- compcor(n1 = length(data1[data1[,group] == as.character(data1[,which(colnames(data1) == group)])[1], variable]),
                    r1 = var_cor[,1],
                    n2 = length(data1[data1[,group] != as.character(data1[,which(colnames(data1) == group)])[1], variable]),
                    r2 = var_cor[,3])$pval

  # Append fisher transformation and pval adjustments for multiple testing
  p_adjust <- p.adjust(fisher, method = "BH", n = length(fisher))
  var_cor_fish <- cbind(var_cor, fisher, p_adjust)

  # Rename columns
  colnames(var_cor_fish) <- c(paste(variable, "~ Correlations in", group_1_name),
                              "Correlation p-value",
                              paste(variable, "~ Correlations in", group_2_name),
                              "Correlation p-value",
                              "Fisher R to Z",
                              "BH p-value adjustment")

  # Round correlation coefficients
  var_cor_fish[,1] <- round(var_cor_fish[,1], digits = 2)
  var_cor_fish[,3] <- round(var_cor_fish[,3], digits = 2)

  # Significant figures p-values
  var_cor_fish[,2] <- signif(var_cor_fish[,2], digits = 4)
  var_cor_fish[,4] <- signif(var_cor_fish[,4], digits = 4)
  var_cor_fish[,5] <- signif(var_cor_fish[,5], digits = 4)
  var_cor_fish[,6] <- signif(var_cor_fish[,6], digits = 4)


  # Ordering options
  if(ordered == "fisher"){
    var_cor_fish <- var_cor_fish[order(var_cor_fish[,5]),]
  }

  else if(ordered == "g1cor"){
    var_cor_fish <- var_cor_fish[order(var_cor_fish[,1]),]
  }

  else if(ordered == "g1p"){
    var_cor_fish <- var_cor_fish[order(var_cor_fish[,2]),]
  }

  else if(ordered == "g2cor"){
    var_cor_fish <- var_cor_fish[order(var_cor_fish[,3]),]
  }

  else if(ordered == "g2p"){
    var_cor_fish = var_cor_fish[order(var_cor_fish[,4]),]
  }

  else if(ordered == "BH"){
    var_cor_fish = var_cor_fish[order(var_cor_fish[,6]),]
  }

  else {
    var_cor_fish <- var_cor_fish[order(var_cor_fish[,5]),]
    warning("Undefined ordering column.")
  }

  # Return output
  return(as.data.frame(var_cor_fish))

}
