#' Calculates differential correlation statistics between a specified variable and all other variables within a single dataset.
#'
#' Creates a table of pairwise correlations statistics, with separation of two groups for comparison.
#' This can be used as an exploratory tool to investigate correlations between a specific variable of interest and all other variables within the same dataset.
#' The significance of the difference between two correlation coefficients is ......
#'
#' @param variable character string defining the variable of interest for differential correlation
#' @param data dataframe for analysis
#' @param group character string defining the grouping variable for comparative differential correlations
#' @param ordered character string defining which column the table should be ordered by. Choose from \code{g1cor}, \code{g1p}, \code{g2cor}, \code{g2r}, \code{fisher} (default) and \code{BH}.
#'
#'
#' @return a table (or dataframe) with pearson correlation coefficients (r) and associated p-values between variables for each comparative group. Also gives fisher r-to-z p-value and adjusted p-value using Benjamini-Hochberg correction for each correlation pair.
#'
#' @author Emily Mears, \email{mears.emilyrose@gmail.com}, Matthew Grant, \email{mgra576@aucklanduni.ac.nz}
#' @author Ben Day, \email{benjamindayengineer@gmail.com}
#'
#' @examples
#'## Load example dataframes
#'df <- read.csv("example_data/excorr_df1.csv")
#'
#'## Run function
#'cor.var.compare(variable = "NTRpFI", data = df, group = "sex")
#'cor.var.compare(variable = "NTRpFI", data = df, group = "sex", ordered = "fisher")
#'
#' @export
#'

cor.var.compare <- function(variable, data, group, ordered = "fisher"){

  # Define existing function to be used
  compcor <- function(n1, r1, n2, r2){
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
    warning(paste("Some columns contain only zeros and were removed."))
  }

  # Convert group columns to char
  data[,which(colnames(data) == group)] = as.character(data[,which(colnames(data) == group)])

  # Find both group names
  group_1_name = as.character(data[,group][1])
  group_2_name = as.character(data[data[,group] != as.character(data[,group][1]),group][1])

  # Check variable exists in data
  if (!(variable %in% colnames(data)))
    stop("Variable not in dataset")

  # Check if group column contains more than 2 groups
  if(length(unique(as.character(data[data[,group] != as.character(data[,group][1]),group]))) > 1)
    stop("More than 2 groups in grouping variable")


  # Begin pairwise correlations
  cond1_cor_pvals = list()
  for(i in 1:ncol(data[, lapply(data, class) == "numeric"])){
    cond1_cor_pvals[i] = cor.test(method = "pearson", x = data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], variable],
                                  y = data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], lapply(data, class) == "numeric"][,i])$p.value
  }

  cond2_cor_pvals = list()
  for(i in 1:ncol(data[, lapply(data, class) == "numeric"])){
    cond2_cor_pvals[i] = cor.test(method = "pearson",x = data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], variable],
                                  y = data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], lapply(data, class) == "numeric"][,i])$p.value
  }

  # Vectorise cor.test pvalue outputs
  cond1_cor_pvals = unlist(cond1_cor_pvals)
  cond2_cor_pvals = unlist(cond2_cor_pvals)

  cond1_cor <- cor(x = data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], lapply(data, class) == "numeric"],
                   y = data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], variable])
  cond2_cor <- cor(x = data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], lapply(data, class) == "numeric"],
                   y = data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], variable])

  # Combine variable correlations and pvlaues for each group
  var_cor <- cbind(cond1_cor, cond1_cor_pvals, cond2_cor, cond2_cor_pvals)


  # Then append Fisher r to z as a 3rd column and order
  fisher <- compcor(n1 = length(data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], variable]),
                    r1 = var_cor[,1],
                    n2 = length(data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], variable]),
                    r2 = var_cor[,3])$pval

  # Append fisher transformation and pval adjustments for multiple testing
  p_adjust <- p.adjust(fisher, method = "BH", n = length(fisher))
  var_cor_fish <- cbind(var_cor, fisher, p_adjust)

  # Rename columns
  colnames(var_cor_fish) <- c(paste(variable, "~ correlation coefficient in", group_1_name),
                              paste(group_1_name, "correlation p-value"),
                              paste(variable, "~ correlation coefficient in", group_2_name),
                              paste(group_2_name, "correlation p-value"),
                              "Fisher R to Z p-value",
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


