#' Differential correlations between a single variable and all variables within a dataset
#'
#' @variable the variable of interest 
#'
#' @data the dataset to be used
#'
#' @group the variable to distinguish between correlation groups
#'
#' @author Emily Mears, \email{mears.emilyrose@gmail.com}
#'
#' @export
#'

cor.var.compare <- function(variable, data, group){
  
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
  
  data[,which(colnames(data) == group)] = as.character(data[,which(colnames(data) == group)])
  
  group_1_name = as.character(data[,group][1])
  group_2_name = as.character(data[data[,group] != as.character(data[,group][1]),group][1])
  
  data[which(as.character(data[,which(colnames(data) == group)]) == as.character(data[,which(colnames(data) == group)])[1]), which(colnames(data) == group)] = 1
  data[which(as.character(data[,which(colnames(data) == group)]) != as.character(data[,which(colnames(data) == group)])[1]), which(colnames(data) == group)] = 2
  
  data[,group] = as.numeric(data[,group])
  
  cond1_cor_pvals = list()
  for(i in 1:ncol(data[, lapply(data, class) == "numeric"])){
    cond1_cor_pvals[i] = cor.test(data[data[,group] == 1, variable], data[data[,group] == 1, lapply(data, class) == "numeric"][,i])$p.value
  }
  
  cond2_cor_pvals = list()
  for(i in 1:ncol(data[, lapply(data, class) == "numeric"])){
    cond2_cor_pvals[i] = cor.test(data[data[,group] == 2, variable], data[data[,group] == 2, lapply(data, class) == "numeric"][,i])$p.value
  }
  cond1_cor_pvals = unlist(cond1_cor_pvals)
  cond2_cor_pvals = unlist(cond2_cor_pvals)
  
  cond1_cor <- cor(data[data[,group] == 1, lapply(data, class) == "numeric"],data[data[,group] == 1, variable])
  cond2_cor <- cor(data[data[,group] == 2, lapply(data, class) == "numeric"],data[data[,group] == 2, variable])
  var_cor <- cbind(cond1_cor, cond1_cor_pvals, cond2_cor, cond2_cor_pvals)
  # then append Fisher r to z as a 3rd column and order
  fisher <- compcor(n1 = length(data[data[,group] == 1,variable]), r1 = var_cor[,1], 
                    n2 = length(data[data[,group] == 2,variable]), r2 = var_cor[,3])$pval
  var_cor_fish <- cbind(var_cor, fisher)
  # then append a column that adjusts p-values for multiple testing
  p_adjust <- p.adjust(fisher, method = "BH", n = length(fisher))
  var_cor_fish <- cbind(var_cor_fish, p_adjust)
  var_cor_fish <- var_cor_fish[order(var_cor_fish[,5]),]
  colnames(var_cor_fish) <- c(paste(variable, "~ Correlations in", group_1_name),
                              "Correlation p-value",
                              paste(variable, "~ Correlations in", group_2_name),
                              "Correlation p-value",
                              "Fisher R to Z",
                              "BH p-value adjustment")
  return(as.data.frame(var_cor_fish))
}
