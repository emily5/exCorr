#' Differential correlation visualisation
#'
#' Creates two correlation plots between two variables, separated by a grouping variable. Requires datasets with same samples.
#'
#' @param var_1 character string defining the first variable for differential correlation
#' @param var_2 character string defining the second variable for differential correlation
#' @param data1 dataframe containing var_1
#' @param data2 dataframe containing var_2
#' @param sample_col character string defining the sample or identifying column in both dataframes
#' @param group character string defining the grouping variable for comparative differential correlations
#' @param group_1_title character string defining the first plot title
#' @param group_2_title character string defining the second plot title
#'
#' @return two correlation plots for each group with correlation coefficents and p-values
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
#'multi.gg.cor.vis(var_1 = "NTRpFI", var_2 = "X70J8B6", df1, df2, group = "sex")
#'multi.gg.cor.vis(var_1 = "NTRpFI", var_2 = "X70J8B6", df1, df2, group = "sex", group_1_title = "Group 1 Correlation", group_2_title = "Group 2 Correlation")
#'
#'
#' @export
#'

multi.gg.cor.vis <- function(var_1, var_2, data1, data2, sample_col, group,
                              group_1_title = "Group 1 Correlation", group_2_title = "Group 2 Correlation") {

  ##--- First check data inputs
  # Check if group column contains more than 2 groups
  if(length(unique(as.character(data1[data1[,group] != as.character(data1[,group][1]),group]))) > 1) {
    stop("More than 2 groups in data1 grouping variable")}
  if(length(unique(as.character(data2[data2[,group] != as.character(data2[,group][1]),group]))) > 1) {
    stop("More than 2 groups in data2 grouping variable")}

  # Stop function if sample columns do not match
  if (any(data1[, sample_col] != data2[, sample_col]) == TRUE)
    stop("Sample column must be identical across datasets (and in the same order).")

  # Check if length of datasets is same
  if(max(lengths(data1)) != max(lengths(data2))) {
    stop("Datasets must contain same number of rows")
  }

  # Check variable exists in data
  if (!(var_1 %in% colnames(data1)))
    stop("Var_1 not in dataset 1")
  else if (!(var_2 %in% colnames(data2)))
    stop("Var_2 not in dataset 2")

  # Convert group columns to char
  data1[, which(colnames(data1) == group)] = as.character(data1[, which(colnames(data1) == group)])
  data2[, which(colnames(data2) == group)] = as.character(data2[, which(colnames(data2) == group)])

  #Extract the names for the two levels of the factor used, so this information can be used down-stream in the legend
  group_1_name = as.character(data2[,group][1])
  group_2_name = as.character(data2[data2[,group] != as.character(data2[,group][1]),group][1])


  ##--- Define functions
  # function to append stats to PolyQ plots
  ggplottrend <- function (fit) {

    require(ggplot2)

    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm") +
      labs(title = "paste here")

  }


  # function to visualise multiple plots
  multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    require(grid)

    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    if (is.null(layout)) {
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots == 1) {
      print(plots[[1]])

    } else {
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      for (i in 1:numPlots) {
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }

  }



  ##--- Fit Models

  model1 <- lm(data1[data1[,group] == as.character(data1[,which(colnames(data1) == group)])[1], var_1] ~ data2[data2[,group] == as.character(data2[,which(colnames(data2) == group)])[1], var_2]) # for condition 1
  model2 <- lm(data1[data1[,group] != as.character(data1[,which(colnames(data1) == group)])[1], var_1] ~ data2[data2[,group] != as.character(data2[,which(colnames(data2) == group)])[1], var_2]) # for condition 2

  p1 <- ggplottrend(model1) +
    labs(title = group_1_title) +
    xlab(var_2) +
    ylab(var_1) +
    labs(subtitle = paste("r =", signif(cor.test(method = "pearson", x = data2[data2[,group] == as.character(data2[,which(colnames(data2) == group)])[1], var_2],
                                                 y = data1[data1[,group] == as.character(data1[,which(colnames(data1) == group)])[1], var_1])$estimate[[1]],2),
                          "    ",
                          "P =", signif(cor.test(method = "pearson", x = data2[data2[,group] == as.character(data2[,which(colnames(data2) == group)])[1], var_2],
                                                 y = data1[data1[,group] == as.character(data1[,which(colnames(data1) == group)])[1], var_1])$p.value[[1]],4))) +      ######
    stat_smooth(method = "lm", col = "blue")

  p2 <- ggplottrend(model2) +
    labs(title = group_2_title) +
    xlab(var_2) +
    ylab(var_1) +
    labs(subtitle = paste("r =", signif(cor.test(method = "pearson", x = data2[data2[,group] != as.character(data2[,which(colnames(data2) == group)])[1], var_2],
                                                 y = data1[data1[,group] != as.character(data1[,which(colnames(data1) == group)])[1], var_1])$estimate[[1]],2),
                          "    ",
                          "P =", signif(cor.test(method = "pearson", x = data2[data2[,group] != as.character(data2[,which(colnames(data2) == group)])[1], var_2],
                                                 y = data1[data1[,group] != as.character(data1[,which(colnames(data1) == group)])[1], var_1])$p.value[[1]],4))) +      ######
    stat_smooth(method = "lm", col = "red")

  # visualise both plots
  multiplot(p1, p2, cols = 2)
}
