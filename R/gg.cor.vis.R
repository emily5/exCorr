#' Differential correlation visualisation
#'
#' Creates two correlation plots between two variables, separated by a grouping variable containing two groups.
#'
#' @param var_1 character string indicating the first variable for differential correlation
#' @param var_2 character string indicating the second variable for differential correlation
#' @param data dataframe for analysis
#' @param group character string defining the grouping variable for comparative differential correlations
#' @param group_1_title character string defining the first plot title
#' @param group_2_title character string defining the title second plot title
#'
#' @return dual correlation plot for each group with correlation coefficents and p-values
#'
#' @author Emily Mears, \email{mears.emilyrose@gmail.com}, Matthew Grant
#' @author Ben Day, \email{benjamindayengineer@gmail.com}
#'
#' @examples
#'## Load example dataframes
#'df <- read.csv("example_data/excorr_df1.csv")
#'
#'## Run function
#'gg.cor.vis(var_1 = "NTRpFI", var_2 = "X70J8B6", df, group = "sex")
#'gg.cor.vis(var_1 = "NTRpFI", var_2 = "X70J8B6", df, group = "sex", group_1_title = "Group 1 Correlation", group_2_title = "Group 2 Correlation")
#'
#' @export
#'

gg.cor.vis <- function(var_1, var_2, data, group, group_1_title = "Group 1 Correlation", group_2_title = "Group 2 Correlation"){

  #Changes the input factor into a character string, which can be easily manipulated
  data[,which(colnames(data) == group)] = as.character(data[,which(colnames(data) == group)])

  #Extract the names for the two levels of the factor used, so this information can be used down-stream in the legend
  group_1_name = as.character(data[,group][1])
  group_2_name = as.character(data[data[,group] != as.character(data[,group][1]),group][1])

  # Check variable exists in data
  if (!(var_1 %in% colnames(data)))
    stop("Variable not in dataset")
  else if (!(var_2 %in% colnames(data)))
    stop("Variable not in dataset")

  # Check if group column contains more than 2 groups
  if(length(unique(as.character(data[data[,group] != as.character(data[,group][1]),group]))) > 1) {
    stop("More than 2 groups in grouping variable")}

  # function to append stats to PolyQ plots
  ggplottrend <- function (fit) {

    require(ggplot2)

    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm") +
      labs(title = "paste here")
      #geom_text(aes(label = label), vjust = 0.2)                                     #####
      #geom_text(aes(label = signif(fit$fitted.values, digits = 4)), vjust = 0.2)      ######
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

  # FIT MODELS

  model1 <- lm(data[data[,group] == as.character(data[,which(colnames(data) == group)])[1] , var_1] ~ data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], var_2]) # for condition 1
  model2 <- lm(data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], var_1] ~ data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], var_2]) # for condition 2

  p1 <- ggplottrend(model1)+
    labs(title = group_1_title)+
    xlab(var_2)+
    ylab(var_1)+
    labs(subtitle = paste("r =", signif(cor.test(data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], var_2], data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], var_1])$estimate[[1]],2),
                          "    ",
                          "P =", signif(cor.test(data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], var_2], data[data[,group] == as.character(data[,which(colnames(data) == group)])[1], var_1])$p.value[[1]],4))) +      ######
    stat_smooth(method = "lm", col = "blue")


  p2 <- ggplottrend(model2)+
    labs(title = group_2_title)+
    xlab(var_2)+
    ylab(var_1)+
    labs(subtitle = paste("r =", signif(cor.test(data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], var_2], data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], var_1])$estimate[[1]],2),
                          "    ",
                          "P =", signif(cor.test(data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], var_2], data[data[,group] != as.character(data[,which(colnames(data) == group)])[1], var_1])$p.value[[1]],4))) +       #######
    stat_smooth(method = "lm", col = "red")

  # visualise both plots
  multiplot(p1, p2, cols = 2)
}
