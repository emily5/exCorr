#' Differential correlation visualisation
#'
#' Creates two correlation plots between two variables, separated by a grouping variable
#'
#' @var_1 the first variable for corre;ation analyses
#'
#' @var_2 the second variable for correlation analyses
#'
#' @data the dataset to be used
#'
#' @group the variable to distinguish between correlation groups
#'
#' @author Emily Mears, \email{mears.emilyrose@gmail.com}
#'
#' @export
#'

gg.corr.vis <- function(var_1, var_2, data, group){
  
    #Changes the input factor into a character string, which can be easily manipulated
  data[,which(colnames(data) == group)] = as.character(data[,which(colnames(data) == group)])
  
  #Extract the names for the two levels of the factor used, so this information can be used down-stream in the legend
  group_1_name = as.character(data[,group][1])
  group_2_name = as.character(data[data[,group] != as.character(data[,group][1]),group][1])
  
  #Change the character string occuring first into 1 and the character strings that do not equal the first into 2.
  data[which(as.character(data[,which(colnames(data) == group)]) == as.character(data[,which(colnames(data) == group)])[1]), which(colnames(data) == group)] = as.numeric(1)
  data[which(as.character(data[,which(colnames(data) == group)]) != as.character(data[,which(colnames(data) == group)])[1]), which(colnames(data) == group)] = as.numeric(2)
  
  # function to append stats to PolyQ plots
  ggplottrend <- function (fit) {
    
    require(ggplot2)
    
    ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
      geom_point() +
      stat_smooth(method = "lm")+
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

  # FIT MODELS
  
  model1 <- lm(data[data[,group] == 1 , var_1] ~ data[data[,group] == 1, var_2]) # for condition 1
  model2 <- lm(data[data[,group] == 2, var_1] ~ data[data[,group] == 2, var_2]) # for condition 2
  
  p1 <- ggplottrend(model1)+
    labs(title = "Control correlation")+
    xlab(var_2)+
    ylab(var_1)+
    labs(subtitle = paste("r =", signif(cor.test(data[data[,group] == 1, var_2], data[data[,group] == 1, var_1])$estimate[[1]],2), 
                          "    ",
                          "P =", signif(cor.test(data[data[,group] == 1, var_2], data[data[,group] == 1, var_1])$p.value[[1]],4)))+
    stat_smooth(method = "lm", col = "blue")
  
  
  p2 <- ggplottrend(model2)+
    labs(title = "OVT73 correlation")+
    xlab(var_2)+
    ylab(var_1)+
    labs(subtitle = paste("r =", signif(cor.test(data[data[,group] == 2, var_2], data[data[,group] == 2, var_1])$estimate[[1]],2),
                          "    ",
                          "P =", signif(cor.test(data[data[,group] == 2, var_2], data[data[,group] == 2, var_1])$p.value[[1]],4)))+
    stat_smooth(method = "lm", col = "red")
  
  # visualise both plots
  multiplot(p1,p2,cols = 2)
}

