% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi.gg.cor.vis.R
\name{multi.gg.cor.vis}
\alias{multi.gg.cor.vis}
\title{Differential correlation visualisation}
\usage{
multi.gg.cor.vis(var_1, var_2, data1, data2, sample_col, group,
  group_1_title = "Group 1 Correlation",
  group_2_title = "Group 2 Correlation")
}
\arguments{
\item{var_1}{character string defining the first variable for differential correlation}

\item{var_2}{character string defining the second variable for differential correlation}

\item{data1}{dataframe containing var_1}

\item{data2}{dataframe containing var_2}

\item{sample_col}{character string defining the sample or identifying column in both dataframes}

\item{group}{character string defining the grouping variable for comparative differential correlations}

\item{group_1_title}{character string defining the first plot title}

\item{group_2_title}{character string defining the second plot title}
}
\value{
two correlation plots for each group with correlation coefficents and p-values
}
\description{
Creates two correlation plots between two variables, separated by a grouping variable. Requires datasets with same samples.
}
\examples{
## Load example dataframes
df1 <- read.csv("example_data/excorr_df1.csv")
df2 <- read.csv("example_data/excorr_df2.csv")

## Run function
multi.gg.cor.vis(var_1 = "NTRpFI", var_2 = "X70J8B6", df1, df2, group = "sex")
multi.gg.cor.vis(var_1 = "NTRpFI", var_2 = "X70J8B6", df1, df2, group = "sex", group_1_title = "Group 1 Correlation", group_2_title = "Group 2 Correlation")


}
\author{
Emily Mears, \email{mears.emilyrose@gmail.com}, Matthew Grant, \email{mgra576@aucklanduni.ac.nz}

Ben Day, \email{benjamindayengineer@gmail.com}
}
