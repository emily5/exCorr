% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.cor.R
\name{data.cor}
\alias{data.cor}
\title{Calculates differential correlation statistics for all the possible variable-variable combinations within a single dataset.}
\usage{
data.cor(data, group, ordered = "fisher", limit = NA)
}
\arguments{
\item{data}{dataframe for analysis}

\item{group}{character string defining the grouping variable for comparative differential correlations}

\item{ordered}{character string defining the variable to order output by. Choose from \code{g1cor}, \code{g1p}, \code{g2cor}, \code{g2r}, \code{fisher} (default) and \code{BH}.}

\item{limit}{numeric input to limit number of output correlation pairings}
}
\value{
a table (or dataframe) with correlation coefficients, p-values, fisher r-to-z p-value and adjusted p-value using Benjamini-Hochberg corretion for each correlation pair, with separation of two groups for comparison.
}
\description{
Creates a table of pairwise correlations statistics for all variables within a dataframe, with separation of two groups for comparison. This can be used as an exploratory tool to investigate differential correlations within a dataset.
Requires grouping variable to compare two groups.
}
\examples{
## Load example dataframes
df <- read.csv("example_data/excorr_df2.csv")

## Run function
data.cor(df, group = "sex")
data.cor(df, group = "sex", ordered = "fisher", limit = 20)

}
\author{
Emily Mears, \email{mears.emilyrose@gmail.com}, Matthew Grant, \email{mgra576@aucklanduni.ac.nz}

Ben Day, \email{benjamindayengineer@gmail.com}
}
