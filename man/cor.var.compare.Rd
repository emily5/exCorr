% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor.var.compare.R
\name{cor.var.compare}
\alias{cor.var.compare}
\title{Calculates differential correlation statistics between a specified variable and all other variables within a single dataset.}
\usage{
cor.var.compare(variable, data, group, ordered = "fisher")
}
\arguments{
\item{variable}{character string defining the variable of interest for differential correlation}

\item{data}{dataframe for analysis}

\item{group}{character string defining the grouping variable for comparative differential correlations}

\item{ordered}{character string defining which column the table should be ordered by. Choose from \code{g1cor}, \code{g1p}, \code{g2cor}, \code{g2r}, \code{fisher} (default) and \code{BH}.}
}
\value{
a table (or dataframe) with pearson correlation coefficients (r) and associated p-values between variables for each comparative group. Also gives fisher r-to-z p-value and adjusted p-value using Benjamini-Hochberg correction for each correlation pair.
}
\description{
Creates a table of pairwise correlations statistics, with separation of two groups for comparison.
This can be used as an exploratory tool to investigate correlations between a specific variable of interest and all other variables within the same dataset.
The significance of the difference between two correlation coefficients is ......
}
\examples{
## Load example dataframes
df <- read.csv("example_data/excorr_df1.csv")

## Run function
cor.var.compare(variable = "NTRpFI", data = df, group = "sex")
cor.var.compare(variable = "NTRpFI", data = df, group = "sex", ordered = "fisher")

}
\author{
Emily Mears, \email{mears.emilyrose@gmail.com}, Matthew Grant, \email{mgra576@aucklanduni.ac.nz}

Ben Day, \email{benjamindayengineer@gmail.com}
}
