% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/NanoStringQualityMetrics.R
\name{flagSamplesCount}
\alias{flagSamplesCount}
\title{flagSamplesCount}
\usage{
flagSamplesCount(rccSet, method = c("cutoffByMMAD", "cutoffByVar"), cutoff,
  reference = c("blank", "negatives"), stringency = 2.5, maxMiss)
}
\arguments{
\item{rccSet}{NanoString ExpressionSet object}

\item{method}{scalar character, specifying the method for outlier detection; either cutoffByMMAD or cutoffByVar}

\item{cutoff}{scalar numeric, cutoff in method to determine outliers}

\item{reference}{scalar character, reference to determine limit of detection, passed on to \code{\link{lodAssess}}}

\item{stringency}{passed on to \code{\link{lodAssess}}}

\item{maxMiss}{scalar numeric, indicating the fraction of genes accepted to be missed}
}
\value{
A numeric vector giving the indices of samples with outlier values based on overall counts
}
\description{
Flag samples based on overall counts
}
\details{
Outliers are determined based on method and cutoff. Samples defined as outliers in the following categories
will be flagged: the sum of counts of endogeneous genes, the ratio of sums for positive controls to the sums of
endogenous genes, the number of genes below the lower limit of detection
}
\author{
Dorothee Nickles
}
