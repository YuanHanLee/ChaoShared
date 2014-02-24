\name{ChaoShared.Basic}
\alias{ChaoShared.Basic}
\title{Basic data information for shiny}
\description{
Basic data information for shiny app with a list form by data.frame
}
\usage{
ChaoShared.Basic(data, datatype = c("abundance", "incidence"), nboot = 200)
}
\arguments{
  \item{data}{a numerical matrix or a data frame with two columns. Each column represent the species abundances or incidence frequencies of each community. If datatype = "incidence", then the input format of first entry should be total number of sampling units, and followed by species incidence frequencies in each column (See examples).}
  \item{datatype}{the data type of input data. That is individual-based abundance data (datatype = "abundance") or presence/absence sample-base incidence data (datatype = "incidence").}
  \item{nboot}{the number of bootstrap resampling times, default is 200.}
}