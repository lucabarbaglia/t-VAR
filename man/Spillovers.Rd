% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spillovers.R
\name{Spillovers}
\alias{Spillovers}
\title{Spillovers}
\usage{
Spillovers(fit, dec.type = "Generalized", h = 1)
}
\arguments{
\item{fit}{t-VAR fit obtained from the function Large.tVAR}

\item{dec.type}{decomposition for forecast error varaince decomposition: "Cholesky", "Spectral", "Generalized", "GeneralizedIRF" or "Generalized_Lanne". Default "Generalized".}

\item{h}{forecast horizon for variance decomposition. Default 1.}
}
\value{
A list containing:
\item{"spill"}{matrix of volatility spillovers.}
\item{"spill_index"}{volatility spillover index.}
}
\description{
Function to compute volatility spillovers and 
its index based on forecast error variance decomposition
}
