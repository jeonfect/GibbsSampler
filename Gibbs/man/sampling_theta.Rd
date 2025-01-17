% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GibbsSampler.R
\name{sampling_theta}
\alias{sampling_theta}
\title{Sampling method for observation means.}
\usage{
sampling_theta(x, mu, Sigma, V, z)
}
\arguments{
\item{x}{observation \verb{n x p}}

\item{mu}{matrix of cluster means \verb{k x p}}

\item{Sigma}{matrix of cluster variances \verb{p x p x n}}

\item{V}{matrix of variance for each observation x \verb{p x p x n}}

\item{z}{vector of each observation's cluster allocation, length \code{n}}
}
\value{
matrix \verb{n x p} of parameter values for observation x
}
\description{
Sampling method for observation means.
}
\details{
\code{theta_i} has conditional distribution \verb{N(G_i(y_i-mu_\{z_i\})+mu_\{z_i\},G_iV_i)} where \verb{G_i=Sig_\{z_i\}(V_i+Sig_\{z_i\})^\{-1\}}. All sampling can be done through \code{rRxNorm()} from \code{mniw} package.
}
\references{
Martin Lysy and Bryan Yates (2019). mniw: The Matrix-Normal Inverse-Wishart Distribution. R package version 1.0. https://CRAN.R-project.org/package=mniw
}
