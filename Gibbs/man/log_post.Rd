% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/condPDF.R
\name{log_post}
\alias{log_post}
\title{Unnormalized complete data log-posterior}
\usage{
log_post(mu, Sigma, rho, y, V, theta, z, Omega, v_k)
}
\arguments{
\item{mu}{matrix \verb{k x p} of group mean}

\item{Sigma}{matrix \verb{p x p x k} of within-group variance}

\item{rho}{vector length \code{k} of prior group membership probability}

\item{y}{matrix \verb{n x p} of observation}

\item{V}{matrix \verb{p x p x n} of observation variance}

\item{theta}{matrix \verb{n x p} of observation mean}

\item{z}{vector length \code{n} of group membership}

\item{Omega}{matrix \verb{p x p x k} of between row precision}

\item{v_k}{vector length \code{k} of prior contribution to each cluster}
}
\value{
integer value of log posterior for given parameters
}
\description{
Unnormalized complete data log-posterior
}
