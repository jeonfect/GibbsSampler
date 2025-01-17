% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GibbsSampler.R
\name{gibbs}
\alias{gibbs}
\title{Gibbs Sampler for model-based clustering}
\usage{
gibbs(x, k, m_iter, init, V, burn_in, theta_out = FALSE, z_out = FALSE)
}
\arguments{
\item{x}{matrix \verb{n x p} of observations}

\item{k}{number of clusters}

\item{m_iter}{integer number of iterations used for Gibbs sampler}

\item{init}{list of initial values of parameters(\code{mu},\code{Sigma},\code{theta},\code{rho},\code{z}). If not specified, this function assigns initial values based on \code{x}}

\item{V}{martrix \verb{p x p x N} of variances corresponding to each observation}

\item{burn_in}{positive integer number of iterations that are only used for throwing out garbage values. If not specified, 10\% of \code{m_iter} is used as default}

\item{theta_out}{logical value whether to save every \code{theta} produced through all iterations.. Default is \code{FALSE}}

\item{z_out}{logical value whether to save every \code{z} updated through all iterations. Default is \code{FALSE}.}
}
\value{
a list with elements:
\describe{
\item{\code{mu}}{matrix \verb{n x p x m_iter} of mean values for clutsers}
\item{\code{sigma}}{matrix \verb{p x p x k x m_iter} of within-group variance}
\item{\code{rho}}{matrix \verb{n x k x m_iter} of proportion of group membership}
\item{\code{lambda}}{matrix \verb{n x k} of average probability of group membership}
\item{\code{z}}{vector length \code{n} of posterior group membership if \code{z_out=FALSE}. matrix \verb{n x m_iter} is produced otherwise}
\item{\code{theta}}{matrix \verb{n x p} of observation means if \code{theta_out=FALSE}, matrix \verb{n x p x m_iter} is produced otherwise}
\item{\code{burn}}{integer number of burn-in periods}
\item{\code{m_iter}}{integer number of total iterations}
}
}
\description{
Gibbs Sampler for model-based clustering
}
