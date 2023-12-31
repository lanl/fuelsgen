% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contour_pairs.R
\name{contour_pairs}
\alias{contour_pairs}
\title{Plot MCMC posteriors}
\usage{
contour_pairs(
  theta.list,
  truth = NULL,
  cols = c("cornflowerblue", "darkorange", "forestgreen"),
  rho.min = 0,
  rho.max = 10,
  rho.prior = "unif",
  mu.min = 0,
  mu.max = 3,
  mu.prior.mean = 1.5,
  mu.prior.sd = 0.5,
  var.min = 0,
  var.max = 0.05,
  var.prior.type = "gamma",
  var.prior.params = c(1, 0.001),
  lambda.max = 0.5,
  lambda.prior.type = "gamma",
  lambda.prior.params = c(1, 0.25),
  labels = c("n=10", "n=25")
)
}
\arguments{
\item{theta.list}{list of matrices containing posterior samples}

\item{truth}{optional vector of known parameter values}

\item{cols}{colors for posterior density estimate plots}

\item{labels}{legend names for posterior density estimates}
}
\description{
Pairs plot of posterior density estimates and prior distributions
}
