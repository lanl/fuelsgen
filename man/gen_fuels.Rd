% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R
\name{gen_fuels}
\alias{gen_fuels}
\title{Generate fuel layouts}
\usage{
gen_fuels(
  dimX,
  dimY,
  density,
  radius,
  sd_radius,
  height = NULL,
  sd_height = 0,
  heterogeneity,
  heterogeneity.scale = 1,
  sp.cov.locs = NULL,
  sp.cov.vals = NULL,
  sp.cov.scale = NULL,
  reps = 1,
  GP.init.size = 1000,
  seed = NULL,
  logis.scale = 0.217622,
  parallel = F,
  verbose = T
)
}
\arguments{
\item{dimX}{(0,Inf): The size of the domain in X direction (meters)}

\item{dimY}{(0,Inf): The size of the domain in X direction (meters)}

\item{density}{(0,Inf): The mean number of fuels items per unit area (square meter)}

\item{radius}{(0,Inf): The Average radius of a circular fuel element}

\item{sd_radius}{(0,Inf): The standard deviation of fuel element radii}

\item{height}{(0,Inf): The Average height of a fuel element. Default is NULL indicating we should not sample heights.}

\item{sd_height}{(0,Inf): The standard deviation of fuel element heights. Default is 0 indicating no variance in heights.}

\item{heterogeneity}{(0,Inf): How heterogeneous the fuel maps will be. Larger values result in more heterogeneity.}

\item{heterogeneity.scale}{(0,Inf): How big of an effect should the heterogeneity have? Scaling factor.}

\item{reps}{(integer): The number of fuel layouts to generate}

\item{seed}{(integer): Optional random seed for reproducibility}
}
\description{
Generate fuel layouts on a rectangular domain using a handful of parameters to a generative model.
}
