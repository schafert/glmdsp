% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/threshold_nb.R
\name{abco_nb}
\alias{abco_nb}
\title{Title}
\usage{
abco_nb(
  y,
  D = 1,
  useAnom = FALSE,
  useObsSV = FALSE,
  nsave = 1000,
  nburn = 1000,
  nskip = 4,
  mcmc_params = list("mu", "yhat", "evol_sigma_t2", "r", "dhs_phi", "dhs_phi2",
    "dhs_mean", "gamma", "omega"),
  r_init = NULL,
  r_sample = "int_mh",
  step = 1,
  mu_init = NULL,
  mu_sample = TRUE,
  prior_r = rlang::expr(log(1 + x^2/100)),
  evol0_sample = TRUE,
  evol0_sd = 10,
  thresh_sample = TRUE,
  thresh_init = NULL,
  thresh_proposal = "unif",
  thresh_proposal_sd = 0.1,
  omega_trans = rlang::expr(log(omega^2 + 1e-04)),
  chol0 = NULL,
  computeDIC = TRUE,
  offset = 0,
  verbose = TRUE,
  cp_thres = 0.5,
  seed = NULL
)
}
\arguments{
\item{seed}{}
}
\description{
Title
}
