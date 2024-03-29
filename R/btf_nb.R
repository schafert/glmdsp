#' MCMC Sampler for Bayesian Trend Filtering with Negative Binomial Observation Model
#'
#' @description
#' Run the MCMC for Bayesian trend filtering with a penalty on first (D = 1) or
#' second (D = 2) differences of the log conditional expectation
#' of a Negative Binomial distribution. The penalty is determined by the prior
#' on the evolution errors, currently only the following options are available:
#' \itemize{
#' \item the dynamic horseshoe prior ('DHS');
#' \item the static horseshoe prior ('HS');
#' }
#' In each case, the evolution error is a scale mixture of Gaussians.
#' Sampling is accomplished with a (parameter-expanded) Gibbs sampler,
#' mostly relying on a dynamic linear model representation.
#'
#' @inheritParams dspCP::btf
#' @param y the length \code{T} vector of time series observations; expected to be in the support of a Negative Binomial distribution
#' @param mcmc_params character list of parameters for which we store the MCMC output;
#' must be one or more of:
#' \itemize{
#' \item 'mu' (conditional mean)
#' \item 'yhat' (posterior predictive distribution)
#' \item 'evol_sigma_t2' (evolution error variance)
#' \item 'r' (overdispersion)
#' \item 'dhs_phi' (DHS AR(1) coefficient)
#' \item 'dhs_mean' (DHS AR(1) unconditional mean)
#' }
#' Defaults to everything.
#' @param r_init numeric; initial value (defaults to 5) of MCMC sampling for overdispersion parameter;
#' must be an integer if r_sample is 'int_mh' or `NULL`
#' @param r_sample character string specifying sampling strategy for overdispersion;
#' must be one of 'int_mh' (integer random walk Metropolis-Hastings; default),
#' 'slice' (slice sampling), 'mh' (uniform random walk Metropolis-Hastings), or
#' `NULL` (not sampled, fixed at `r_init`)
#' @param step numeric; step length of proposal distribution for 'int_mh' and 'mh'
#' values of `r_sample` (defaults to 1)
#' @param mu_init length \code{T} vector of initial values for the log conditional expectation or
#' `NULL` (default) for data driven starting value (`log(y + 1) + noise`)
#' @param mu_sample logical; If `TRUE` (default), the log conditional mean is sampled,
#' else it is fixed at `mu_init`
#' @param prior_r expression of `x` for the proportional log prior density of the overdispersion;
#' defaults to half-Cauchy prior, `log(1 + x^2/100)`.
#' @param evol0_sample logical; If `TRUE` (default), the prior variance of the initial `D` values of mu is sampled,
#' else it is fixed at `evol0_sd`
#' @param evol0_sd numeric; initial value (defaults to 10) or the prior standard deviation
#' for the initial `D` values of mu
#' @param sigma_e numeric; scale value (defaults to `1/sqrt(T)`) for half-Cauchy prior on global variance parameter
#' @param chol0 logical; If anything except `NULL` (the default), the Cholesky term of the log volatility is precomputed
#' @param offset length \code{T} vector of offset values for the log conditional expectation
#' @param verbose logical; If TRUE (the default), time remaining is printed to console
#' @param seed optional seed for random number generation for reproducible results
#'
#' @return
#' A named list with the following. One named element for each of the parameters specified in `mcmc_params` containing
#' the `nsave` posterior draws as a vector for scalar parameters and matrix for multivariate parameters.
#' An element named `loglike` containing a vector of the `nsave` evaluations of the log
#' likelihood of the data `y` with a Negative Binomial distribution defined by the current iteration conditional expectation + offset and overdispersion.
#' If `computeDIC` is `TRUE`, a named element for the computed value of DIC and effective parameters p_d for the data `y`
#' with a Negative Binomial distribution defined by the posterior mean of the conditional expectation + offset and posterior
#' mean overdispersion.
#'
#' @export
#'
#' @examples
#' beta <- bumps_sim()
#' y <- rnbinom(n = length(beta), size = 5, mu = beta)
#'
#' # Need to typically run sampler for longer than specified below
#' fit <- btf_nb(
#'   y = y,
#'   D = 2,
#'   nburn = 5000,
#'   evol0_sample = FALSE,
#'   verbose = FALSE,
#'   sigma_e = 1,
#'   chol0 = TRUE
#' )
btf_nb = function(y, evol_error = 'DHS', D = 2,
                  nsave = 1000, nburn = 1000, nskip = 4,
                  mcmc_params = list("mu", "yhat","evol_sigma_t2", "r", "dhs_phi", "dhs_mean"),
                  r_init = NULL, r_sample = "int_mh", step = 1,
                  mu_init = NULL, mu_sample = TRUE,
                  prior_r = rlang::expr(log(1 + x^2/100)), # half-Cauchy r prior
                  evol0_sample = TRUE, evol0_sd = 10, # should the initial variances for for the state vector be sampled and the fixed sd if not sampled
                  sigma_e = 1/sqrt(Nt), # What is prior variance for half-cauchy on tau
                  chol0 = NULL, computeDIC = TRUE, # TODO: bad default argument for chol0 but...
                  offset = 0, verbose = TRUE,
                  seed = NULL){

  set.seed(seed)

  # Convert to upper case:
  evol_error = toupper(evol_error)

  # For D = 0, return special case:
  # if(D == 0){
  #   return(btf0(y = y, evol_error = evol_error,
  #               nsave = nsave, nburn = nburn, nskip = nskip,
  #               mcmc_params = mcmc_params,
  #               computeDIC = computeDIC,
  #               verbose = verbose))
  # }

  # Time points (in [0,1])
  Nt = length(y); t01 = seq(0, 1, length.out=Nt);

  # Initialize bandsparse locations
  loc = dspCP::t_create_loc(length(y)-D, 1)
  # loc_obs = t_create_loc(length(y), D)

  # Begin by checking for missing values, then imputing (for initialization)
  is.missing = which(is.na(y)); any.missing = (length(is.missing) > 0)

  # Impute the active "data"
  y = approxfun(t01, y, rule = 2)(t01)

  # Initial overdispersion parameter
  if(is.null(r_init)){
    r <- 5 ## TODO: is there a better overdispersion starting value?
  }else{
    r <- r_init
  }

  # Compute the Cholesky term (uses random variances for a more conservative sparsity pattern)
  if(!is.null(chol0)) chol0 = dspCP::initChol.spam(T = Nt, D = D)

  # Initialize the conditional mean, mu:
  # mu <- matrix(log(y + 1), ncol = 1) #this initialization breaks things when omega is 0...

  if(is.null(mu_init)){
    mu <-  matrix(log(y + 1) + 0.1*rnorm(Nt), ncol = 1)
  }else{
    mu <- matrix(mu_init, ncol = 1)
  }


  # Initial SD (implicitly assumes a constant mean)

  if(is.null(r_sample) || r_sample == "int_mh"){
    eta_t <- pgdraw::pgdraw(y + r, mu + offset - log(r))
  }else{
    eta_t <- BayesLogit::rpg(Nt, y + r, mu + offset - log(r))
  }

  sigma_et <- sqrt(1/eta_t)

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  # And the initial states:
  mu0 = as.matrix(mu[1:D,])

  # Initialize the evolution error variance parameters:
  evolParams = dspCP::initEvolParams(omega, evol_error = evol_error)

  # Initial variance parameters:
  evolParams0 = dspCP::initEvol0(mu0)

  if(!evol0_sample){
    evolParams0$sigma_w0 <- rep(evol0_sd, D)
  }

  # For HS MCMC comparisons:
  # evolParams$dhs_phi = 0

  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu = array(NA, c(nsave, Nt))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, Nt))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post_r = array(NA, c(nsave))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, Nt))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  post_loglike = numeric(nsave)

  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler

  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = stats::rnbinom(length(is.missing), size = r, mu = mu[is.missing])

    # Sample the overdispersion:
    if(!is.null(r_sample)){
      r <- sample_r(y, r, exp(mu + offset), r_sample, step = step, prior_r = prior_r)
    }

    # Sample auxiliary variables
    if(is.null(r_sample) || r_sample == "int_mh"){
      eta_t <- pgdraw::pgdraw(y + r, mu + offset - log(r))
    }else{
      eta_t <- BayesLogit::rpg(Nt, y + r, mu + offset - log(r))
    }

    sigma_et <- sqrt(1/eta_t)

    # Sample the states:
    if(mu_sample){
      mu = sampleBTF_nb(y, r, offset = offset,
                        eta_t = eta_t, obs_sigma_t2 = sigma_et^2,
                        evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2),
                        D = D, chol0 = chol0)
    }

    # Compute the evolution errors:
    omega = diff(mu, differences = D)

    # And the initial states:
    mu0 = as.matrix(mu[1:D,])

    # Sample the initial variance parameters:
    ## is this used anywhere?! yes in the sampleBTF
    if(evol0_sample) evolParams0 = dspCP::sampleEvol0(mu0, evolParams0, A = 1)

    # Evolution error variance + params:
    evolParams = dspCP::sampleEvolParams(omega = omega, evolParams = evolParams,
                                  sigma_e = sigma_e, evol_error = evol_error,
                                  loc = loc)

    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        # Save the MCMC samples:
        if(!is.na(match('mu', mcmc_params)) || computeDIC) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = rnbinom(Nt, size = r, mu = exp(mu + offset))
        if(!is.na(match('r', mcmc_params)) || computeDIC) post_r[isave] = r
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnbinom(y, size = r, mu = exp(mu + offset), log = TRUE))

        # And reset the skip counter:
        skipcount = 0
      }
    }
    if(verbose) dspCP::computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }

  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post_r
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean

  # Also include the log-likelihood:
  mcmc_output$loglike = post_loglike

  if(computeDIC){
    # Log-likelihood evaluated at posterior means:
    loglike_hat = sum(stats::dnbinom(y,
                              size = mean(post_r),
                              mu = exp(colMeans(post_mu + offset)),
                              log = TRUE))

    # Effective number of parameters (Note: two options)
    p_d = c(2*(loglike_hat - mean(post_loglike)),
            2*var(post_loglike))
    # DIC:
    DIC = -2*loglike_hat + 2*p_d

    # Store the DIC and the effective number of parameters (p_d)
    mcmc_output$DIC = DIC; mcmc_output$p_d = p_d
  }

  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

  return (mcmc_output);
}

sampleBTF_nb <- function(y, r, offset, eta_t, obs_sigma_t2, evol_sigma_t2,
                         D = 1, chol0 = NULL){

  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')

  if(any(is.na(y))) stop('y cannot contain NAs')

  Nt = length(y)

  # Linear term:
  linht = eta_t*(log(r) - offset + 0.5*(y - r)/eta_t)

  # Quadratic terms and solutions are computed differently, depending on D:

  if(D == 0){
    # Special case: no differencing

    # Posterior SDs and posterior means:
    postSD = 1/sqrt(1/obs_sigma_t2 + 1/evol_sigma_t2)
    postMean = (linht)*postSD^2

    # Sample the states:
    mu = rnorm(n = Nt, mean = postMean, sd = postSD)

  } else {
    # All other cases: positive integer differencing (D = 1 or D = 2)

    # Quadratic term (D = 1 or D = 2)
    QHt_Matrix = dspCP::build_Q(obs_sigma_t2 = obs_sigma_t2, evol_sigma_t2 = evol_sigma_t2, D = D)

    if(!is.null(chol0)){
      # New sampler, based on spam package:

      # Sample the states:
      mu = matrix(spam::rmvnorm.canonical(n = 1,
                                    b = linht,
                                    Q = spam::as.spam.dgCMatrix(as(QHt_Matrix, "generalMatrix")),
                                    Rstruct = chol0))
    } else {
      # Original sampler, based on Matrix package:

      # Cholesky of Quadratic term:
      chQht_Matrix = Matrix::chol(QHt_Matrix) #this is upper triangular

      # Sample the states:
      mu = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(Nt)))
      # mu = matrix(chol2inv(chQht_Matrix)%*%(linht) + t(chQht_Matrix)%*%rnorm(Nt)) # TODO: check if this is doing what I want

    }
  }

  # And return the states:
  mu
}


#' Title
#'
#' @param y
#' @param d.prev
#' @param mu
#' @param r_sample
#' @param step
#' @param lambda_r
#' @param prior_r
#'
#' @return
#' @export
#'
#' @examples
sample_r <- function(y, d.prev, mu, r_sample, step = 1, lambda_r = 10, prior_r = rlang::expr(1 + x^2/100))
{

  # browser()

  if(r_sample == "int_mh"){
    # browser()

    d.new = d.prev

    rw.lower = max(d.prev - step, 1);
    rw.upper = d.prev + step;
    rw.grid  = rw.lower:rw.upper;
    rw.n     = length(rw.grid)
    rw.p     = rep(1/rw.n, rw.n);

    d.prop = sample(rw.grid, 1, prob=rw.p);

    ltarget = sum(dnbinom(y,
                          size = d.prop,
                          mu = mu,
                          log = TRUE)) + #proposal likelihood
      dpois(d.prop, lambda_r, log = TRUE) - # proposal prior
      sum(dnbinom(y,
                  size = d.prev,
                  mu = mu,
                  log = TRUE)) - # previous likelihood
      dpois(d.prev, lambda_r, log = TRUE) # previous prior
    lppsl = log(ifelse(d.prev - step <= 1, 1/(d.prev + step), 1/(2*step+1))) -
      log(ifelse(d.prop - step <= 1, 1/(d.prop + step), 1/(2*step+1)));
    lratio = ltarget + lppsl # (llik.prop + lprior.prop + ll.J(d.prev)) - (llik.prev + + lprior.prev + ll.J(d.prop))

    if (runif(1) < exp(lratio)) d.new = d.prop

    return(d.new)
  }

  ## Kowal Slice sampler
  if(r_sample == "slice"){

    # e0 <- 1/sigma_c^2 # variance of the half cauchy (10^2)

    r = dspCP::uni.slice(d.prev, g = function(x){
      sum(dnbinom(y, size = x, mu = mu, log = TRUE)) + eval(prior_r)
    }, lower = 0, upper = 1000)

    return(r)
  }

  ## BayesLogit MH
  if(r_sample == "mh"){

    r = d.prev

    if(d.prev > 1) rstar = runif(1,d.prev-1,d.prev+1)
    else rstar = runif(1,0,2)
    ll = sum(dnbinom(y, d.prev, mu = mu, log=TRUE)) + # previous likelihood
      dgamma(1/d.prev, 1, 1, log = TRUE) # previous prior
    llstar = sum(dnbinom(y, rstar, mu = mu, log=TRUE)) + # proposal likelihood
      dgamma(1/rstar, 1, 1, log = TRUE) # proposal prior
    lalpha = llstar - ll
    ## TODO: add prior for this sampler

    lu = log(runif(1))

    # if(is.na(lalpha)){browser()}

    if(lu < lalpha) r = rstar

    return(r)

  }

}

#' Generate Trend Functions
#'
#' Two helper functions for calculating trend functions used to test various trend
#' filtering methods. The functions are modified variations of the wavelet test functions bumps and doppler.
#'
#' @param n integer with default to 1024 of the length of time series; increasing n infills
#'
#' @return vector of length `n`
#' @export
#'
#' @examples
#'
#' beta_bumps <- bumps_sim()
#' y <- rnbinom(n = length(beta_bumps), size = 5, mu = beta_bumps)
#'
#'
bumps_sim <- function(n=1024){

  x <- (0:(n-1.))/n

  pos <- c(0.1, 0.15, 0.25, 0.4, 0.65, 0.78)
  hgt <- c(4, 3, 5, 4.2, 4.3, 5.1)
  wth <- 2*c(0.01, 0.012, 0.02, 0.06, 0.02, 0.016)
  z <- rep(1, n)
  for(j in 1:length(pos)) z <- z+hgt[j]/(1+abs((x-pos[j]))/wth[j])^4

  return(z)

}


#' @export
#'
#' @examples
#' beta_dop <- bumps_sim()
#' y <- rnbinom(n = length(beta_dop), size = 5, mu = exp(beta_dop))
#'
#' @rdname bumps_sim
dop_sim <- function(n = 1024){
  x <- (0:(n-1))/n

  ifelse(x <= 1/6,
         -sqrt(x*(1-x)),
         sqrt(x*(1-x))*sin((2*pi*1.05)/(x+.05)))
}
