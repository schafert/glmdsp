#' Title
#'
#' @param y
#' @param D
#' @param useAnom
#' @param useObsSV
#' @param nsave
#' @param nburn
#' @param nskip
#' @param mcmc_params
#' @param r_init
#' @param r_sample
#' @param step
#' @param mu_init
#' @param mu_sample
#' @param prior_r
#' @param evol0_sample
#' @param evol0_sd
#' @param thresh_sample
#' @param thresh_init
#' @param thresh_proposal
#' @param thresh_proposal_sd
#' @param omega_trans
#' @param chol0
#' @param computeDIC
#' @param offset
#' @param verbose
#' @param cp_thres
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
abco_nb <- function(y, D = 1, useAnom=FALSE, useObsSV = FALSE,
                    nsave = 1000, nburn = 1000, nskip = 4,
                    mcmc_params = list("mu", "yhat","evol_sigma_t2", "r", "dhs_phi",
                                       "dhs_phi2", "dhs_mean", "gamma", "omega"),
                    r_init = NULL, r_sample = "int_mh", step = 1,
                    mu_init = NULL, mu_sample = TRUE,
                    prior_r = rlang::expr(log(1 + x^2/100)), # half-Cauchy r prior
                    evol0_sample = TRUE, evol0_sd = 10, # should the initial variances for for the state vector be sampled and the fixed sd if not sampled
                    thresh_sample = TRUE, thresh_init = NULL, # should the threshold be sampled? Is it initialized?
                    thresh_proposal = "unif", thresh_proposal_sd = 0.1, # threshold proposal
                    omega_trans = rlang::expr(log(omega^2+0.0001)), # transformation for threshold
                    chol0 = NULL, computeDIC = TRUE,
                    offset = 0,
                    verbose = TRUE, cp_thres = 0.5, seed = NULL){

  set.seed(seed)

  # Time points (in [0,1])
  Nt = length(y); t01 = seq(0, 1, length.out=Nt)
  evol_error = 'DHS'

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
  if(!is.null(chol0)) chol0 = initChol.spam(T = Nt, D = D)
  loc = dspCP::t_create_loc(length(y)-D, 1)
  loc_obs = dspCP::t_create_loc(length(y), D)

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

  if (useAnom){
    ## TODO: add in useAnom later
    stop('useAnom sampler not available yet')
    # zeta = t_sampleBTF(y[-(1:D)]-mu[-(1:D)], obs_sigma_t2 = sigma_e[-(1:D)]^2, evol_sigma_t2 = 0.01*sigma_e[-(1:D)]^2, D = 0, loc_obs)
    # zeta_params = t_initEvolZeta_ps(zeta)
  }

  # Compute the evolution errors:
  omega = diff(mu, differences = D)

  evolParams = t_initEvolParams_no2(log(y + 0.0001), D, omega, omega_trans)
  lower_b = evolParams$lower_b
  upper_b = evolParams$upper_b

  if(!is.null(thresh_init)){
    evolParams$r <- thres_init
  }

  # And the initial states:
  mu0 = as.matrix(mu[1:D,])

  # Initial variance parameters:
  evolParams0 = dspCP::initEvol0(mu0)

  if(!evol0_sample){
    evolParams0$sigma_w0 <- rep(evol0_sd, D)
  }

  #Array to store omega
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', length(mcmc_params)); names(mcmc_output) = mcmc_params
  if(!is.na(match('mu', mcmc_params))) post_mu = array(NA, c(nsave, Nt))
  if(!is.na(match('zeta', mcmc_params))) post_zeta = array(NA, c(nsave, Nt-D))
  if(!is.na(match('yhat', mcmc_params))) post_yhat = array(NA, c(nsave, Nt))
  if(!is.na(match('r', mcmc_params)) || computeDIC) post_r = array(NA, c(nsave))
  if(!is.na(match('h', mcmc_params))) post_h = array(NA, c(nsave, Nt-D))
  if(!is.na(match('gamma', mcmc_params)) && evol_error == "DHS") post_gamma = numeric(nsave)
  if(!is.na(match('zeta_sigma_t2', mcmc_params))) post_zeta_sigma_t2 = array(NA, c(nsave, Nt-D))
  if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2 = array(NA, c(nsave, Nt))
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi = numeric(nsave)
  if(!is.na(match('dhs_phi2', mcmc_params)) && evol_error == "DHS") post_dhs_phi2 = numeric(nsave)
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean = numeric(nsave)
  if(!is.na(match('omega', mcmc_params))) post_omega = array(NA, c(nsave, Nt-D))
  post_loglike = numeric(nsave)

  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting

  # Run the MCMC:
  if(verbose) timer0 = proc.time()[3] # For timing the sampler
  for(nsi in 1:nstot){

    # Impute missing values, if any:
    if(any.missing) y[is.missing] = stats::rnbinom(length(is.missing), size = r, mean = exp(mu[is.missing]))

    # Sample the overdispersion:
    if(!is.null(r_sample)){
      ## TODO: probably need to add the anomaly detection stuff here if going to add that in
      r <- sample_r(y, r, exp(mu + offset),
                    r_sample, step = step, prior_r = prior_r)
    }

    # Sample auxiliary variables
    if(r_sample == "int_mh"){
      eta_t <- pgdraw::pgdraw(y + r, mu + offset - log(r))
    }else{
      eta_t <- BayesLogit::rpg(Nt, y + r, mu + offset - log(r))
    }

    sigma_et <- sqrt(1/eta_t)

    # Sample the states:
    if (useAnom){
      ## TODO: add in useAnom later
      stop('useAnom sampler not available yet')
      # mu = t_sampleBTF_nb(y-c(rep(0, D), zeta), obs_sigma_t2 = sigma_e^2, evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2), D = D, loc_obs)
      # zeta = t_sampleBTF_nb(y[-(1:D)]-mu[-(1:D)], obs_sigma_t2 = sigma_et[-(1:D)]^2, evol_sigma_t2 = zeta_params$sigma_wt^2, D = 0)
      # zeta_params = t_sampleEvolZeta_ps(zeta, zeta_params)
    } else{
      mu = t_sampleBTF_nb(y, r, offset = offset,
                          obs_sigma_t2 = sigma_et^2,
                          eta_t = eta_t,
                          evol_sigma_t2 = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2),
                          D = D, chol0 = chol0)
    }

    # Compute the evolution errors:
    omega = diff(mu, differences = D)

    # And the initial states:
    mu0 = as.matrix(mu[1:D,])

    # Sample the initial variance parameters:
    if(evol0_sample) evolParams0 = dspCP::sampleEvol0(mu0, evolParams0, A = 1)

    # Evolution error variance + params:
    evolParams <- t_sampleEvolParams2(omega = omega, evolParams = evolParams,
                                      D = D, sigma_e = 1/sqrt(Nt),
                                      lower_b = lower_b, upper_b = upper_b,
                                      loc = loc,
                                      thresh_sample = thresh_sample,
                                      thresh_proposal = thresh_proposal,
                                      thresh_proposal_sd = thresh_proposal_sd,
                                      omega_trans = omega_trans)

    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1

      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1

        if(!is.na(match('mu', mcmc_params))) post_mu[isave,] = mu
        if(!is.na(match('yhat', mcmc_params))) post_yhat[isave,] = rnbinom(Nt, size = r, mu = exp(mu + offset))
        if(!is.na(match('r', mcmc_params)) || computeDIC) post_r[isave] = r
        if (useAnom) {
          if(!is.na(match('zeta', mcmc_params))) post_zeta[isave,] = zeta
          if(!is.na(match('zeta_sigma_t2', mcmc_params))) post_zeta_sigma_t2[isave,] = zeta_params$sigma_wt^2
        }
        if(!is.na(match('h', mcmc_params))) post_h[isave,] = evolParams$ht
        if(!is.na(match('omega', mcmc_params))) post_omega[isave,] = omega
        if(!is.na(match('gamma', mcmc_params)) && evol_error == "DHS") post_gamma[isave] = evolParams$r
        if(!is.na(match('evol_sigma_t2', mcmc_params))) post_evol_sigma_t2[isave,] = c(evolParams0$sigma_w0^2, evolParams$sigma_wt^2)
        if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") post_dhs_phi[isave] = evolParams$dhs_phi
        if(!is.na(match('dhs_phi2', mcmc_params)) && evol_error == "DHS") post_dhs_phi2[isave] = evolParams$dhs_phi2
        if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") post_dhs_mean[isave] = evolParams$dhs_mean
        post_loglike[isave] = sum(dnbinom(y, size = r, mu = exp(mu + offset), log = TRUE))

        skipcount = 0
      }
    }
    if(verbose) computeTimeRemaining(nsi, timer0, nstot, nrep = 1000)
  }
  if(!is.na(match('zeta', mcmc_params))) mcmc_output$zeta = post_zeta
  if(!is.na(match('mu', mcmc_params))) mcmc_output$mu = post_mu
  if(!is.na(match('yhat', mcmc_params))) mcmc_output$yhat = post_yhat
  if(!is.na(match('r', mcmc_params))) mcmc_output$r = post_r
  if(!is.na(match('gamma', mcmc_params)) && evol_error == "DHS") mcmc_output$gamma = post_gamma
  if(!is.na(match('h', mcmc_params))) mcmc_output$h = post_h
  if(!is.na(match('zeta_sigma_t2', mcmc_params))) mcmc_output$zeta_sigma_t2 = post_zeta_sigma_t2
  if(!is.na(match('evol_sigma_t2', mcmc_params))) mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(!is.na(match('dhs_phi', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi = post_dhs_phi
  if(!is.na(match('dhs_phi2', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_phi2 = post_dhs_phi2
  if(!is.na(match('dhs_mean', mcmc_params)) && evol_error == "DHS") mcmc_output$dhs_mean = post_dhs_mean
  if(!is.na(match('omega', mcmc_params))) mcmc_output$omega = post_omega

  if(verbose) print(paste('Total time: ', round((proc.time()[3] - timer0)), 'seconds'))

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


  mcmc_output$cp = identify_cp_nb(D, mcmc_output, omega_trans, cp_thres)
  return (mcmc_output)
}

t_sampleBTF_nb <- function(y, r, offset, eta_t, obs_sigma_t2,
                           evol_sigma_t2,
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
                                    Q = spam::as.spam.dgCMatrix(as(QHt_Matrix, "dgCMatrix")),
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

identify_cp_nb = function(D, mcmc_output, omega_trans, cp_thres= 0.5){
  cp_list = rep(0, length(mcmc_output$omega[1,])+D)
  for (j in 1:length(mcmc_output$omega[,1])){
    nz_list = which(eval(omega_trans, list(omega = mcmc_output$omega[j,])) > mcmc_output$gamma[j])+D
    for (k in nz_list){
      cp_list[k] = cp_list[k]+1
    }
  }
  cp_list = cp_list / length(mcmc_output$omega[,1])

  which(cp_list >= cp_thres)
}

dic_nb <- function(y, mcmc_output){

  niter <- length(mcmc_output$r)

  post_loglike <- c()

  for(j in 1:niter){
    post_loglike[j] <- sum(stats::dnbinom(y, size = mcmc_output$r[j], mu = exp(mcmc_output$mu[j,]), log = TRUE))
  }

  loglike_hat = sum(stats::dnbinom(y,
                            size = mean(mcmc_output$r),
                            mu = exp(colMeans(mcmc_output$mu)),
                            log = TRUE))

  # Effective number of parameters (Note: two options)
  p_d = c(2*(loglike_hat - mean(post_loglike)),
          2*var(post_loglike))
  # DIC:
  DIC = -2*loglike_hat + 2*p_d

  # Store the DIC and the effective number of parameters (p_d)
  list(DIC = DIC, p_d = p_d)
}

t_initEvolParams_no2 = function(y, D, omega, omega_trans = rlang::expr(log(omega^2+0.0001))){

  # "Local" number of time points
  yoffset = tcrossprod(rep(1,length(omega)),
                       apply(as.matrix(omega), 2,
                             function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
  ht = log(omega^2+0.0001)
  n = length(omega)

  # Initialize the AR(1) model to obtain unconditional mean and AR(1) coefficient
  arCoefs = arima(ht, c(1,0,0), method="ML")$coef
  dhs_mean = arCoefs[2]; dhs_phi = arCoefs[1]
  dhs_phi2 = msm::rtnorm(1, mean = -0.5, sd = .5, upper = 0, lower = -1)

  # Initialize the SD of log-vol innovations simply using the expectation:
  sigma_eta_t = matrix(pi, nrow = n-1, ncol = 1)
  sigma_eta_0 = rep(pi) # Initial value

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  diff_y = eval(omega_trans, list(omega = diff(y, differences=D)))
  lower_b = log(0.0001)
  upper_b = max(diff_y)
  r = runif(1, min=lower_b, max = upper_b)

  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, dhs_phi2 = dhs_phi2, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, r=r, lower_b = lower_b, upper_b = upper_b)
}

t_sampleEvolParams2 = function(omega, evolParams, D = 1, sigma_e = 1,
                               lower_b, upper_b, loc,
                               prior_dhs_phi = c(20,1), alphaPlusBeta = 1,
                               thresh_proposal = "unif", thresh_proposal_sd = 0.1,
                               thresh_sample = TRUE, omega_trans = rlang::expr(log(omega^2+0.0001))){
  # Store the DSP parameters locally:
  ht = evolParams$ht; dhs_mean = evolParams$dhs_mean; dhs_phi = evolParams$dhs_phi; dhs_phi2 = evolParams$dhs_phi2
  sigma_eta_t = evolParams$sigma_eta_t; sigma_eta_0 = evolParams$sigma_eta_0; r=evolParams$r

  # "Local" number of time points
  n = length(ht); p = 1
  st = (log(omega^2 + 0.0001) >= r)

  # Sample the log-volatilities using AWOL sampler
  ht = t_sampleLogVols(h_y = omega, h_prev = ht, h_mu = dhs_mean, h_phi=dhs_phi, h_phi2 = dhs_phi2, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_st = st, loc = loc)

  # Compute centered log-vols for the samplers below:
  ht_tilde = ht - dhs_mean

  # Sample AR(1) parameters
  # Note: dhs_phi = 0 means non-dynamic HS, while dhs_phi = 1 means RW, in which case we don't sample either
  coef = dspCP::t_sampleAR1(h_yc = ht_tilde, h_phi = dhs_phi, h_phi2 = dhs_phi2, h_sigma_eta_t = sigma_eta_t, h_st = st, prior_dhs_phi = prior_dhs_phi)
  dhs_phi = coef[1]
  dhs_phi2 = coef[2]

  # Sample the unconditional mean(s), unless dhs_phi = 1 (not defined)
  dhs_mean = dspCP::t_sampleLogVolMu(h = ht, h_mu = dhs_mean, h_phi = dhs_phi, h_phi2 = dhs_phi2, h_sigma_eta_t = sigma_eta_t, h_sigma_eta_0 = sigma_eta_0, h_st = st, h_log_scale = log(sigma_e^2));

  # Sample the evolution error SD of log-vol (i.e., Polya-Gamma mixing weights)
  eta_t = ht_tilde[-1] - (dhs_phi+ dhs_phi2*st[-n])*ht_tilde[-n]       # Residuals
  #sigma_eta_t = matrix(1/sqrt(rpg(num = (n-1), h = alphaPlusBeta, z = eta_t)), nc = 1) # Sample
  #sigma_eta_0 = 1/sqrt(rpg(num = 1, h = 1, z = ht_tilde[1]))                # Sample the inital
  sigma_eta_t = 1/sqrt(pgdraw::pgdraw(alphaPlusBeta, eta_t))
  sigma_eta_0 = 1/sqrt(pgdraw::pgdraw(1, ht_tilde[1]))

  if(thresh_sample){
    r = t_sampleR_mh2(h_yc = ht_tilde, h_phi = dhs_phi, h_phi2 = dhs_phi2,
                     h_sigma_eta_t = sigma_eta_t,
                     h_sigma_eta_0 = sigma_eta_0, h_st = st, h_r = r,
                     lower_b = lower_b, upper_b = upper_b, omega = omega, D = D,
                     proposal = thresh_proposal, proposal_sd = thresh_proposal_sd,
                     omega_trans = omega_trans)
  }

  # Evolution error SD:
  sigma_wt = exp(ht/2)

  # Return the same list, but with the new values
  list(sigma_wt = sigma_wt, ht = ht, dhs_mean = dhs_mean, dhs_phi = dhs_phi, dhs_phi2 = dhs_phi2, sigma_eta_t = sigma_eta_t, sigma_eta_0 = sigma_eta_0, r=r)
}

#----------------------------------------------------------------------------
#' Sample the threshold parameter
#'
#' Compute one draw of thee threshold parameter in th TAR(1) model with Gaussian innovations
#' and time-dependent innovation variances. The sampler utilizes metropolis hasting to draw
#' from uniform prior.
#'
#' @param h_yc the \code{T} vector of centered log-volatilities
#' (i.e., the log-vols minus the unconditional means \code{dhs_mean})
#' @param h_phi the \code{1} vector of previous AR(1) coefficient(s)
#' @param h_phi2 the \code{1} vector of previous penalty coefficient(s)
#' @param h_sigma_eta_t the \code{T} vector of log-vol innovation standard deviations
#' @param h_sigma_eta_0 the \code{1} vector of initial log-vol innovation standard deviations
#' @param h_st the \code{T} vector of indicators on whether each time-step exceed the estimated threshold
#' @param h_r \code{1} the previous draw of the threshold parameter
#' @param upper_b the upper bound in the uniform prior of the threshold variable
#' @param lower_b the lower bound in the uniform prior of the threshold variable
#' @param omega \code{T} vector of evolution errors
#' @param D the degree of differencing (one or two)
#' @param proposal string indicating whether or not to propose new value from random walk "rw"
#' @param proposal_sd numeric value of standard deviation of random walk proposal
#' @param omega_trans expression for the transformation of omega used to decide thresholding
#'
#' @return the sampled threshold value \code{r}
t_sampleR_mh2 = function(h_yc, h_phi, h_phi2, h_sigma_eta_t, h_sigma_eta_0,
                        h_st, h_r, lower_b, upper_b, omega, D,
                        proposal = "unif", proposal_sd = 0.1,
                        omega_trans = rlang::expr(log(omega^2+0.0001))){
  n = length(h_yc);

  y_ar = h_yc[-1]/h_sigma_eta_t
  x_ar = h_yc[-n]/h_sigma_eta_t

  if(proposal == "rw"){
    new_cand = h_r + proposal_sd*rnorm(1)
  }else{
    new_cand = runif(1, lower_b, upper_b)
  }

  num = sum(dnorm(y_ar, (h_phi + h_phi2*(eval(omega_trans)[-n] >= new_cand))*x_ar, 1, log=TRUE))
  den = sum(dnorm(y_ar, (h_phi + h_phi2*(eval(omega_trans)[-n] >= h_r))*x_ar, 1, log=TRUE))

  if ((num > den) & (new_cand >= lower_b) & (new_cand <= upper_b)){ # needs to be in support of prior
    r = new_cand
  }else{
    prob = runif(1, 0, 1)
    if ((prob < exp(num-den)) & (new_cand >= lower_b) & (new_cand <= upper_b)){
      r = new_cand
    } else{
      r = h_r
    }
  }

  r
}

t_sampleLogVols = function(h_y, h_prev, h_mu, h_phi, h_phi2, h_sigma_eta_t, h_sigma_eta_0, h_st, loc){

  # Compute dimensions:
  n = length(h_prev); p = 1

  # Mixture params: mean, variance, and weights
  # Kim, Shephard, Chib (1998) 7-component mixture:
  #m_st  = c(-11.40039, -5.24321, -9.83726, 1.50746,  -0.65098, 0.52478,  -2.35859)
  #v_st2 = c(5.795960,  2.613690, 5.179500, 0.167350, 0.640090, 0.340230, 1.262610)
  #q     = c(0.007300,  0.105560, 0.000020, 0.043950, 0.340010, 0.245660, 0.257500)

  # Omori, Chib, Shephard, Nakajima (2007) 10-component mixture:
  m_st  = c(1.92677, 1.34744, 0.73504, 0.02266, -0.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000)
  v_st2 = c(0.11265, 0.17788, 0.26768, 0.40611,  0.62699,  0.98583,  1.57469,  2.54498,  4.16591,   7.33342)
  q     = c(0.00609, 0.04775, 0.13057, 0.20674,  0.22715,  0.18842,  0.12047,  0.05591,  0.01575,   0.00115)

  # Add an offset: common for all times, but distict for each j=1,...,p
  #yoffset = tcrossprod(rep(1,n),
  #                     apply(as.matrix(h_y), 2,
  #                           function(x) any(x^2 < 10^-16)*max(10^-8, mad(x)/10^6)))
  yoffset = 10^-8

  # This is the response in our DLM, log(y^2)
  ystar = log(h_y^2 + yoffset)

  # Sample the mixture components
  # z = sapply(ystar-h_prev, ncind, m_st, sqrt(v_st2), q)
  #z = sapply(ystar-h_prev, ncind, m_st, sqrt(v_st2), q)
  z = c(dspCP::draw_indicators_generic(ystar-h_prev, rep(0, length(ystar)), length(ystar),
                                q, m_st, sqrt(v_st2), length(q)))

  # Subset mean and variances to the sampled mixture components; (n x p) matrices
  m_st_all = m_st[z]
  v_st2_all = v_st2[z]

  # Joint AWOL sampler for j=1,...,p:

  # Constant (but j-specific) mean
  #h_mu_all = (1-h_phi-h_phi2*h_st)*h_mu
  #h_mu_all[1] = h_mu
  #h_mu_all = tcrossprod(rep(1,n), h_mu)

  # Constant (but j-specific) AR(1) coef
  #h_phi_all = tcrossprod(rep(1,n), h_phi)
  #h_phi2_all = tcrossprod(rep(1,n), h_phi2)

  # Linear term:
  linht = (ystar - m_st_all - h_mu)/v_st2_all

  # Evolution precision matrix (n x p)
  evol_prec_mat = rep(0, n);
  evol_prec_mat[1] = 1/h_sigma_eta_0^2;
  evol_prec_mat[-1] = 1/h_sigma_eta_t^2;

  # Lagged version, with zeros as appropriate (needed below)
  evol_prec_lag_mat = rep(0, n)
  evol_prec_lag_mat[1:(n-1)] = evol_prec_mat[-1]

  # Diagonal of quadratic term:
  Q_diag = 1/v_st2_all +  evol_prec_mat + (h_phi+h_phi2*h_st)^2*evol_prec_lag_mat

  # Off-diagonal of quadratic term:
  Q_off = (-(h_phi+h_phi2*h_st)*evol_prec_lag_mat)[-n]

  # Quadratic term:
  #QHt_Matrix = bandSparse(n*p, k = c(0,1), diag = list(Q_diag, Q_off), symm = TRUE)

  # Cholesky:
  #chQht_Matrix = Matrix::chol(QHt_Matrix)

  # Sample the log-vols:
  #hsamp = h_mu_all + matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(length(linht))), nr = n)
  rd = zrnorm(length(linht))
  hsamp = h_mu + dspCP::sample_mat(loc$r, loc$c, c(Q_diag, Q_off, Q_off), length(Q_diag), length(loc$r), c(linht), rd, 1)

  # Return the (uncentered) log-vols
  hsamp
}
