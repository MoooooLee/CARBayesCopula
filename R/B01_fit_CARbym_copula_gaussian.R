#### Function of fitting Bivariate BYM CAR Model - Binomial ####
#### with Gaussian Copula ####

fit_CARbym_copula_gaussian <- function(formula,
                                    data = NULL,
                                    trials,
                                    W,
                                    burnin = 5000,
                                    n_sample = 10000,
                                    thin = 5,
                                    prior_beta_mean = NULL,
                                    prior_beta_var = NULL,
                                    prior_tau2 = NULL,
                                    prior_sigma2 = NULL,
                                    tau2_AB = NULL,
                                    sigma2_AB = NULL,
                                    verbose = TRUE) {
  # 2 Data preparing --------------------------------------------------------

  # 2.1 Frame Object --------------------------------------------------------

  a <- common_verbose(verbose)

  frame_results <- common_frame(formula, data, "binomial")
  K <- frame_results$n
  p <- frame_results$p
  X <- frame_results$X
  X_standardised <- frame_results$X.standardised
  X_sd <- frame_results$X.sd
  X_mean <- frame_results$X.mean
  X_indicator <- frame_results$X.indicator
  offset <- frame_results$offset
  Y <- frame_results$Y
  which_miss <- frame_results$which.miss
  n_miss <- frame_results$n.miss
  J <- ncol(Y)
  N_all <- K * J
  Y_DA <- Y


  # 2.2 Create missing list -------------------------------------------------

  if(n_miss>0){
    miss_locator <- array(NA, c(n_miss, 2))
    colnames(miss_locator) <- c("row", "column")
    locations <- which(t(which_miss)==0)
    miss_locator[ ,1] <- ceiling(locations/J)
    miss_locator[ ,2] <- locations - (miss_locator[ ,1]-1) * J
  }else{}


  # 2.3 Check input parameters ----------------------------------------------
  # 2.3.1 Check and format trials -------------------------------------------

  if(ncol(trials)!=J) stop("trials has the wrong number of columns_", call_=FALSE)
  if(nrow(trials)!=K) stop("trials has the wrong number of rows_", call_=FALSE)
  if(sum(is.na(trials))>0) stop("the numbers of trials has missing 'NA' values_", call_=FALSE)
  if(!is.numeric(trials)) stop("the numbers of trials has non-numeric values_", call_=FALSE)
  int_check <- N_all-sum(ceiling(trials)==floor(trials))
  if(int_check > 0) stop("the numbers of trials has non-integer values_", call_=FALSE)
  if(min(trials)<=0) stop("the numbers of trials has zero or negative values_", call_=FALSE)

  failures <- trials - Y
  failures_DA <- trials - Y_DA
  if(sum(Y>trials, na.rm=TRUE)>0){
    stop("the response variable has larger values that the numbers of trials_", call_=FALSE)
  }

  # if(!is.logical(Copula_Inv)) stop("the indicator of Inverse of Copula should be logical value_", call_=FALSE)

  # 2.3.2 Check W -----------------------------------------------------------

  if(!is.matrix(W)) stop("W is not a matrix_", call_=FALSE)
  if(ceiling(N_all/K)!= floor(N_all/K)){
    stop("The number of data points divided by the number of rows in W is not a whole number_",
         call_=FALSE)
  }


  # 2.3.3 Check and format rho ----------------------------------------------

  # if(is.null(rho))
  # {
  #   rho <- runif(1)
  #   fix_rho <- FALSE
  # }else
  # {
  #   fix_rho <- TRUE
  # }
  # if(!is.numeric(rho) ) stop("rho is fixed but is not numeric_", call_=FALSE)
  # if(rho<0 ) stop("rho is outside the range [0, 1]_", call_=FALSE)
  # if(rho>1 ) stop("rho is outside the range [0, 1]_", call_=FALSE)


  # 2.3.4 add Check the kendall's tau (alpha) -------------------------------

  if(is.null(tau2_AB))
  {
    tau2_AB <- runif(1, min = -1, max = 1)
    fix_tau2_AB <- FALSE
  }else
  {
    fix_tau2_AB <- TRUE
  }
  tau2_AB <- tau2_AB
  if(!is.numeric(tau2_AB) ) stop("tau2_AB is fixed but is not numeric_",
                                 call. = FALSE)
  if(tau2_AB < -1 | tau2_AB > 1) stop("tau2_AB is outside the range of correlation [-1,1].",
                                      call. = FALSE)

  if(is.null(sigma2_AB))
  {
    sigma2_AB <- runif(1, min = -1, max = 1)
    fix_sigma2_AB <- FALSE
  }else
  {
    fix_sigma2_AB <- TRUE
  }
  sigma2_AB <- sigma2_AB
  if(!is.numeric(sigma2_AB) ) stop("sigma2_AB is fixed but is not numeric_",
                                   call. = FALSE)
  if(sigma2_AB < -1 | sigma2_AB > 1) stop("tau2_AB is outside the range of correlation [-1,1].",
                                          call. = FALSE)


  # 2.3.5 Check and format prior --------------------------------------------

  if(is.null(prior_beta_mean)) prior_beta_mean <- rep(0, p)
  if(is.null(prior_beta_var)) prior_beta_var <- rep(100000, p)

  if(is.null(prior_tau2)) prior_tau2 <- c(1, 0.01)
  if(is.null(prior_sigma2)) prior_sigma2 <- c(1, 0.01)
  # if(is.null(prior_Sigma_df)) prior_Sigma_df <- J+1
  # if(is.null(prior_Sigma_scale)) prior_Sigma_scale <- diag(rep(1,J)) / 1000

  # if(is.null(prior_Tau_df)) prior_Tau_df <- J+1
  # if(is.null(prior_Tau_scale)) prior_Tau_scale <- diag(rep(1,J)) / 1000

  common_prior_beta_check(prior_beta_mean, prior_beta_var, p)
  # common_prior_varmat_check(prior_Sigma_scale, J)
  # common_prior_varmat_check(prior_Tau_scale, J)


  # 2.3.6 Separate beta into blocks -----------------------------------------

  block_temp <- common_betablock(p)
  beta_beg  <- block_temp[[1]]
  beta_fin <- block_temp[[2]]
  n_beta_block <- block_temp[[3]]
  list_block <- as.list(rep(NA, n_beta_block*2))
  for(r in 1:n_beta_block)
  {
    list_block[[r]] <- beta_beg[r]:beta_fin[r]-1
    list_block[[r+n_beta_block]] <- length(list_block[[r]])
  }


  # 2.3.7 Checking input of MCMC parameters ---------------------------------

  common_burnin_nsample_thin_check(burnin, n_sample, thin)


  # 3 Initial parameter values ----------------------------------------------

  beta <- array(NA, c(p, J))
  for(i in 1:J) {
    mod_glm <- glm(cbind(Y[ ,i], failures[ ,i])~X_standardised-1,
                   offset=offset[ ,i],
                   family="quasibinomial")
    beta_mean <- mod_glm$coefficients
    beta_sd <- sqrt(diag(summary(mod_glm)$cov.scaled))
    beta[ ,i] <- rnorm(n=p, mean=beta_mean, sd=beta_sd)
  }

  theta_hat <- Y / trials
  theta_hat[theta_hat==0] <- 0.01
  theta_hat[theta_hat==1] <- 0.99
  res_temp <- log(theta_hat / (1 - theta_hat)) - X_standardised %*% beta - offset
  res_sd <- sd(res_temp, na.rm=TRUE)/5

  phi_vec <- rnorm(n=N_all, mean=0, sd=res_sd)
  phi <- matrix(phi_vec, nrow=K, byrow=TRUE)
  tau2 <- diag(cov(phi))
  # Tau <- cov(phi)
  # Tau_inv <- cov(Tau)

  theta_vec <- rnorm(n=N_all, mean=0, sd=res_sd)
  theta <- matrix(theta_vec, nrow=K, byrow=TRUE)
  sigma2 <- diag(cov(theta))
  # Sigma <- cov(phi)
  # Sigma_inv <- solve(Sigma)

  regression <- X_standardised %*% beta
  lp <- regression + phi + theta + offset
  prob <- exp(lp)  / (1 + exp(lp))


  # 4 Set up MCMC quantities ------------------------------------------------

  # 4.1 Store samples -------------------------------------------------------

  n_keep <- floor((n_sample - burnin)/thin)
  samples_beta <- array(NA, c(n_keep, J*p))

  samples_phi <- array(NA, c(n_keep, N_all))
  # samples_Tau <- array(NA, c(n_keep, J, J))
  samples_tau2 <- array(NA, c(n_keep, J))

  samples_theta <- array(NA, c(n_keep, N_all))
  # samples_Sigma <- array(NA, c(n_keep, J, J))
  samples_sigma2 <- array(NA, c(n_keep, J))

  # if(!fix_rho) samples_rho <- array(NA, c(n_keep, 1))
  if(!fix_tau2_AB) samples_tau2_AB <- array(NA, c(n_keep, 1))
  if(!fix_sigma2_AB) samples_sigma2_AB <- array(NA, c(n_keep, 1))

  samples_loglike <- array(NA, c(n_keep, N_all))
  samples_fitted <- array(NA, c(n_keep, N_all))
  if(n_miss>0) samples_Y <- array(NA, c(n_keep, n_miss))


  # 4.2 Metropolis quantities -----------------------------------------------

  accept_beta <- rep(0,2*J)
  accept_phi <- rep(0,2)
  accept_theta <- rep(0,2)

  accept_tau2 <- rep(0,2)
  accept_sigma2 <- rep(0,2)

  # accept_rho <- rep(0,2)
  accept_tau2_AB <- rep(0,2)
  accept_sigma2_AB <- rep(0,2)

  proposal_sd_beta <- rep(0.01, J)
  proposal_sd_phi <- 0.1
  proposal_sd_theta <- 0.1

  proposal_sd_tau2 <- res_sd^2 / 2
  proposal_sd_sigma2 <- res_sd^2 / 2

  # proposal_sd_rho <- 0.02

  # Sigma_post_df <- prior_Sigma_df + K
  # Tau_post_df <- prior_Tau_df + K
  proposal_sd_tau2_AB <- 0.2
  proposal_sd_sigma2_AB <- 0.2


  # 5 Set up spatial quantities ---------------------------------------------

  # 5.1 CAR quantities ------------------------------------------------------

  W_quants <- common_Wcheckformat(W)
  W <- W_quants$W
  W_triplet <- W_quants$W.triplet
  n_triplet <- W_quants$n.triplet
  W_triplet_sum <- W_quants$W.triplet.sum
  n_neighbours <- W_quants$n.neighbours
  W_begfin <- W_quants$W.begfin

  Wstar <- diag(apply(W,1,sum)) - W
  Q <- 1 * Wstar + diag(rep(1-1,K))

  # 5.1.2 Create Identity Matrix for sampling theta -------------------------

  I <- diag(1, nrow = K, ncol = K)

  I_quants <- common_Wcheckformat(I)
  I <- I_quants$W
  I_triplet <- I_quants$W.triplet
  I_n_triplet <- I_quants$n.triplet
  I_triplet_sum <- I_quants$W.triplet_sum
  I_triplet_sum <- as.vector(I_triplet_sum)
  I_n_neighbours <- I_quants$n.neighbours
  I_begfin <- I_quants$W.begfin

  Istar <- diag(apply(I, 1, sum)) - I
  I_Q <- 0 * Istar + diag(rep(1-0,K))


  # 5.2 Create the determinant ----------------------------------------------

  # if(!fix_rho){
  #   Wstar_eigen <- eigen(Wstar)
  #   Wstar_val <- Wstar_eigen$values
  #   det_Q <- sum(log((rho * Wstar_val + (1-rho))))
  # }else{}


  # 5.3 Check for islands ---------------------------------------------------

  W_list<- mat2listw(W, style = "B")
  W_nb <- W_list$neighbours
  W_islands <- n.comp.nb(W_nb)
  islands <- W_islands$comp.id
  islands_all <- rep(islands,J)
  n_islands <- max(W_islands$nc)


  # 5.4 Make matrix versions of the variables -------------------------------

  Y_vec <- as.numeric(t(Y))
  trials_vec <- as.numeric(t(trials))


  # 6 Run Bayesian Model ----------------------------------------------------
  # 6.1 Start Timer ---------------------------------------------------------

  if(verbose){
    cat("Generating",
        n_keep,
        "post burnin and thinned (if requested) samples_\n",
        sep = " ")
    progressBar <- txtProgressBar(style = 3)
    percentage_points<-round((1:100/100)*n_sample)
  }else{
    percentage_points<-round((1:100/100)*n_sample)
  }


  # 6.2 Create MCMC Samples -------------------------------------------------
  j <- 1
  for(j in 1:n_sample){

    # 6.2.1 Sample from Y - data augmentation ---------------------------------

    if(n_miss>0){
      Y_DA[miss_locator] <- rbinom(n=n_miss,
                                   size=trials[miss_locator],
                                   prob=prob[miss_locator])
      failures_DA <- trials - Y_DA
    }else{}


    # 6.2.2 Sample from beta --------------------------------------------------

    offset_temp <- phi + theta + offset
    for(r in 1:J){
      temp_beta <- binomialbetaupdateRW(X = X_standardised, nsites = K, p = p,
                                        beta = beta[ ,r], offset = offset_temp[ ,r],
                                        y = Y_DA[ ,r], failures = failures_DA[ ,r],
                                        prior_meanbeta = prior_beta_mean,
                                        prior_varbeta = prior_beta_var,
                                        nblock = n_beta_block,
                                        beta_tune = proposal_sd_beta[r],
                                        block_list = list_block)
      beta[ ,r] <- temp_beta[[1]]
      accept_beta[r] <- accept_beta[r] + temp_beta[[2]]
      accept_beta[(r+J)] <- accept_beta[(r+J)] + n_beta_block
    }
    regression <- X_standardised %*% beta


    # 6.2.3 Sample from phi ---------------------------------------------------

    den_offset <- 1 * W_triplet_sum + 1 - 1
    phi_offset <- regression + theta + offset
    temp_phi <- binomialmcarupdateGaussianRW(Wtriplet = W_triplet,
                                             Wbegfin = W_begfin,
                                             nsites = K,
                                             nvar = J,
                                             phi = phi,
                                             Y = Y_DA,
                                             failures = failures_DA,
                                             phioffset = phi_offset,
                                             denoffset = den_offset,
                                             tau2 =  tau2,
                                             alpha = tau2_AB,
                                             rho = 1,
                                             phi_tune = proposal_sd_phi)
    phi <- temp_phi[[1]]
    for(r in 1:J){
      phi[ ,r] <- phi[ ,r] - mean(phi[ ,r])
    }
    accept_phi[1] <- accept_phi[1] + temp_phi[[2]]
    accept_phi[2] <- accept_phi[2] + K


    # 6.2.4 Sample from theta --------------------------------------------------

    den_offset <- 0 * W_triplet_sum + 1 - 0
    theta_offset <- regression + phi + offset
    temp_theta <- binomialmcarupdateGaussianRW(Wtriplet = W_triplet,
                                               Wbegfin = W_begfin,
                                               nsites = K,
                                               nvar = J,
                                               phi = theta,
                                               Y = Y_DA,
                                               failures = failures_DA,
                                               phioffset = theta_offset,
                                               denoffset = den_offset,
                                               tau2 =  sigma2,
                                               alpha = sigma2_AB,
                                               rho = 0,
                                               phi_tune = proposal_sd_theta)
    theta <- temp_theta[[1]]
    for(r in 1:J){
      theta[ ,r] <- theta[ ,r] - mean(theta[ ,r])
    }
    accept_theta[1] <- accept_theta[1] + temp_theta[[2]]
    accept_theta[2] <- accept_theta[2] + K


    # 6.2.4 Sample from tau2 --------------------------------------------------

    den_offset <- 1 * W_triplet_sum + 1 - 1
    temp_tau2 <- Gaussiantau2updateRW(Wtriplet = W_triplet,
                                      Wbegfin = W_begfin,
                                      nsites = K,
                                      nvar = J,
                                      phi = phi,
                                      denoffset = den_offset,
                                      rho = 1,
                                      prior_tau2 = prior_tau2,
                                      tau2_tune = proposal_sd_tau2,
                                      tau2 = tau2,
                                      alpha = tau2_AB)
    tau2 <- temp_tau2[[1]]
    accept_tau2[1] <- accept_tau2[1] + temp_tau2[[2]]
    accept_tau2[2] <- accept_tau2[2] + 1


    # 6.2.4 Sample from sigma2 --------------------------------------------------

    den_offset <- 0 * W_triplet_sum + 1 - 0
    temp_sigma2 <- Gaussiantau2updateRW(Wtriplet = W_triplet,
                                        Wbegfin = W_begfin,
                                        nsites = K,
                                        nvar = J,
                                        phi = theta,
                                        denoffset = den_offset,
                                        rho = 0,
                                        prior_tau2 = prior_sigma2,
                                        tau2_tune = proposal_sd_sigma2,
                                        tau2 = sigma2,
                                        alpha = sigma2_AB)
    sigma2 <- temp_sigma2[[1]]
    accept_sigma2[1] <- accept_sigma2[1] + temp_sigma2[[2]]
    accept_sigma2[2] <- accept_sigma2[2] + 1

    # 6.2.5 Sample from rho ---------------------------------------------------

    # if(!fix_rho)
    # {
    #   temp3 <- ClaytonrhoupdateRW(Wtriplet = W_triplet,
    #                               Wtripletsum = W_triplet_sum,
    #                               Wbegfin = W_begfin,
    #                               nsites = K, nvar = J,
    #                               phi = phi,
    #                               denoffset = den_offset,
    #                               rho = rho,
    #                               rho_tune = proposal_sd_rho,
    #                               tau2 = tau2, alpha = alpha,
    #                               Inverse = Copula_Inv)
    #
    #   rho <- temp3[[1]]
    #   # Q <- temp3[[2]]
    #   # det_Q <- temp3[[3]]
    #   accept_rho[1] <- accept_rho[1] + temp3[[2]]
    #   accept_rho[2] <- accept_rho[2] + 1
    # }else
    # {}


    # 6.2.6 Sample from tau2_AB  -------------------------------------------------

    den_offset <- 1 * W_triplet_sum + 1 - 1
    if(!fix_tau2_AB)
    {
      temp_tau2_AB <- GaussianalphaupdateRW(Wtriplet = W_triplet,
                                            Wbegfin = W_begfin,
                                            nsites = K, nvar = J,
                                            phi = phi,
                                            denoffset = den_offset,
                                            rho = 1,
                                            alpha_tune = proposal_sd_tau2_AB,
                                            tau2 = tau2,
                                            alpha = tau2_AB)

      tau2_AB <- temp_tau2_AB[[1]]
      accept_tau2_AB[1] <- accept_tau2_AB[1] + temp_tau2_AB[[2]]
      accept_tau2_AB[2] <- accept_tau2_AB[2] + 1
    }else
    {}

    # 6.2.7 Sample from sigma2_AB -------------------------------------------------

    den_offset <- 0 * W_triplet_sum + 1 - 0
    if(!fix_sigma2_AB)
    {
      temp_sigma2_AB <- GaussianalphaupdateRW(Wtriplet = W_triplet,
                                              Wbegfin = W_begfin,
                                              nsites = K, nvar = J,
                                              phi = theta,
                                              denoffset = den_offset,
                                              rho = 0,
                                              alpha_tune = proposal_sd_sigma2_AB,
                                              tau2 = sigma2,
                                              alpha = sigma2_AB)

      sigma2_AB <- temp_sigma2_AB[[1]]
      accept_sigma2_AB[1] <- accept_sigma2_AB[1] + temp_sigma2_AB[[2]]
      accept_sigma2_AB[2] <- accept_sigma2_AB[2] + 1
    }else
    {}

    # 6.3 Calculate the deviance ----------------------------------------------

    lp <- regression + phi + offset
    prob <- exp(lp)  / (1 + exp(lp))
    fitted <- trials * prob
    loglike <- dbinom(x=as.numeric(t(Y)),
                      size=as.numeric(t(trials)),
                      prob=as.numeric(t(prob)),
                      log=TRUE)


    # 6.4 Save the results ----------------------------------------------------

    if(j > burnin & (j-burnin)%%thin==0){
      ele <- (j - burnin) / thin
      samples_beta[ele, ] <- as.numeric(beta)
      samples_phi[ele, ] <- as.numeric(t(phi))
      samples_theta[ele, ] <- as.numeric(t(theta))

      samples_tau2[ele, ] <- as.numeric(tau2)
      samples_sigma2[ele, ] <- as.numeric(sigma2)
      if(!fix_tau2_AB){
        samples_tau2_AB[ele, ] <- as.numeric(tau2_AB)
      }
      if(!fix_sigma2_AB){
        samples_sigma2_AB[ele, ] <- as.numeric(sigma2_AB)
      }
      # samples_Tau[ele, , ] <- Tau
      # samples_Sigma[ele, , ] <- Sigma

      # if(!fix_rho) samples_rho[ele, ] <- rho
      # if(!fix_alpha) samples_alpha[ele, ] <- alpha

      samples_loglike[ele, ] <- loglike
      samples_fitted[ele, ] <- as.numeric(t(fitted))
      if(n_miss>0) samples_Y[ele, ] <- Y_DA[miss_locator]
    }else{}


    # 6.5 Self tune the acceptance probabilities ------------------------------

    if(ceiling(j/100)==floor(j/100) & j < burnin)
    {
      #### Update the proposal sds
      for(r in 1:J)
      {
        if(p>2)
        {
          proposal_sd_beta[r] <- common_acceptrates1(accept_beta[c(r, (r+J))],
                                                     proposal_sd_beta[r],
                                                     40, 50)
        }else
        {
          proposal_sd_beta[r] <- common_acceptrates1(accept_beta[c(r, (r+J))],
                                                     proposal_sd_beta[r],
                                                     30, 40)
        }
      }

      proposal_sd_phi <- common_acceptrates1(accept_phi,
                                             proposal_sd_phi,
                                             40, 50)

      proposal_sd_theta <- common_acceptrates1(accept_theta,
                                               proposal_sd_theta,
                                               40, 50)

      proposal_sd_tau2 <- common_acceptrates1(accept_tau2,
                                              proposal_sd_tau2,
                                              40, 50)
      proposal_sd_sigma2 <- common_acceptrates1(accept_sigma2,
                                                proposal_sd_sigma2,
                                                40, 50)

      proposal_sd_tau2_AB <- common_acceptrates1(accept_tau2_AB,
                                                 proposal_sd_tau2_AB,
                                                 20, 30)
      proposal_sd_sigma2_AB <- common_acceptrates1(accept_sigma2_AB,
                                                   proposal_sd_sigma2_AB,
                                                   20, 30)

      # if(!fix_rho)
      # {
      #   proposal_sd_rho <- common_acceptrates2(accept_rho[1:2],
      #                                          proposal_sd_rho,
      #                                          40, 50, 0.5)
      # }

      # if(!fix_alpha){
      #   proposal_sd_alpha <- common_acceptrates1(accept_alpha,
      #                                            proposal_sd_alpha,
      #                                            40, 50)
      # }

      accept_beta <- rep(0,2*J)

      accept_phi <- rep(0,2)
      accept_theta <- rep(0,2)

      accept_tau2 <- rep(0,2)
      accept_sigma2 <- rep(0,2)

      accept_tau2_AB <- rep(0,2)
      accept_sigma2_AB <- rep(0,2)

      # accept_rho <- rep(0,2)
      # accept_alpha <- rep(0,2)

    }else
    {}

    # 6.6 Print progress to the console ---------------------------------------

    if(j %in% percentage_points & verbose)
    {
      setTxtProgressBar(progressBar, j/n_sample)
    }
  }

  # End Timer ---------------------------------------------------------------

  ##### end timer
  if(verbose)
  {
    cat("\nSummarising results.")
    close(progressBar)
  }else
  {}


  # 7 Summaries and save the results ----------------------------------------

  # 7.1 Compute the acceptance rates ----------------------------------------

  accept_beta <- 100 * sum(accept_beta[1:J]) / sum(accept_beta[(J+1):(2*J)])
  accept_phi <- 100 * accept_phi[1] / accept_phi[2]
  accept_theta <- 100 * accept_theta[1] / accept_theta[2]

  accept_tau2 <- 100*accept_tau2[1] / accept_tau2[2]
  accept_sigma2 <- 100*accept_sigma2[1] / accept_sigma2[2]

  # # accept_Tau <- 100
  # accept_Sigma <- 100

  if(!fix_tau2_AB){
    accept_tau2_AB <- 100 * accept_tau2_AB[1] / accept_tau2_AB[2]
  }else{
    accept_tau2_AB <- NA
  }

  if(!fix_sigma2_AB){
    accept_sigma2_AB <- 100*accept_sigma2_AB[1] / accept_sigma2_AB[2]
  }else
  {
    accept_sigma2_AB <- NA
  }

  accept_final <- c(accept_beta,
                    accept_phi,
                    accept_theta,
                    accept_tau2,
                    accept_sigma2,
                    accept_tau2_AB,
                    accept_sigma2_AB)
  names(accept_final) <- c("beta", "phi", "theta", "tau2","sigma2", "tau2_AB", "sigma2_AB")


  # 7.2 Compute the fitted deviance -----------------------------------------

  mean_beta <- matrix(apply(samples_beta, 2, mean),
                      nrow=p, ncol=J, byrow=F)
  mean_phi <- matrix(apply(samples_phi, 2, mean),
                     nrow=K, ncol=J, byrow=T)
  mean_theta <- matrix(apply(samples_theta, 2, mean),
                       nrow=K, ncol=J, byrow=T)

  mean_logit <- X_standardised %*% mean_beta + mean_phi + mean_theta + offset
  mean_prob <- exp(mean_logit)  / (1 + exp(mean_logit))
  fitted_mean <- trials * mean_prob
  deviance_fitted <- -2 * sum(dbinom(x=as.numeric(t(Y)),
                                     size=as.numeric(t(trials)),
                                     prob=as.numeric(t(mean_prob)),
                                     log=TRUE),
                              na.rm=TRUE)


  # 7.3 Model fit criteria --------------------------------------------------

  modelfit <- common_modelfit(samples_loglike, deviance_fitted)


  # 7.4 transform the parameters back to the original covariate scale -------

  samples_beta_orig <- samples_beta
  for(r in 1:J){
    samples_beta_orig[ ,((r-1)*p+1):(r*p)] <- common_betatransform(
      samples_beta[ ,((r-1)*p+1):(r*p)],
      X_indicator,
      X_mean,
      X_sd,
      p, FALSE)
  }


  # 7.5 Create a summary object ---------------------------------------------

  samples_beta_orig <- mcmc(samples_beta_orig)
  summary_beta <- t(apply(samples_beta_orig,
                          2,
                          quantile,
                          c(0.5, 0.025, 0.975)))
  summary_beta <- cbind(summary_beta,
                        rep(n_keep, p),
                        rep(accept_beta,J*p),
                        effectiveSize(samples_beta_orig),
                        geweke.diag(samples_beta_orig)$z)
  col_name <- rep(NA, p*(J-1))

  if(is.null(colnames(Y))){
    for(r in 1:J){
      col_name[((r-1)*p+1):(r*p)] <- paste("Variable ",
                                           r,
                                           " - ",
                                           colnames(X), sep="")
    }
  }else{
    for(r in 1:J){
      col_name[((r-1)*p+1):(r*p)] <- paste(colnames(Y)[r],
                                           " - ",
                                           colnames(X), sep="")
    }
  }
  rownames(summary_beta) <- col_name
  colnames(summary_beta) <- c("Median", "2.5%", "97.5%",
                              "n_sample", "%accept", "n_effective",
                              "geweke.diag")

  summary_hyper <- array(NA, c(2*(J+1) ,7))
  for (r in 1:J) {
    summary_hyper[r, 1:3] <- quantile(samples_tau2[,r], c(0.5, 0.025, 0.975))
    summary_hyper[r, 4:7] <- c(n_keep,
                               accept_tau2,
                               effectiveSize(samples_tau2[,r]),
                               geweke.diag(samples_tau2[,r])$z)
    summary_hyper[J+1+r, 1:3] <- quantile(samples_sigma2[,r], c(0.5, 0.025, 0.975))
    summary_hyper[J+1+r, 4:7] <- c(n_keep,
                                   accept_sigma2,
                                   effectiveSize(samples_sigma2[,r]),
                                   geweke.diag(samples_sigma2[,r])$z)
  }

  if(!fix_tau2_AB)
  {
    summary_hyper[(J+1), 1:3] <- quantile(samples_tau2_AB, c(0.5, 0.025, 0.975))
    summary_hyper[(J+1), 4:7] <- c(n_keep,
                                   accept_tau2_AB,
                                   effectiveSize(samples_tau2_AB),
                                   geweke.diag(samples_tau2_AB)$z)
  }else
  {
    summary_hyper[(J+1), 1:3] <- c(tau2_AB, tau2_AB, tau2_AB)
    summary_hyper[(J+1), 4:7] <- rep(NA, 4)
  }

  if(!fix_sigma2_AB)
  {
    summary_hyper[(J+2+J), 1:3] <- quantile(samples_sigma2_AB, c(0.5, 0.025, 0.975))
    summary_hyper[(J+2+J), 4:7] <- c(n_keep,
                                     accept_sigma2_AB,
                                     effectiveSize(samples_sigma2_AB),
                                     geweke.diag(samples_sigma2_AB)$z)
  }else
  {
    summary_hyper[(J+2+J), 1:3] <- c(sigma2_AB, sigma2_AB, sigma2_AB)
    summary_hyper[(J+2+J), 4:7] <- rep(NA, 4)
  }

  # summary_hyper[1:J, 1] <- diag(apply(samples_Sigma, c(2,3), quantile, c(0.5)))
  # summary_hyper[1:J, 2] <- diag(apply(samples_Sigma, c(2,3), quantile, c(0.025)))
  # summary_hyper[1:J, 3] <- diag(apply(samples_Sigma, c(2,3), quantile, c(0.975)))
  #
  # summary_hyper[J+1, 1] <- apply(samples_Sigma, c(2,3), quantile, c(0.5))[1,2]
  # summary_hyper[J+1, 2] <- apply(samples_Sigma, c(2,3), quantile, c(0.025))[1,2]
  # summary_hyper[J+1, 3] <- apply(samples_Sigma, c(2,3), quantile, c(0.975))[1,2]
  #
  # summary_hyper[1:(J+1), 4] <- n_keep
  # summary_hyper[1:(J+1), 5] <- accept_Sigma
  # summary_hyper[1:J, 6] <- diag(apply(samples_Sigma, c(2,3), effectiveSize))
  # summary_hyper[J+1, 6] <- apply(samples_Sigma, c(2,3), effectiveSize)[1,2]
  #
  # for(r in 1:J){
  #   summary_hyper[r, 7] <- geweke.diag(samples_Sigma[ ,r,r])$z
  # }
  # summary_hyper[J+1, 7] <- geweke.diag(samples_Sigma[ ,1,2])$z

  # if(!fix_rho)
  # {
  #   summary_hyper[(J+1), 1:3] <- quantile(samples_rho, c(0.5, 0.025, 0.975))
  #   summary_hyper[(J+1), 4:7] <- c(n_keep, accept_rho, effectiveSize(samples_rho), geweke.diag(samples_rho)$z)
  # }else
  # {
  #   summary_hyper[(J+1), 1:3] <- c(rho, rho, rho)
  #   summary_hyper[(J+1), 4:7] <- rep(NA, 4)
  # }
  #

  summary_results <- rbind(summary_beta, summary_hyper)

  rownames(summary_results)[((J*p)+1): nrow(summary_results)] <- c(paste(rep("tau2",J), 1:J, sep = "-"),"tau2-12", paste(rep("sigma2",J), 1:J, sep = "-"),"sigma2-12")

  summary_results[ , 1:3] <- round(summary_results[ , 1:3], 4)
  summary_results[ , 4:7] <- round(summary_results[ , 4:7], 1)


  # 7.6 Create the fitted values and residuals ------------------------------

  fitted_values <- matrix(apply(samples_fitted, 2, mean),
                          nrow=K, ncol=J,
                          byrow=T)
  response_residuals <- Y - fitted_values
  pearson_residuals <- response_residuals /sqrt(fitted_values * (1 - mean_prob))
  residuals <- list(response=response_residuals,
                    pearson=pearson_residuals)


  # 7.7 Compile and return the results --------------------------------------

  model_string <- c("Likelihood model - Binomial (logit link function)",
                    "\nRandom effects model - Gaussian Copula BYM MCAR\n")
  if(fix_tau2_AB) samples_tau2_AB=NA
  if(fix_sigma2_AB) samples_sigma2_AB=NA
  if(n_miss==0) samples_Y = NA
  samples <- list(beta=samples_beta_orig,
                  phi=mcmc(samples_phi),
                  theta=mcmc(samples_theta),
                  tau2=samples_tau2,
                  sigma2=samples_sigma2,
                  tau2_AB=mcmc(samples_tau2_AB),
                  sigma2_AB=mcmc(samples_sigma2_AB),
                  fitted=mcmc(samples_fitted),
                  Y=mcmc(samples_Y),
                  loglike = mcmc(samples_loglike))
  results <- list(summary_results=summary_results,
                  samples=samples,
                  fitted_values=fitted_values,
                  residuals=residuals,
                  modelfit=modelfit,
                  accept=accept_final,
                  localised_structure=NULL,
                  formula=formula,
                  model=model_string,
                  X=X)


  # 8 Finish by stating the time taken --------------------------------------

  if(verbose)
  {
    b<-proc.time()
    cat("Finished in ", round(b[3]-a[3], 1), "seconds.\n")
    results <- c(results, elapsed = round(b[3]-a[3], 1))
  }else
  {}

  class(results) <- "CARBayes"

  return(results)
}


# End of this section -----------------------------------------------------
