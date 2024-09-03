#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix matrixMultiply(NumericMatrix mat1, NumericMatrix mat2) {
  // Check dimensions
  if (mat1.ncol() != mat2.nrow()) {
    stop("Incompatible matrix dimensions for multiplication.");
  }

  int n = mat1.nrow(), k = mat2.ncol(), m = mat1.ncol();
  NumericMatrix result(n, k);

  // Perform matrix multiplication
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      double sum = 0;
      for (int l = 0; l < m; ++l) {
        sum += mat1(i, l) * mat2(l, j);
      }
      result(i, j) = sum;
    }
  }

  return result;
}

// [[Rcpp::export]]
NumericVector getDiagonal(NumericMatrix mat) {
  int n = std::min(mat.nrow(), mat.ncol()); // The length of the diagonal is the minimum of the number of rows and columns
  NumericVector diagonal(n);

  for (int i = 0; i < n; ++i) {
    diagonal[i] = mat(i, i); // Elements where row index equals column index
  }

  return diagonal;
}


// [[Rcpp::export]]
NumericVector linpredcompute(NumericMatrix X, const int nsites, const int p,
                             NumericVector beta, NumericVector offset){
  //Create new objects
  // Compute the linear predictor
  NumericVector linpred(nsites);
  double temp;


  //  Compute the linear predictor via a double for loop
  for(int j = 0; j < nsites; j++)
  {
    temp = 0;

    for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];

    linpred[j] = temp + offset[j];
  }


  // Return the result
  return linpred;
}

// [[Rcpp::export]]
List binomialbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta,
                          NumericVector offset, NumericVector y,  NumericVector failures,
                          NumericVector prior_meanbeta, NumericVector prior_varbeta,
                          const int nblock, double beta_tune, List block_list){
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);

  // Create two beta vectors
  NumericVector beta_old(p);
  NumericVector beta_new(p);
  for(int g=0; g<p; g++)
  {
    beta_old[g] = beta[g];
    beta_new[g] = beta[g];
  }

  // Update each block in turn
  for(int r=0; r<nblock; r++)
  {
    // Determine the block to update
    IntegerVector idx = block_list[r];
    int len = block_list[(nblock+r)];

    // Propose a value
    for(int g=0; g<len; g++)
    {
      beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }


    // Compute the acceptance ratio - full conditionals
    oldlikebit = 0;
    newlikebit=0;
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);
    for(int j = 0; j < nsites; j++)
    {
      p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
      p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
      oldlikebit = oldlikebit +  y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
      newlikebit = newlikebit +  y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
    }
    likebit = newlikebit - oldlikebit;

    priorbit = 0;
    for(int g = 0; g < len; g++)
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }


    // Accept or reject hte proposal
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance)
    {
      for(int g=0; g<len; g++)
      {
        beta_old[idx[g]] = beta_new[idx[g]];
      }
      accept = accept + 1;
    }
    else
    {
      for(int g=0; g<len; g++)
      {
        beta_new[idx[g]] = beta_old[idx[g]];
      }
    }
  }


  // Compute the acceptance probability and return the value
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;
}

// [[Rcpp::export]]
List binomialmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                          const int nsites,  const int nvar, NumericMatrix phi,
                          NumericMatrix Y, NumericMatrix failures,
                          NumericMatrix phioffset, NumericVector denoffset,
                          NumericMatrix Sigmainv, double rho, double phi_tune){
  // Update the spatially correlated random effects
  //Create new objects
  NumericMatrix fcprec(nvar, nvar);
  int rowstart=0, rowend=0, accept=0;
  NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
  NumericVector diffcurrent(nvar), diffprop(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar);
  NumericVector pold(nvar), pnew(nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;

  //  Update each random effect in turn
  for(int j = 0; j < nsites; j++)
  {
    // Calculate the prior precision and mean
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = denoffset[j] * Sigmainv(_,r);
    }

    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
    fcmean = rho * sumphi / denoffset[j];


    // Generate the proposal distribution mean and propose a value
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = rnorm(1, phi(j,r), phi_tune)[0];
    }


    // Compute the prior ratio
    diffcurrent = phi(j,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));
      quadprop[r] = sum(diffprop * fcprec(_,r));
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);

    // Likelihood ratio
    pold = exp(phioffset(j,_) + phi(j,_)) / (1 + exp(phioffset(j,_) + phi(j,_)));
    pnew = exp(phioffset(j,_) + propphi) / (1 + exp(phioffset(j,_) + propphi));
    oldlikebit = sum(Y(j,_) * log(pold) + failures(j,_) * log(1 - pold));
    newlikebit = sum(Y(j,_) * log(pnew) + failures(j,_) * log(1 - pnew));


    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance)
    {
      phi(j,_) = propphi;
      accept = accept + 1;
    }
    else
    {
    }
  }


  // Return the results
  List out(2);
  out[0] = phi;
  out[1] = accept;
  return out;
}

// [[Rcpp::export]]
List mcarrhoupdateRW(Rcpp::NumericMatrix Wstar,
                     Rcpp::NumericVector Wstar_val,
                     int K,
                     int J,
                     double rho,
                     Rcpp::NumericMatrix Q,
                     double det_Q,
                     double proposal_sd_rho,
                     Rcpp::NumericMatrix phi,
                     Rcpp::NumericMatrix Sigma_inv) {
  // Update the strength factor of spatial effect - rho
  //Create new objects
  double a=0, b=1, proposal_rho;
  NumericMatrix Q_prop(K, K);
  double det_Q_prop=0;
  double logprob_current, logprob_proposal, hastings, acceptance;
  int accept=0;
  // Sampling from a truncated normal distribution
  proposal_rho = r_truncnorm(rho, proposal_sd_rho, a, b);

  // Compute Q.prop and its determinant component
  for (int i = 0; i < K; ++i) {
    for (int j = 0; j < K; ++j) {
      Q_prop(i, j) = proposal_rho * Wstar(i, j) + (i == j) * (1 - proposal_rho);
    }
    det_Q_prop += std::log(proposal_rho * Wstar_val[i] + (1 - proposal_rho));
  }

  // Computing the log-probability of the current and proposed values
  NumericMatrix ele_sum_current = matrixMultiply(matrixMultiply(matrixMultiply(transpose(phi), Q), phi), Sigma_inv);
  NumericMatrix ele_sum_prop = matrixMultiply(matrixMultiply(matrixMultiply(transpose(phi), Q_prop), phi), Sigma_inv);

  logprob_current = 0.5 * J * det_Q - 0.5 * sum(getDiagonal(ele_sum_current));
  logprob_proposal = 0.5 * J * det_Q_prop - 0.5 * sum(getDiagonal(ele_sum_prop));

  // Compute the Hastings ratio
  hastings = std::log(d_truncnorm(rho, proposal_rho, proposal_sd_rho, a, b,0)) - std::log(d_truncnorm(proposal_rho, rho, proposal_sd_rho, a, b));

  // Compute the acceptance probability
  acceptance = std::exp(logprob_proposal - logprob_current + hastings);

  if (acceptance > R::runif(0, 1)) {
    rho = proposal_rho;
    det_Q = det_Q_prop;
    Q = Q_prop;
    accept = accept + 1;
  }else
  {
  }
  // Return the results
  List out(4);
  out[0] = rho;
  out[1] = Q;
  out[2] = det_Q;
  out[3] = accept;
  return out;
}


// [[Rcpp::export]]
List gaussianmcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                          const int nsites,  const int nvar, NumericMatrix phi,
                          NumericMatrix phioffset, NumericVector denoffset,
                          NumericMatrix Sigmainv, double rho, NumericVector nu2,
                          double phi_tune){
  // Update the spatially correlated random effects
  //Create new objects
  NumericMatrix fcprec(nvar, nvar);
  int rowstart=0, rowend=0, accept=0;
  NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
  NumericVector diffcurrent(nvar), diffprop(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar);
  NumericVector lpold(nvar), lpnew(nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;

  //  Update each random effect in turn
  for(int j = 0; j < nsites; j++)
  {
    // Calculate the prior precision and mean
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = denoffset[j] * Sigmainv(_,r);
    }

    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
    fcmean = rho * sumphi / denoffset[j];

    // Generate a proposal value
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = rnorm(1, phi(j,r), phi_tune)[0];
    }


    // Compute the prior ratio
    diffcurrent = phi(j,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));
      quadprop[r] = sum(diffprop * fcprec(_,r));
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);

    // Likelihood ratio
    lpold = pow((phioffset(j,_) - phi(j,_)),2);
    lpnew = pow((phioffset(j,_) - propphi),2);
    oldlikebit = 0.5 * sum(lpold / nu2);
    newlikebit = 0.5 * sum(lpnew / nu2);

    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit + oldlikebit - newlikebit);
    if(runif(1)[0] <= acceptance)
    {
      phi(j,_) = propphi;
      accept = accept + 1;
    }
    else
    {
    }
  }

  List out(2);
  out[0] = phi;
  out[1] = accept;
  return out;
}

