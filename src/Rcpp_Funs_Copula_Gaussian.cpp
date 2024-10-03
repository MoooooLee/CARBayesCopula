#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>

#include <cmath>  // Include for std::atanh
using namespace Rcpp;

// This file contains the following functions:
//[[Rcpp::depends(RcppArmadillo, RcppDist)]]

// [[Rcpp::export]]
double GaussianDensity(NumericVector CDF_q, double alpha, double nvar) {
  double logden;

  double x = CDF_q[0];
  double y = CDF_q[1];

  double term1 = -0.5 * log(1-alpha*alpha);
  double term2_numerator = 2 * alpha * x * y - alpha * alpha * (x * x + y * y);
  double term2_dominator = 2 * (1 - alpha * alpha);

  logden = term1 + term2_numerator / term2_dominator;

  return logden;
}


// [[Rcpp::export]]
List binomialmcarupdateGaussianRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                                  const int nsites, const int nvar,
                                  NumericMatrix phi, NumericMatrix Y, NumericMatrix failures,
                                  NumericMatrix phioffset, NumericVector denoffset,
                                  NumericVector tau2,
                                  double alpha, double rho, double phi_tune){
  // Update the spatially correlated random effects

  // Create new objects
  NumericVector fcsd(nvar);
  int rowstart=0, rowend=0, accept=0;
  NumericVector sumphi(nvar), fcmean(nvar), propphi(nvar);
  NumericVector diffcurrent(nvar), diffprop(nvar);
  NumericVector LCDF_current(nvar), CDF_current(nvar), LPDF_current(nvar), CDF_current_q(nvar);
  NumericVector LCDF_prop(nvar), CDF_prop(nvar), LPDF_prop(nvar), CDF_prop_q(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar);
  NumericVector pold(nvar), pnew(nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
  // double prod_CDF_current_q, prod_CDF_prop_q;

  // Update each random effect in turn
  for(int j = 0; j<nsites; j++)
  {
    for(int r = 0; r<nvar; r++){
      fcsd[r] = sqrt(tau2[r]/denoffset[j]);
    }

    rowstart = Wbegfin(j,0)-1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0, nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      sumphi += Wtriplet(l,2) * phi((Wtriplet(l,1) - 1),_);
    }
    fcmean = rho * sumphi / denoffset[j];

    // Generate the proposal distribution mean and purpose a value
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = rnorm(1, phi(j,r), phi_tune)[0];
    }

    // Compute the prior ratio
    diffcurrent = phi(j,_) - fcmean;
    diffprop = propphi - fcmean;

    for(int r=0; r<nvar; r++)
    {
      // LCDF_current[r] = R::pnorm(phi(j,r), fcmean[r], fcsd[r], true, true);
      // CDF_current[r] = R::pnorm(phi(j,r), fcmean[r], fcsd[r], true, false);
      LPDF_current[r] = R::dnorm(phi(j,r), fcmean[r], fcsd[r], true);
      CDF_current_q[r] = (phi(j,r) - fcmean[r]) / fcsd[r];

      // LCDF_prop[r] = R::pnorm(propphi[r], fcmean[r], fcsd[r], true, true);
      // CDF_prop[r] = R::pnorm(propphi[r], fcmean[r], fcsd[r], true, false);
      LPDF_prop[r] = R::dnorm(propphi[r], fcmean[r], fcsd[r], true);
      CDF_prop_q[r] = (propphi[r] - fcmean[r]) / fcsd[r];
    }

    // oldpriorbit = ((-alpha - 1) * sum(LCDF_current) + (- 1 / alpha - 2) * log(sum(pow(CDF_current,- alpha)) - 1) + sum(LPDF_current));
    // newpriorbit = ((-alpha - 1) * sum(LCDF_prop) + (- 1 / alpha - 2) * log(sum(pow(CDF_prop, - alpha)) - 1) + sum(LPDF_prop));

    // prod_CDF_current_q = 1;
    // prod_CDF_prop_q = 1;
    // for(int r=0; r<nvar; r++)
    // {
    //   prod_CDF_current_q *= CDF_current_q[r];
    //   prod_CDF_prop_q *= CDF_prop_q[r];
    // }


    // oldpriorbit = -1/(2*(1-pow(alpha,2))) * ((pow(alpha,2) * sum(pow(CDF_current_q,2))) - 2*alpha*prod_CDF_current_q) + sum(LPDF_current) ;
    // newpriorbit = -1/(2*(1-pow(alpha,2))) * ((pow(alpha,2) * sum(pow(CDF_prop_q,2))) - 2*alpha*prod_CDF_prop_q) + sum(LPDF_prop) ;

    oldpriorbit = GaussianDensity(CDF_current_q, alpha, nvar) + sum(LPDF_current);
    newpriorbit = GaussianDensity(CDF_prop_q, alpha, nvar) + sum(LPDF_prop);

    // Likelihood ratio
    pold = exp(phioffset(j,_) + phi(j,_)) / (1 + exp(phioffset(j,_) + phi(j,_)));
    pnew = exp(phioffset(j,_) + propphi) / (1 + exp(phioffset(j,_) + propphi));
    oldlikebit = sum(Y(j,_) * log(pold) + failures(j,_) * log(1 - pold));
    newlikebit = sum(Y(j,_) * log(pnew) + failures(j,_) * log(1 - pnew));

    // Accept or reject the value
    acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
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
List Gaussiantau2updateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                          const int nsites, const int nvar,
                          NumericMatrix phi,
                          NumericVector denoffset,
                          double rho, NumericVector prior_tau2,
                          double tau2_tune, NumericVector tau2, double alpha){
  // Create new objects
  NumericVector fcsd(nvar), fcsd_prop(nvar);
  int rowstart=0, rowend=0, accept=0;
  NumericVector sumphi(nvar), fcmean(nvar), proptau2(nvar);
  double a = prior_tau2[0];
  double b = prior_tau2[1];
  NumericVector LCDF_current(nvar), CDF_current(nvar), LPDF_current(nvar), CDF_current_q(nvar);
  NumericVector LCDF_prop(nvar), CDF_prop(nvar), LPDF_prop(nvar), CDF_prop_q(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar);
  NumericVector pold(nvar), pnew(nvar);
  NumericVector oldlikebit(nsites), newlikebit(nsites);
  NumericVector oldhasting(nvar), newhasting(nvar);
  double oldpriorbit, newpriorbit, hasting, acceptance;
  // double prod_CDF_current_q, prod_CDF_prop_q;

  //Update each random effect in turn
  //Generate the proposal distribution mean and purpose a value
  for(int r=0; r<nvar; r++)
  {
    proptau2[r] = rtruncnorm(1, tau2[r], tau2_tune, 0, 1000)[0];
  }

  for(int j=0; j<nsites; j++)
  {
    for(int r=0; r<nvar; r++)
    {
      fcsd(r) = sqrt(tau2(r) / denoffset(j));
      fcsd_prop(r) = sqrt(proptau2(r) / denoffset(j));
    }
    rowstart = Wbegfin(j,0)-1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0, nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
    }
    fcmean = rho * sumphi / denoffset[j];

    //Likelihood Ratio
    for(int r=0; r<nvar; r++)
    {
      // LCDF_current[r] = R::pnorm(phi(j,r), fcmean[r], fcsd[r], true, true);
      // CDF_current[r] = R::pnorm(phi(j,r), fcmean[r], fcsd[r], true, false);
      LPDF_current[r] = R::dnorm(phi(j,r), fcmean[r], fcsd[r], true);
      CDF_current_q[r] = (phi(j,r) - fcmean[r]) / fcsd[r];


      // LCDF_prop[r] = R::pnorm(phi(j,r), fcmean[r], fcsd_prop[r], true, true);
      // CDF_prop[r] = R::pnorm(phi(j,r), fcmean[r], fcsd_prop[r], true, false);
      LPDF_prop[r] = R::dnorm(phi(j,r), fcmean[r], fcsd_prop[r], true);
      CDF_prop_q[r] = (phi(j,r) - fcmean[r]) / fcsd_prop[r];

    }
    // oldlikebit[j] = ((-alpha - 1) * sum(LCDF_current) + (- 1 / alpha - 2) * log(sum(pow(CDF_current,- alpha)) - 1) + sum(LPDF_current));
    // newlikebit[j] = ((-alpha - 1) * sum(LCDF_prop) + (- 1 / alpha - 2) * log(sum(pow(CDF_prop, - alpha)) - 1) + sum(LPDF_prop));

    // prod_CDF_current_q = 1;
    // prod_CDF_prop_q = 1;
    // for(int r=0; r<nvar; r++)
    // {
    //   prod_CDF_current_q *= CDF_current_q[r];
    //   prod_CDF_prop_q *= CDF_prop_q[r];
    // }

    oldlikebit[j] = GaussianDensity(CDF_current_q, alpha, nvar) + sum(LPDF_current);
    newlikebit[j] = GaussianDensity(CDF_prop_q, alpha, nvar) + sum(LPDF_prop);
    // oldlikebit[j] = (-1/(2*(1-pow(alpha,2)))) * ((pow(alpha,2) * sum(pow(CDF_current_q,2))) - 2 * alpha * prod_CDF_current_q) + sum(LPDF_current) ;
    // newlikebit[j] = (-1/(2*(1-pow(alpha,2)))) * ((pow(alpha,2) * sum(pow(CDF_prop_q,2))) - 2 * alpha * prod_CDF_prop_q) + sum(LPDF_prop) ;

  }

  //Prior Ratio
  oldpriorbit = sum(-(a+1)*log(tau2) - b/tau2);
  newpriorbit = sum(-(a+1)*log(proptau2) - b/proptau2);

  //Hasting
  for(int r = 0; r < nvar; r++)
  {
    oldhasting[r] = d_truncnorm(tau2[r], proptau2[r], tau2_tune, 0, 1000, true);
    newhasting[r] = d_truncnorm(proptau2[r],tau2[r], tau2_tune, 0, 1000, true);
  }
  hasting = sum(oldhasting-newhasting);

  // Accept or reject the value
  acceptance = exp(newpriorbit - oldpriorbit + sum(newlikebit - oldlikebit) + hasting);

  if(runif(1)[0] <= acceptance)
  {
    tau2= proptau2;
    accept = accept + 1;
  }
  else
  {
  }

  //Return the results
  List out(2);
  out[0] = tau2;
  out[1] = accept;
  return out;
}

// [[Rcpp::export]]
List GaussianrhoupdateRW(NumericMatrix Wtriplet, NumericVector Wtripletsum,
                         NumericMatrix Wbegfin,
                         const int nsites, const int nvar,
                         NumericMatrix phi,
                         NumericVector denoffset,
                         double rho, double rho_tune,
                         NumericVector tau2, double alpha){
  // Create new objects
  NumericVector sumphi(nvar), fcsd(nvar), fcsd_prop(nvar), fcmean(nvar), fcmean_prop(nvar);
  int rowstart=0, rowend=0, accept=0;
  NumericVector LCDF_current(nvar), CDF_current(nvar), LPDF_current(nvar), CDF_current_q(nvar);
  NumericVector LCDF_prop(nvar), CDF_prop(nvar), LPDF_prop(nvar), CDF_prop_q(nvar);
  NumericVector oldlikebit(nsites), newlikebit(nsites), propdenoffset(nsites);
  double oldhasting, newhasting, proprho;
  double hasting, acceptance;
  // double prod_CDF_current_q, prod_CDF_prop_q;
  //double oldpriorbit, newpriorbit;

  //Update each random effect in turn
  //Generate the proposal distribution mean and purpose a value
  proprho = r_truncnorm(rho, rho_tune, 0, 1);
  for(int j =0; j<nsites; j++)
  {
    propdenoffset[j] = proprho * Wtripletsum[j] + 1 - proprho;
  }

  for(int j = 0; j < nsites; j++)
  {
    for(int r = 0; r < nvar; r++)
    {
      fcsd[r] = sqrt(tau2[r] / denoffset[j]);
      fcsd_prop[r] = sqrt(tau2[r] / propdenoffset[j]);
    }
    rowstart = Wbegfin(j,0)-1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0, nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
    }
    fcmean = rho * sumphi / denoffset[j];
    fcmean_prop = proprho * sumphi / propdenoffset[j];

    //Likelihood Ratio
    for(int r=0; r<nvar; r++)
    {
      // LCDF_current[r] = R::pnorm(phi(j,r), fcmean[r], fcsd[r], true, true);
      // CDF_current[r] = R::pnorm(phi(j,r), fcmean[r], fcsd[r], true, false);
      LPDF_current[r] = R::dnorm(phi(j,r), fcmean[r], fcsd[r], true);
      CDF_current_q[r] = (phi(j,r) - fcmean[r]) / fcsd[r];

      // LCDF_prop[r] = R::pnorm(phi(j,r), fcmean_prop[r], fcsd_prop[r], true, true);
      // CDF_prop[r] = R::pnorm(phi(j,r), fcmean_prop[r], fcsd_prop[r], true, false);
      LPDF_prop[r] = R::dnorm(phi(j,r), fcmean_prop[r], fcsd_prop[r], true);
      CDF_prop_q[r] = (phi(j,r) - fcmean_prop[r]) / fcsd_prop[r];
    }
    // oldlikebit[j] = ((-alpha - 1) * sum(LCDF_current) + (- 1 / alpha - 2) * log(sum(pow(CDF_current,- alpha)) - 1) + sum(LPDF_current));
    // newlikebit[j] = ((-alpha - 1) * sum(LCDF_prop) + (- 1 / alpha - 2) * log(sum(pow(CDF_prop, - alpha)) - 1) + sum(LPDF_prop));

    // prod_CDF_current_q = 1;
    // prod_CDF_prop_q = 1;
    // for(int r=0; r<nvar; r++)
    // {
    //   prod_CDF_current_q *= CDF_current_q[r];
    //   prod_CDF_prop_q *= CDF_prop_q[r];
    // }


    oldlikebit[j] = GaussianDensity(CDF_current_q, alpha, nvar) + sum(LPDF_current);
    newlikebit[j] = GaussianDensity(CDF_prop_q, alpha, nvar) + sum(LPDF_prop);
    // oldlikebit[j] = (-1/(2*(1-pow(alpha,2)))) * ((pow(alpha,2) * sum(pow(CDF_current_q,2))) - 2 * alpha * prod_CDF_current_q) + sum(LPDF_current) ;
    // newlikebit[j] = (-1/(2*(1-pow(alpha,2)))) * ((pow(alpha,2) * sum(pow(CDF_prop_q,2))) - 2 * alpha * prod_CDF_prop_q) + sum(LPDF_prop) ;
  }
  //Prior Ratio - Assuming unif

  //Hasting
  oldhasting = d_truncnorm(rho, proprho, rho_tune, 0, 1, 1);
  newhasting = d_truncnorm(proprho, rho, rho_tune, 0, 1, 1);
  hasting = oldhasting-newhasting;

  // Accept or reject the value
  acceptance = exp(sum(newlikebit - oldlikebit) + hasting);

  if(runif(1)[0] <= acceptance)
  {
    rho= proprho;
    accept = accept + 1;
  }
  else
  {
  }

  //Return the results
  List out(2);
  out[0] = rho;
  out[1] = accept;
  return out;
}

// [[Rcpp::export]]
List GaussianalphaupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin,
                           const int nsites, const int nvar,
                           NumericMatrix phi,
                           NumericVector denoffset,
                           double rho, double alpha_tune,
                           NumericVector tau2, double alpha){

  // Create new objects
  NumericVector fcsd(nvar);
  int rowstart=0, rowend=0, accept=0;
  NumericVector sumphi(nvar), fcmean(nvar);
  NumericVector LCDF_current(nvar), CDF_current(nvar), LPDF_current(nvar), CDF_current_q(nvar);
  NumericVector oldlikebit(nsites), newlikebit(nsites);
  double oldhasting, newhasting;
  double propalpha, hasting, acceptance;
  // double oldpriorbit, newpriorbit;
  // double prod_CDF_current_q;

  //Update each random effect in turn
  //Generate the proposal distribution mean and purpose a value
  double trans_alpha = std::atanh(alpha);
  double prop_trans_alpha = rnorm(1, trans_alpha, alpha_tune)[0];
  propalpha = tanh(prop_trans_alpha);

  // propalpha = r_truncnorm(alpha, alpha_tune, -1, 1);
  for(int j=0; j<nsites; j++)
  {
    for(int r=0; r<nvar; r++)
    {
      fcsd(r) = sqrt(tau2(r) / denoffset(j));
    }
    rowstart = Wbegfin(j,0)-1;
    rowend = Wbegfin(j,1);
    sumphi = rep(0, nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      sumphi += Wtriplet(l, 2) * phi((Wtriplet(l,1) - 1),_);
    }
    fcmean = rho * sumphi / denoffset[j];

    //Likelihood Ratio
    for(int r=0; r<nvar; r++)
    {
      // LCDF_current[r] = R::pnorm(phi(j,r), fcmean[r], fcsd[r], true, true);
      // CDF_current[r] = R::pnorm(phi(j,r), fcmean[r], fcsd[r], true, false);
      // LPDF_current[r] = R::dnorm(phi(j,r), fcmean[r], fcsd[r], true);
      CDF_current_q[r] = (phi(j,r) - fcmean[r]) / fcsd[r];
    }

    // oldlikebit[j] = ((-alpha - 1) * sum(LCDF_current) + (- 1 / alpha - 2) * log(sum(pow(CDF_current,- alpha)) - 1) + sum(LPDF_current));
    // newlikebit[j] = ((-propalpha - 1) * sum(LCDF_current) + (- 1 / propalpha - 2) * log(sum(pow(CDF_current, - propalpha)) - 1) + sum(LPDF_current));

    // prod_CDF_current_q = 1;
    // for(int r=0; r<nvar; r++)
    // {
    //   prod_CDF_current_q *= CDF_current_q[r];
    // }

    oldlikebit[j] = GaussianDensity(CDF_current_q, alpha, nvar);
    newlikebit[j] = GaussianDensity(CDF_current_q, propalpha, nvar);
    // oldlikebit[j] = -1/(2*(1-pow(alpha,2))) * ((pow(alpha,2) * sum(pow(CDF_current_q,2))) - 2*alpha*prod_CDF_current_q)-1/2*log(1-pow(alpha,2));
    // newlikebit[j] = -1/(2*(1-pow(propalpha,2))) * ((pow(propalpha,2) * sum(pow(CDF_current_q,2))) - 2*propalpha*prod_CDF_current_q)-1/2*log(1-pow(propalpha,2)) ;
    //

  }
  //Prior Ratio - Assume unif prior

  //Hasting
  oldhasting = log(1-propalpha*propalpha);
  newhasting = log(1-alpha*alpha);
  hasting = oldhasting - newhasting;
  // hasting = 0;
  // Accept or reject the value
  acceptance = exp(sum(newlikebit - oldlikebit) + hasting);

  if(runif(1)[0] <= acceptance)
  {
    alpha= propalpha;
    accept = accept + 1;
  }
  else
  {
  }

  //Return the results
  List out(2);
  out[0] = alpha;
  out[1] = accept;
  return out;
}



