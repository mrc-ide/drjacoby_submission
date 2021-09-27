#include <Rcpp.h>
using namespace Rcpp;

// Euclidean distance between two sets of coordiantes (x1, y1), (x2, y2)
double dist(double x1, double y1, double x2, double y2){
  return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  // Distances
  double d1 = dist(params["x3"], params["y3"], 0, 1);
  double d2 = dist(params["x3"], params["y3"], 1, 0);
  double d3 = dist(params["x3"], params["y3"], params["x4"], params["y4"]);
  
  // sum log-likelihood over all data
  // unpack data
  std::vector<double> y = Rcpp::as< std::vector<double> >(data["y"]);
  double ret = 0.0;
  ret += R::dnorm(d1, y[0], 0.1, true);
  ret += R::dnorm(d2, y[1], 0.1, true);
  ret += R::dnorm(d3, y[2], 0.1, true);
    
  // return as SEXP
  return Rcpp::wrap(ret);
}

// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  double ret = 0.0;
  ret += R::dnorm(params["x3"], 0, 10, true);
  ret += R::dnorm(params["y3"], 0, 10, true);
  ret += R::dnorm(params["x4"], 0, 10, true);
  ret += R::dnorm(params["y4"], 0, 10, true);
  return Rcpp::wrap(ret);
}


// NOTE: Do not edit this function name
// [[Rcpp::export]]  
SEXP create_xptr(std::string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  // NOTE: If your loglikelihood function is not called "loglike" please edit:
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  // NOTE: If your logprior function is not called "logprior" please edit:
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 
  
  stop("cpp function %i not found", function_name);
}
