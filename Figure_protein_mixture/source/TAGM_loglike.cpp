#include <Rcpp.h>
using namespace Rcpp;

//------------------------------------------------
// helper function for printing a single value or series of values
template<typename TYPE>
void print(TYPE x) {
  Rcpp::Rcout << x << "\n";
}

template<typename TYPE, typename... Args>
void print(TYPE first, Args... args) {
  Rcpp::Rcout << first << " ";
  print(args...);
}

//------------------------------------------------
// helper function for printing contents of a vector or set
template<class TYPE>
void print_vector(const TYPE &x) {
  for (auto i : x) {
    Rcpp::Rcout << i << " ";
  }
  Rcpp::Rcout << "\n";
}

//------------------------------------------------
// helper function for printing contents of a matrix
template<class TYPE>
void print_matrix(const std::vector<std::vector<TYPE>> &x) {
  for (int i = 0; i < x.size(); ++i) {
    print_vector(x[i]);
  }
  Rcpp::Rcout << "\n";
}

//------------------------------------------------
// helper function for printing contents of a 3D array
template<class TYPE>
void print_array(const std::vector<std::vector<std::vector<TYPE>>> &x) {
  for (int i = 0; i < x.size(); ++i) {
    Rcpp::Rcout << "--- slice " << i+1 << " ---\n";
    print_matrix(x[i]);
  }
  Rcpp::Rcout << "\n";
}

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
std::vector<std::vector<double>> rcpp_to_matrix_double(Rcpp::List x) {
  int nrow = int(x.size());
  std::vector< std::vector<double> > x_mat(nrow);
  for (int i = 0; i < nrow; ++i) {
    x_mat[i] = Rcpp::as<std::vector<double>>(x[i]);
  }
  return x_mat;
}

//------------------------------------------------
// multivariate gamma function of x with dimension p, returned in log space
double lmvgamma_func(double x, double p) {
  double ret = 0.25*p*(p - 1)*log(M_PI);
  for (int i = 0; i < p; ++i) {
    ret += lgamma(x - 0.5*i);
  }
  return ret;
}

//------------------------------------------------
// calculate Cholesky decomposition of positive definite matrix sigma
void cholesky(std::vector<std::vector<double>> &chol, const std::vector<std::vector<double>> &sigma) {
  
  for (int i = 0; i < int(sigma.size()); ++i) {
    for (int j = 0; j < (i+1); ++j) {
      chol[i][j] = sigma[i][j];
      if (i == j) {
        if (i > 0) {
          for (int k = 0; k < i; ++k) {
            chol[i][i] -= chol[i][k]*chol[i][k];
          }
        }
        chol[i][i] = sqrt(chol[i][i]);
      } else {
        if (j > 0) {
          for (int k = 0; k < j; ++k) {
            chol[i][j] -= chol[i][k]*chol[j][k];
          }
        }
        chol[i][j] /= chol[j][j];
      }
    }
  }
  
}

//------------------------------------------------
// return (in log space) determinant of positive definite matrix. chol is the
// Cholesky decomposition of the target matrix.
double log_determinant(const std::vector<std::vector<double>> &chol) {
  double ret = 0;
  for (int i = 0; i < int(chol.size()); i++) {
    ret += 2*log(chol[i][i]);
  }
  return ret;
}

//------------------------------------------------
// return diagonal matrix of dimension d with elements x
std::vector<std::vector<double>> diag(int d, double x) {
  std::vector<std::vector<double>> ret(d, std::vector<double>(d));
  for (int i = 0; i < d; i++) {
    ret[i][i] = x;
  }
  return ret;
}

//------------------------------------------------
// find the inverse of the matrix M by Gauss-Jordan elimination. Note that M is
// modified inside the function, and therefore cannot be passed in by reference.
std::vector<std::vector<double>> inverse(std::vector<std::vector<double>> M) {
  
  int d = int(M.size());
  std::vector<std::vector<double>> IM = diag(d, 1.0);
  
  double temp;
  for (int j = 0; j < d; j++) {
    for (int i = j; i < d; i++) {
      if (i == j) {
        temp = M[i][i];
        for (int k = 0; k < d; k++) {
          M[i][k] /= temp;
          IM[i][k] /= temp;
        }
      } else {
        if (M[i][j]!=0) {
          temp = M[i][j];
          for (int k = 0; k < d; k++) {
            M[i][k] -= temp * M[j][k];
            IM[i][k] -= temp * IM[j][k];
          }
        }
      }
    }
  }
  for (int j = (d - 1); j > 0; j--) {
    for (int i = (j - 1); i >= 0; i--) {
      temp = M[i][j];
      for (int k = 0; k < d; k++) {
        M[i][k] -= temp * M[j][k];
        IM[i][k] -= temp * IM[j][k];
      }
    }
  }
  
  return IM;
}

//------------------------------------------------
// density of multivariate normal distribution given mean vector mu, double
// logdet, and matrix chol_inverse. logdet is the logarithm of the determinant
// of the covariance matrix, chol_inverse is the inverse of the Cholesky
// decomposition of the covariance matrix.
double dmnorm1(const std::vector<double> &x,
               const std::vector<double> &mu,
               double logdet,
               const std::vector< std::vector<double> > &chol_inverse) {
  
  int d = int(x.size());
  double ret = -0.5*d*log(2*M_PI) - 0.5*logdet;
  double tmp;
  for (int i = 0; i < d; i++) {
    tmp = 0;
    for (int j = 0; j < (i + 1); j++) {
      tmp += (x[j] - mu[j]) * chol_inverse[i][j];
    }
    ret += -0.5*tmp*tmp;
  }
  
  return ret;
}

//------------------------------------------------
// equivalent to dmnorm1, but computes determinants etc. internally rather than
// as inputs. More convenient in terms of inputs, but less efficient if the same
// input matrices will be used a large number of times.
double dmnorm2(const std::vector<double> &x,
               const std::vector<double> &mu,
               const std::vector<std::vector<double>> &sigma) {
  
  // calculate determinants etc.
  int d = int(x.size());
  std::vector<std::vector<double>> sigma_chol(d, std::vector<double>(d));
  cholesky(sigma_chol, sigma);
  double logdet = log_determinant(sigma_chol);
  std::vector<std::vector<double>> sigma_chol_inverse = inverse(sigma_chol);
  
  // run dmnorm1 function
  double ret = dmnorm1(x, mu, logdet, sigma_chol_inverse);
  
  return ret;
}

//------------------------------------------------
// density of inverse Wishart distribution on matrix sigma given scale matrix
// psi and degrees of freedom nu. Sigma and psi are input pre-transformed to
// save time when running this function many times with the same inputs.
// sigma_inv is the inverse of sigma. sigma_chol and psi_chol are the Cholesky
// decompositions of sigma and psi, respectively.
double dinvwish1(const std::vector<std::vector<double>> &sigma_inv,
                 const std::vector<std::vector<double>> &sigma_chol,
                 const std::vector<std::vector<double>> &psi,
                 const std::vector<std::vector<double>> &psi_chol,
                 double nu) {
  
  // get basic properties
  int d = sigma_inv.size();
  
  // get matrix determinants
  double sigma_logdet = log_determinant(sigma_chol);
  double psi_logdet = log_determinant(psi_chol);
  
  // calculate probability density
  double ret = 0.5*nu*psi_logdet - 0.5*nu*d*log(2.0) - lmvgamma_func(0.5*nu, d) - 0.5*(nu + d + 1)*sigma_logdet;
  for (int i = 0; i < d; i++) {
    for (int j = 0; j < d; j++) {
      ret -= 0.5 * psi[i][j] * sigma_inv[i][j];
    }
  }
  
  return ret;
}

//------------------------------------------------
// equivalent to dinvwish1, but computes inverse and Cholesky decomposition
// matrices internally rather than as inputs. More convenient in terms of
// inputs, but less efficient if the same input matrices will be used a large
// number of times.
double dinvwish2(const std::vector<std::vector<double>> &sigma,
                 const std::vector<std::vector<double>> &psi,
                 double nu) {
  
  // get basic properties
  int d = sigma.size();
  
  // get matrix inverses and cholesky decompositions
  std::vector<std::vector<double>> sigma_inv = inverse(sigma);
  std::vector<std::vector<double>> sigma_chol(d, std::vector<double>(d));
  std::vector<std::vector<double>> psi_chol(d, std::vector<double>(d));
  cholesky(sigma_chol, sigma);
  cholesky(psi_chol, psi);
  
  // calculate and return
  return dinvwish1(sigma_inv, sigma_chol, psi, psi_chol, nu);
}

//------------------------------------------------
// [[Rcpp::export]]
SEXP loglike(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc) {
  
  // unpack misc
  int K = misc["K"];
  std::vector<double> t_prob = Rcpp::as< std::vector<double> >(misc["t_prob"]);
  
  // unpack data
  std::vector<double> x = Rcpp::as< std::vector<double> >(data["x"]);
  std::vector<double> y = Rcpp::as< std::vector<double> >(data["y"]);
  int n = x.size();
  
  // unpack parameters
  std::vector<double> mu_x(K);
  std::vector<double> mu_y(K);
  std::vector<double> sigma_x(K);
  std::vector<double> sigma_y(K);
  std::vector<double> rho(K);
  int pi = 0;
  for (int k = 0; k < K; ++k) {
    mu_x[k] = params[pi++];
    mu_y[k] = params[pi++];
    sigma_x[k] = params[pi++];
    sigma_y[k] = params[pi++];
    rho[k] = params[pi++];
  }
  std::vector<double> w_gamma(K);
  double w_gamma_sum = 0.0;
  for (int k = 0; k < K; ++k) {
    w_gamma[k] = params[pi++];
    w_gamma_sum += w_gamma[k];
  }
  std::vector<double> w(K);
  for (int k = 0; k < K; ++k) {
    w[k] = w_gamma[k] / w_gamma_sum;
  }
  double eps = params[pi++];
  
  // sum log-likelihood over all data
  double ret = 0.0;
  std::vector<double> logd(K);
  std::vector<double> const1(K);
  std::vector<double> const2(K);
  std::vector<double> px2(K);
  std::vector<double> pxy(K);
  std::vector<double> py2(K);
  for (int k = 0; k < K; ++k) {
    const1[k] = -log(2*PI*sigma_x[k]*sigma_y[k]) -0.5*log(1 - rho[k]*rho[k]);
    const2[k] = -0.5/double(1.0 - rho[k]*rho[k]);
    px2[k] = 1.0/double(sigma_x[k]*sigma_x[k]);
    pxy[k] = 1.0/double(sigma_x[k]*sigma_y[k]);
    py2[k] = 1.0/double(sigma_y[k]*sigma_y[k]);
  }
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < K; ++k) {
      double dx = x[i] - mu_x[k];
      double dy = y[i] - mu_y[k];
      logd[k] = const1[k] + const2[k] * (dx*dx*px2[k] - 2*rho[k]*dx*dy*pxy[k] + dy*dy*py2[k]);
    }
    double tmp = eps * t_prob[i];
    for (int k = 0; k < K; ++k) {
      tmp += (1 - eps) * w[k] * exp(logd[k]);
    }
    ret += log(tmp);
  }
  
  // return as SEXP
  return Rcpp::wrap(ret);
}

//------------------------------------------------
// [[Rcpp::export]]
SEXP logprior(Rcpp::NumericVector params, Rcpp::List misc) {
  
  // unpack misc
  int K = misc["K"];
  std::vector<double> mu_0 = Rcpp::as<std::vector<double>>(misc["mu_0"]);
  double lambda_0 = Rcpp::as<double>(misc["lambda_0"]);
  double nu_0 = Rcpp::as<double>(misc["nu_0"]);
  std::vector<std::vector<double>> S_0 = rcpp_to_matrix_double(misc["S_0"]);
  double beta = Rcpp::as<double>(misc["beta"]);
  double u = Rcpp::as<double>(misc["u"]);
  double v = Rcpp::as<double>(misc["v"]);
  
  // unpack parameters
  std::vector<double> mu_x(K);
  std::vector<double> mu_y(K);
  std::vector<double> sigma_x(K);
  std::vector<double> sigma_y(K);
  std::vector<double> rho(K);
  int pi = 0;
  for (int k = 0; k < K; ++k) {
    mu_x[k] = params[pi++];
    mu_y[k] = params[pi++];
    sigma_x[k] = params[pi++];
    sigma_y[k] = params[pi++];
    rho[k] = params[pi++];
  }
  std::vector<double> w_gamma(K);
  for (int k = 0; k < K; ++k) {
    w_gamma[k] = params[pi++];
  }
  double eps = params[pi++];
  
  // get mean and covariance matrices into convenient form
  std::vector<std::vector<double>> mu(K, std::vector<double>(2));
  std::vector<std::vector<std::vector<double>>> sigma_mat(K, std::vector<std::vector<double>>(2, std::vector<double>(2)));
  std::vector<std::vector<std::vector<double>>> sigma_mat_scaled(K, std::vector<std::vector<double>>(2, std::vector<double>(2)));
  for (int k = 0; k < K; ++k) {
    mu[k][0] = mu_x[k];
    mu[k][1] = mu_y[k];
    
    sigma_mat[k][0][0] = sigma_x[k];
    sigma_mat[k][1][1] = sigma_y[k];
    sigma_mat[k][0][1] = sigma_mat[k][1][0] = sigma_x[k]*sigma_y[k]*rho[k];
    
    sigma_mat_scaled[k][0][0] = sigma_mat[k][0][0] / lambda_0;
    sigma_mat_scaled[k][1][1] = sigma_mat[k][1][1] / lambda_0;
    sigma_mat_scaled[k][0][1] = sigma_mat_scaled[k][1][0] = sigma_mat[k][0][1] / lambda_0;
  }
  
  // apply prior on component means and covariances
  double ret = 0.0;
  for (int k = 0; k < K; ++k) {
    ret += dinvwish2(sigma_mat[k], S_0, nu_0);
    ret += dmnorm2(mu[k], mu_0, sigma_mat_scaled[k]);
  }
  
  // apply prior on weights
  for (int k = 0; k < K; ++k) {
    ret += R::dgamma(w_gamma[k], beta, 1.0, true);
  }
  
  // apply prior on error
  ret += R::dbeta(eps, u, v, true);
  
  // deal with overflow, which can happen in very hot chains
  if (ret < (-DBL_MAX)) {
    ret = -DBL_MAX;
  }
  
  return Rcpp::wrap(ret);
}

//------------------------------------------------
// [[Rcpp::export]]  
SEXP create_xptr(std::string function_name) {  
  typedef SEXP (*funcPtr_likelihood)(Rcpp::NumericVector params, Rcpp::List data, Rcpp::List misc);  
  typedef SEXP (*funcPtr_prior)(Rcpp::NumericVector params, Rcpp::List misc);  
  
  if (function_name == "loglike"){
    return(Rcpp::XPtr<funcPtr_likelihood>(new funcPtr_likelihood(&loglike)));
  } 
  if (function_name == "logprior"){
    return(Rcpp::XPtr<funcPtr_prior>(new funcPtr_prior(&logprior)));
  } 
  
  stop("cpp function %i not found", function_name);
}
