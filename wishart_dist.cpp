#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' @title Generate a sequence of values
//' @description Creates a vector containing a sequence of values starting at the initial point and going to the terminal point.
//' @param a An \code{int}, that denotes the starting point.
//' @param b An \code{int}, that denotes the ending point.
//' @return A \code{vector} containing values moving from a to b. There are no restrictions on A's range.
//' @seealso \code{\link{rwishart}} 
//' @author James J Balamuta
//' @examples 
//' #Call with the following data:
//' seq_cpp(3, 5)
//' seq_cpp(5, 3)
// [[Rcpp::export]]
arma::vec seq_cpp(int a, int b){
  int d = abs(b-a)+1;
  
  int inc = ( a < b ? 1 : -1 );
  arma::vec s(d);
  
  for(unsigned int i = 0; i < d; i++){
    s(i) = i*inc + a;
  }
  
  return s;
}


//' @title Generate Random Wishart Distribution
//' @description Creates a random wishart distribution when given degrees of freedom and a sigma matrix. 
//' @usage rwishart(df, sigma)
//' @param df An \code{int}, which gives the degrees of freedom of the Wishart.  (> 0)
//' @param sigma A \code{matrix} with dimensions df x df that provides the covariance matrix. 
//' @return A \code{matrix} that is a Wishart distribution, aka the sample covariance matrix of a Multivariate Normal Distribution
//' @seealso \code{\link{riwishart}} and \code{\link{probitHLM}}
//' @author James J Balamuta
//' @examples 
//' #Call with the following data:
//' rwishart(3, diag(1))
// [[Rcpp::export]]
arma::mat rwishart(unsigned int df, const arma::mat& sigma) {
  // Random Numbers from R! 
  RNGScope scope;
  // Need to know Sigma Matrix dimension...
  unsigned int p = sigma.n_rows;
  // Take the choleski decomp
  arma::mat CC = chol(sigma);
  
  /* Say df = 9, and p = 6
  * then 9, 8, 7, 6, 5, 4 should appear!
  */
  arma::vec dfs = seq_cpp(df,df-p+1);

  // We'll use this vector to implant chisq values on diagonal
  unsigned int dfs_size = dfs.n_elem;
  arma::vec chi_vals(dfs_size);
  for(unsigned int i = 0; i < dfs_size; i++){
    // Obtain a random chisquare value for each degree of freedom.
    chi_vals(i) = sqrt( R::rchisq(dfs(i))  );
  }
  
  //Z matrix
  arma::mat Z(p,p);
  
  if(p > 1){
    //Fill matrix with random normal numbers
    Z.imbue( norm_rand );
    //Take an upper triangle matrix (includes diag)
    Z = arma::trimatu(Z);
  }
  
  //Replace diagonal with chisq values.
  Z.diag() = chi_vals;
  
  return (Z*CC).t() * (Z*CC);
}

//' @title Generate Random Inverse Wishart Distribution
//' @description Creates a random inverse wishart distribution when given degrees of freedom and a sigma matrix. 
//' @usage riwishart(df, sigma)
//' @param df An \code{int} that represents the degrees of freedom.  (> 0)
//' @param sigma A \code{matrix} with dimensions df x df that provides the covariance matrix. 
//' @return A \code{matrix} that is an inverse wishart distribution.
//' @seealso \code{\link{rwishart}}
//' @author James J Balamuta
//' @examples 
//' #Call with the following data:
//' riwishart(3, diag(2))
// [[Rcpp::export]]
arma::mat riwishart(unsigned int df, const arma::mat& S) {
  return rwishart(df,S.i()).i();
}