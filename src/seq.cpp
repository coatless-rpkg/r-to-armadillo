#include <r2arma.h>

//' @title Generate a sequence of values
//' @description Creates a vector containing a sequence of values starting at the initial point and going to the terminal point.
//' @param a An \code{int}, that denotes the starting point.
//' @param b An \code{int}, that denotes the ending point.
//' @return A \code{vector} containing values moving from a to b. There are no restrictions on A's range.
//' @seealso \code{\link{rwishart}} 
//' @author James J Balamuta
//' @examples 
//' # Call with the following data:
//' seq_int(3, 5)
//' seq_int(5, 3)
// [[Rcpp::export]]
arma::vec seq_int(int a, int b){
  int d = abs(b-a)+1;
  
  int inc = ( a < b ? 1 : -1 );
  arma::vec s(d);
  
  for(unsigned int i = 0; i < d; i++){
    s(i) = i*inc + a;
  }
  
  return s;
}