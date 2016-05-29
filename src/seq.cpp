#include <r2arma.h>

//' @title Generate an integer sequence of values
//' @description Creates a vector containing a sequence of values starting at the initial point and going to the terminal point.
//' @param a A \code{long int} that denotes the starting point.
//' @param b A \code{long int} that denotes the ending point.
//' @return A \code{vector} containing integer values from
//' @seealso \code{\link{rwishart}}
//' @author James J Balamuta
//' @backref src/seq.cpp
//' @backref inst/include/seq.h
//' @examples
//' # Call with the following data:
//' seq_int(3, 5)
//' seq_int(5, 3)
// [[Rcpp::export]]
arma::vec seq_int(long int a, long int b){
  long int d = std::abs(b-a)+1;

  return arma::linspace(a, b, d);
}


//' @title Generate a sequence of values ranging from a to b given a fixed amount of points.
//' @description Creates a vector containing a sequence of values starting at the initial point and going to the terminal point.
//' @param from       A \code{long double} that denotes the starting point.
//' @param to         A \code{long double} indicating the end point.
//' @param length_out A \code{long unsigned int} that denotes the number of points to use.
//' @return A \code{vector} containing values moving from \eqn{[from,to]}.
//' @seealso \code{\link{seq}}
//' @author Anthony R. Colombo, James J Balamuta
//' @backref src/seq.cpp
//' @backref inst/include/seq.h
//' @examples
//' # Call with the following data:
//' seq_default(1, 2, 10)
//' seq_default(2, 1, 10)
// [[Rcpp::export]]
arma::vec seq_default(long double from, long double to, long unsigned int length_out){
  return arma::linspace(from, to, length_out);
}

//' @title Generate a sequence of values ranging from -a to a given a fixed amount of points.
//' @description Creates a vector containing a sequence of values starting at the initial point and going to the terminal point.
//' @param a          A \code{long double} that denotes the starting point.
//' @param length_out A \code{long unsigned int} that denotes the number of points to use.
//' @return A \code{vector} of length \code{length_out} that contains values in \eqn{[-a,a]}.
//' @seealso \code{\link{seq}}
//' @author Anthony R. Colombo, James J Balamuta
//' @examples
//' # Call with the following data:
//' seq_default_a(1, 10)
//' seq_default_a(1, 10)
// [[Rcpp::export]]
arma::vec seq_default_a(long double a, long double length_out){
  return seq_default(-a, a, length_out);
}


//' @title Generate a sequence of values ranging from vector start to vector end
//' @description Creates a vector containing a sequence of values starting at the beginning and going to the end.
//' @param along_with A \code{vec} with length \eqn{n}.
//' @return A \code{vector} of length \eqn{n} that contains integer values in \eqn{[0, n - 1]}.
//' @seealso \code{\link{seq}}
//' @author James J Balamuta
//' @backref src/seq.cpp
//' @backref inst/include/seq.h
//' @examples
//' # Call with the following data:
//' seq_along_cpp(1:10)
//' seq_along_cpp(5:10)
// [[Rcpp::export]]
arma::vec seq_along_cpp(const arma::vec& along_with){
  return seq_len_cpp(along_with.n_elem);
}

//' @title Generate a sequence of values ranging from 0 to n-1
//' @description Creates a vector containing a sequence of values starting at the initial point and going to the terminal point.
//' @param length_out An \code{long unsigned int} that denotes the number of points to use.
//' @return A \code{vector} of length \code{length_out} that contains integer values in \eqn{[0, length_out - 1]}.
//' @seealso \code{\link{seq}}
//' @author James J Balamuta
//' @backref src/seq.cpp
//' @backref inst/include/seq.h
//' @examples
//' # Call with the following data:
//' seq_len_cpp(2)
//' seq_len_cpp(4)
// [[Rcpp::export]]
arma::vec seq_len_cpp(long unsigned int length_out){
  return arma::linspace(0, length_out - 1, length_out);
}
