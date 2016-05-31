#include <r2arma.h>

//' Subset Non-connected Regions
//'
//' Replicates the subset functionality of a matrix in R.
//' @param x       A \code{matrix} of dimensions M x N
//' @param row_ind A \code{unsigned int vec} that contains the row indices within \eqn{[0,M-1]}.
//' @param col_ind A \code{unsigned int vec} that contains the column indices within \eqn{[0,N-1]}.
//' @return A \code{vec} with each element listed according to index specification.
//' @author JJB
//' @examples
//' # Generate a Matrix
//' m = matrix(1:12, nrow = 4)
//'
//' # Select Non-connect regions
//' row_index = c(1, 2, 1)
//' col_index = c(2, 2, 3)
//'
//' # Subset in R
//' m[cbind(row_index, col_index)]
//'
//' # Subset with Armadillo
//' get_elements(m, row_index - 1, col_index - 1)
// [[Rcpp::export]]
arma::vec get_elements(const arma::mat& x, const arma::uvec& row_ind, const arma::uvec& col_ind){

  unsigned int n = row_ind.n_elem;

  arma::vec out(n);

  for(unsigned int i = 0; i < n; i++){
    out(i) = x(row_ind[i], col_ind[i]);
  }

  return out;
}

//' Reverse Subset Column
//'
//' Subsets the column by going from high indices to low (the reverse of the supported practice)
//' @param x     A \code{matrix} of dimensions M x N
//' @param start A \code{unsigned int} that indicates the starting column.
//' @param end   A \code{unsigned int} that indicates the ending column.
//' @return A \code{matrix} with matrix rows displayed in reverse order
//' @details Consider a vector x=[[1,2],[3,4]].
//' By setting \code{start=1} and \code{end=0}, the function would output x=[[2,1],[4,1]].
//' Start and end must be valid C++ matrix locations. (e.g. matrix cols start at 0 and not 1)
//' @author JJB
//' @examples
//' x = matrix(c(1,2,3,4), nrow = 2,byrow = TRUE)
//' rev_col_subset(x, 1, 0)
// [[Rcpp::export]]
arma::mat rev_col_subset(arma::mat x, unsigned int start, unsigned int end){
  arma::mat A = arma::mat(x.n_rows, start-end+1);
  for(unsigned int i = 0; i < start-end+1; i++){
    A.col(i) = x.col(start-i);
  }
  return A;
}

//' Reverse Subset Row
//'
//' Subsets the row by going from high indices to low (the reverse of the supported practice)
//' @param x A \code{matrix} of dimensions M x N
//' @param start A \code{unsigned int} that indicates the starting row.
//' @param end A \code{unsigned int} that indicates the ending row.
//' @return x A \code{matrix} with matrix rows displayed in reversed order
//' @details Consider a vector x=[[1,2],[3,4]], the function would output x=[[3,4],[1,2]].
//' Start and end must be valid C++ matrix locations. (e.g. matrix rows start at 0 and not 1)
//' @author JJB
//' @examples
//' x = matrix(c(1,2,3,4), nrow=2,byrow=TRUE)
//' rev_row_subset(x, 1, 0)
// [[Rcpp::export]]
arma::mat rev_row_subset(arma::mat x, unsigned int start, unsigned int end){
  arma::mat A = arma::mat(start-end+1, x.n_cols);
  for(unsigned int i = 0; i < start-end+1; i++){
    A.row(i) = x.row(start-i);
  }
  return A;
}

//' Reverse Armadillo Vector
//'
//' Reverses the order of an Armadillo Vector
//' @param x A \code{column vector} of length N
//' @return A \code{column vector} with its contents reversed.
//' @details Consider a vector x=[1,2,3,4,5], the function would output x=[5,4,3,2,1].
//' @author JJB
//' @examples
//' x = 1:5
//' reverse_vec(x)
// [[Rcpp::export]]
arma::vec reverse_vec(arma::vec x) {
   std::reverse(x.begin(), x.end());
   return x;
}

//' Transform an Armadillo field<vec> to a matrix
//'
//' Unlists vectors in a field and places them into a matrix
//' @param x A \code{field<vec>}.
//' @return A \code{mat} containing the field elements within a column.
//' @author JJB
// [[Rcpp::export]]
arma::mat field_to_matrix(arma::field<arma::vec> x){
  unsigned int nx = x.n_elem;
  unsigned int row;
  if(nx > 0){
    row = x(0).n_elem;
  }else{
    row = 999999999;
  }
  arma::mat A(row,nx);
  for(unsigned int i =0; i<nx; i++){
    A.col(i) = x(i);
  }
  return A;
}

//' Accumulation of Armadillo field<vec>
//'
//' Sums vectors in a field into a single variable.
//' @param x A \code{field<vec>}.
//' @return An \code{mat} containing the field elements within a column.
//' @author JJB
// [[Rcpp::export]]
double sum_field_vec(const arma::field<arma::vec>& x){
  unsigned int nelems = x.n_elem;
  double total_elems = 0;

  for(unsigned int i = 0; i < nelems; i++){
    total_elems += sum(x(i));
  }

  return total_elems;
}
