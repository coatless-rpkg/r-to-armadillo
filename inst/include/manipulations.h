#ifndef R2ARMA_MANIPULATIONS_H
#define R2ARMA_MANIPULATIONS_H

arma::vec get_elements(const arma::mat& x, const arma::uvec& row_ind, const arma::uvec& col_ind);

arma::mat rev_col_subset(arma::mat x, unsigned int start, unsigned int end);

arma::mat rev_row_subset(arma::mat x, unsigned int start, unsigned int end);

arma::vec reverse_vec(arma::vec x);

arma::mat field_to_matrix(arma::field<arma::vec> x);

double sum_field_vec(const arma::field<arma::vec>& x);

#endif
