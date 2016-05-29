#ifndef R2ARMA_SEQ_H
#define R2ARMA_SEQ_H
arma::vec seq_int(long int a, long int b);

arma::vec seq_default(long double from, long double to, long unsigned int length_out);
arma::vec seq_default_a(long double a, long double length_out);

arma::vec seq_along_cpp(const arma::vec& along_with);
arma::vec seq_len_cpp(long unsigned int length_out);


#endif
