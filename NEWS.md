# r2arma 0.0.1.9000

Features:

* Added `seq_along_cpp`, `seq_len_cpp`
* Added `seq_default`, and `seq_default_a` (#1, @arcolombo)

Changes:

* Require C++11 compiler for long support and R 3.3.0
* Rewrote `seq_int` to use `arma::linspace`.

Unit Tests:

* Unit testing framework was added to verify outputs of `seq` functions. 

# r2arma 0.0.0.9000

Features:

* Formally placed functions into a package form with easy `LinkTo:` capabilities. 
