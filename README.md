
[![Travis-CI Build Status](https://travis-ci.org/coatless/r-to-armadillo.svg?branch=master)](https://travis-ci.org/coatless/r-to-armadillo)

R To Armadillo (`r2armadillo`)
==============================

R is an amazing language to explore and communicate statistics with. However, when we care about attain results quickly, R's primary weakness is shown: speed. The objective of this repository is to strengthen this weakness by providing high performing R base functions within armadillo (and potentially eigen). Futhermore, as the code is written within C++, the code is able to be ported to different platforms at ease. Lastly, an additional benefit of this project is a repository containing programming examples by field.

Feel free to use the functions available within the repository per the MIT License.

Project Details
===============

The functions written here act either as a means to replicate existing R functionality within armadillo or as a helpful interface to armadillo's code. These functions are meant to be easily included within existing code bases with minimal dependencies.

At some point in the future, these functions may end up being merged into the `RcppArmadillo` project.

The current implementation has:

-   [Distributions](https://github.com/coatless/r-to-armadillo/blob/master/src/distributions.cpp)
    -   Random Sampling
        -   Wishart
        -   Inverse Wishart
-   [Time Series](https://github.com/coatless/r-to-armadillo/blob/master/src/ts.cpp)
    -   Difference by lag and number of differences
    -   Discrete Integration by lag, number of differences, and - optionally - initial values.
    -   Convert ARMA process to an infinite MA process
    -   Compute the theoretical autocorrelation function (ACF) for an ARMA process.
    -   Convolution (Moving Average) or Recursive (Autoregression) filters.
    -   Calculate Discrete Fourier Transformation (DFT) for Autocovariance Function (ACF)
-   [Sequences](https://github.com/coatless/r-to-armadillo/blob/master/src/seq.cpp)
    -   Generate Integer Sequence from `a` to `b`
-   [Manipulations](https://github.com/coatless/r-to-armadillo/blob/master/src/manipulations.cpp)
    -   Vector
        -   Reverse vector (e.g. `1:3` =&gt; `3:1`)
    -   Matrix
        -   Reverse subset for row/column (e.g. `0:3` vs. `3:0`)
    -   Field
        -   Convert `field<vec>` to `mat`
        -   Total sum of elements in `field`

Contributing
============

Contributions such as translations of R functions to C++ are welcome via a pull request.

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.
