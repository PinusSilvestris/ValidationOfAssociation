// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

arma::vec swap_ranks(arma::vec R){
    int n = R.size();
    return R.transform( [n](int val) { return (n+1 - val); } );
}

double Wn(int n, double C, double i, double j)
{
    double u = (i+0.5)/(n+1);
    double v = (j+0.5)/(n+1);
    double w = pow(n, 0.5) * (C - u * v) / pow(u * v * (1.0 - u) * (1.0 - v), 0.5);
    return w;
}

arma::vec fill_T(arma::vec &Rx, arma::vec &Ry){
    int n = Rx.size();
    arma::vec T = arma::zeros(n);

    for (int i = 0; i < n; i++)
    {
        T((int)Rx[i] - 1) = Ry[i];
    }
    return T;
}

void calculate_copula_part(arma::vec &T, arma::mat &C){

    int n = T.size();

    for (int i = 1; i <= (double)(n + 1) / 2.0; i++)
    {
        for (int j = 0; j < T[i - 1] && j <= (double)(n + 1) / 2.0; j++)
        {
            C(i, j) = C(i - 1, j);
        }
        for (int j = (int)T[i - 1]; j <= (double)(n + 1) / 2.0; j++)
        {
            C(i, j) = C(i - 1, j) + 1.0/n;
        }
    }
}

void fill_matrix(int n, arma::mat &KS, arma::mat &C, bool x_prim, bool y_prim){

    int sign = pow(-1, x_prim+y_prim);

    for (int i = 0; i <= (double)(n + 1) / 2.0; i++)
    {
        for (int j = 0; j <= (double)(n + 1) / 2.0; j++)
        {
            int x = x_prim ? n-i : i;
            int y = y_prim ? n-j : j;
            KS(x, y) = sign*Wn(n, C(i, j), i, j);
        }
    }}

// [[Rcpp::export]]
arma::mat arma_copula(arma::vec Rx, arma::vec Ry) {
    int n = Rx.size();

    // TODO: decrease memory footprint
    int k = (int)(1.0*n/2.0);

    arma::mat Ctab = arma::zeros(n+1, n+1);
    arma::mat Ctabs22 = arma::zeros(n+1, n+1);
    arma::mat Ctabs12 = arma::zeros(n+1, n+1);
    arma::mat Ctabs21 = arma::zeros(n+1, n+1);
    arma::vec Rsx = arma::zeros(n);
    arma::vec Rsy = arma::zeros(n);
    arma::vec T;
    arma::vec Ts22;
    arma::vec Ts12;
    arma::vec Ts21;

    arma::mat KS = arma::zeros(n+1, n+1);

    Rsx = swap_ranks(Rx);
    Rsy = swap_ranks(Ry);

    T = fill_T(Rx, Ry);
    Ts22 = fill_T(Rsx, Rsy);
    Ts12 = fill_T(Rx, Rsy);
    Ts21 = fill_T(Rsx, Ry);

    // calculation of empirical copula on the grid, for ranks and transformed ranks
    calculate_copula_part(T, Ctab);
    calculate_copula_part(Ts22, Ctabs22);
    calculate_copula_part(Ts12, Ctabs12);
    calculate_copula_part(Ts21, Ctabs21);

    // calculation of 𝑄𝑛
    //∗ (function “Wn” defined at the end of this code)

    fill_matrix(n, KS, Ctab, false, false);
    fill_matrix(n, KS, Ctabs22, true, true);
    fill_matrix(n, KS, Ctabs12, false, true);
    fill_matrix(n, KS, Ctabs21, true, false);

    return KS;
}

// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);

    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector,
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}
