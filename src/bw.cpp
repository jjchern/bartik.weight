
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

# include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ComputeAlphaBeta(const arma::vec y, const arma::vec x, const arma::mat WW,
                      const arma::mat weight, const arma::mat Z, const arma::vec G) {

    int n = x.n_elem;

    arma::mat I;
    arma::mat M_W = I.eye(n, n) - WW * inv_sympd(WW.t() * weight * WW) * WW.t() * weight;
    arma::vec xx = M_W * x;
    arma::vec yy = M_W * y;
    arma::mat ZZ = M_W * Z;
    arma::vec Alpha = (diagmat(G) * Z.t() * weight * xx) / as_scalar(G.t() * Z.t() * weight * xx);
    arma::vec Beta = (Z.t() * weight * yy) / (Z.t() * weight * xx);

    return List::create(Alpha, Beta);
}
