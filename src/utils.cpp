#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MappedSparseMatrix;

Eigen::MatrixXd createBeta(const Eigen::Map<Eigen::MatrixXd>& G) {
  int n = G.cols();
  Eigen::VectorXd colSum = G.colwise().sum();

  for (int i = 0; i < n; ++i) {
    if (colSum(i) != 0.0) {
      colSum(i) = 1.0 / colSum(i);
    }
  }

  Eigen::MatrixXd betaG = colSum.asDiagonal();
  
  return betaG;
}


// [[Rcpp::export()]]
SEXP gene_GroupTechEffects(const Eigen::Map<Eigen::MatrixXd> & G,
const Eigen::Map<Eigen::MatrixXd> & BG,
const Eigen::Map<Eigen::MatrixXd> & X
)
{
  Eigen::MatrixXd betaG = createBeta(G);
  Eigen::MatrixXd betaBG = createBeta(BG);
  Eigen::MatrixXd Z = X - X * G * betaG * G.adjoint();
  Eigen::MatrixXd W = X - X * BG * betaBG * BG.adjoint();

  Z = Z.cwiseProduct(Z) * G * betaG;
  W = W.cwiseProduct(W) * G * betaG;
  Eigen::MatrixXd GroupTechEffects = Z-W;
  Eigen::MatrixXd OverallTechEffects = GroupTechEffects.rowwise().sum();

  List result;
  result["GroupTechEffects"] = Rcpp::wrap(GroupTechEffects);
  result["OverallTechEffects"] = Rcpp::wrap(OverallTechEffects);
  return result;
}

