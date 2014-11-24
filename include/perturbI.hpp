#ifndef PERTURBI_HPP
#define PERTURBI_HPP

#include <iostream>
#include <Eigen/Dense>

class perturbI{
  bool eqMatrixStored;
  bool eqMatrixFactorized;
  Eigen::MatrixXd* topU;
  Eigen::MatrixXd* topV;
  Eigen::MatrixXd* bottU;
  Eigen::MatrixXd* bottV;
  Eigen::MatrixXd eqMatrix;
  Eigen::MatrixXd U;
  Eigen::MatrixXd VT;
  double determinant_;
  double logAbsDeterminant_;
  Eigen::PartialPivLU<Eigen::MatrixXd> eqLU;

  void eqMatrixFactorize();
public:
  perturbI();
  perturbI(Eigen::MatrixXd* topU_,Eigen::MatrixXd* topV_, Eigen::MatrixXd* bottU_, Eigen::MatrixXd* bottV_);
  Eigen::MatrixXd solve(const Eigen::MatrixXd &RHS);
  double determinant();
  double logAbsDeterminant();
};

#endif
