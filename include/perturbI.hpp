#ifndef PERTURBI_HPP
#define PERTURBI_HPP

#include <iostream>
#include <Eigen/Dense>

class perturbI{

  bool eqMatrixStored;
  bool eqMatrixFactorized;
  bool calculatedDet;  

  Eigen::MatrixXd* topU;
  Eigen::MatrixXd* topV;
  Eigen::MatrixXd* bottU;
  Eigen::MatrixXd* bottV;
  Eigen::MatrixXd eqMatrix;
  Eigen::MatrixXd U;
  Eigen::MatrixXd VT;
  Eigen::PartialPivLU<Eigen::MatrixXd> eqLU;
  
  double determinant_;
  double logAbsDeterminant_;


  void eqMatrixFactorize();
  void calcDeterminant();
public:
 
  perturbI();
  perturbI(Eigen::MatrixXd* topU_,Eigen::MatrixXd* topV_, Eigen::MatrixXd* bottU_, Eigen::MatrixXd* bottV_);
  Eigen::MatrixXd solve(const Eigen::MatrixXd &RHS);
  double determinant();
  double logAbsDeterminant();


};

#endif
