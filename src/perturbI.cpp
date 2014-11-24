#include "perturbI.hpp"

perturbI::perturbI(){

  topU  = NULL;
  topV  = NULL;
  bottU = NULL;
  bottV = NULL;
  determinant_ = 1;
  logAbsDeterminant_ = 0;
  eqMatrixStored = false;
  eqMatrixFactorized = false;
}
perturbI::perturbI(Eigen::MatrixXd* topU_,Eigen::MatrixXd* topV_, Eigen::MatrixXd* bottU_, Eigen::MatrixXd* bottV_){

  topU = topU_;
  topV = topV_;
  bottU = bottU_;
  bottV = bottV_;

  assert(topU->rows() + bottU->rows() == topV->rows() + bottV->rows());
  assert(topU->cols()  == topV->cols());
  assert(bottU->cols() == bottV->cols());

  int rankTotal = topU->cols() + bottU->cols();
  int blockSize = topU->rows() + bottU->rows();

  U  = Eigen::MatrixXd::Zero(blockSize,rankTotal);
  VT = Eigen::MatrixXd::Zero(rankTotal,blockSize);
  
  U.topLeftCorner(topU->rows(),topU->cols()) = *topU;
  U.bottomRightCorner(bottU->rows(),bottU->cols()) = *bottU;
  VT.topRightCorner(topV->cols(),topV->rows()) = topV->transpose();
  VT.bottomLeftCorner(bottV->cols(),bottV->rows()) = bottV->transpose();
  
  eqMatrix = Eigen::MatrixXd::Identity(rankTotal,rankTotal) + VT * U;
  eqMatrixStored = true;
};

Eigen::MatrixXd eqMatrixFactorize(){


}

Eigen::MatrixXd perturbI::solve(const Eigen::MatrixXd &RHS){
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(eqMatrix);
  return (RHS - U * (lu.solve(VT * RHS)));
}
