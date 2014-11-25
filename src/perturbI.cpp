#include "perturbI.hpp"

perturbI::perturbI(){

  topU  = NULL;
  topV  = NULL;
  bottU = NULL;
  bottV = NULL;
  determinant_       = 1;
  logAbsDeterminant_ = 0;
  eqMatrixStored     = false;
  eqMatrixFactorized = false;
  calculatedDet      = false;
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
  eqMatrixStored     = true;
  eqMatrixFactorized = false;
  calculatedDet      = false;
  determinant_       = 1;
  logAbsDeterminant_ = 0;

};

void perturbI::eqMatrixFactorize(){
  assert(eqMatrixStored == true);
  if (eqMatrixFactorized == false){
    eqLU = Eigen::PartialPivLU<Eigen::MatrixXd>(eqMatrix);
    eqMatrixFactorized = true;
  }
}

Eigen::MatrixXd perturbI::solve(const Eigen::MatrixXd &RHS){
  //Eigen::PartialPivLU<Eigen::MatrixXd> lu(eqMatrix);
  //return (RHS - U * (lu.solve(VT * RHS)));
  eqMatrixFactorize();
  return (RHS - U * (eqLU.solve(VT * RHS)));  
}

void perturbI::calcDeterminant(){
  eqMatrixFactorize();
  if (calculatedDet == false){
    determinant_       = 1;
    logAbsDeterminant_ = 0;
    Eigen::MatrixXd luMatrix = eqLU.matrixLU();
    for (int i = 0; i < luMatrix.rows(); i++){
      determinant_ *= luMatrix(i,i);
      logAbsDeterminant_ += log(fabs(luMatrix(i,i)));
    }
    calculatedDet = true;
  }
}

double perturbI::determinant(){
  calcDeterminant();
  return determinant_;
}

double perturbI::logAbsDeterminant(){
  calcDeterminant();
  return logAbsDeterminant_;
}
