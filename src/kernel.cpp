#include "kernel.hpp"

kernelMatrix::kernelMatrix(){
  numRows = 0;
  numCols = 0;
  kernelData = NULL;
}
kernelMatrix::kernelMatrix(int input_NumRows,int input_NumCols,double(*input_Kernel)(int i, int j, void* kernelData),void* input_KernelData){
  numRows    = input_NumRows;
  numCols    = input_NumCols;
  kernel     = input_Kernel;
  kernelData = input_KernelData;
}

kernelMatrix::~kernelMatrix(){

}

Eigen::MatrixXd kernelMatrix::block(const int min_i,const int min_j,const int blk_NumRows,const int blk_NumCols) const{
  Eigen::MatrixXd result(blk_NumRows,blk_NumCols);
  for (int j = 0; j < blk_NumCols; j++)
    for (int i = 0; i < blk_NumRows; i++)
      result(i,j) = kernel(i + min_i,j + min_j,kernelData);
  return result;
}

Eigen::MatrixXd kernelMatrix::operator*(const Eigen::MatrixXd & X) const{
  assert(X.rows() == numRows);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(numRows,X.cols());
  for (int j = 0; j < X.cols(); j++)
    for (int i = 0; i < numRows; i++)
      for (int k = 0; k < numCols; k++)
	result(i,j) += (*this)(i,k) * X(k,j);
  return result;
}

double kernelMatrix::operator()(const int nRow, const int nCol) const{
  return kernel(nRow,nCol,kernelData);
}

int kernelMatrix::rows() const{
  return numRows;
}

int kernelMatrix::cols() const{
  return numCols;
}
