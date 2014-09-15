#ifndef KERNEL_HPP
#define KERNEL_HPP

#include "Eigen/Dense"

class kernelMatrix{

public:

  kernelMatrix();
  kernelMatrix(int input_NumRows,int input_NumCols,double (*input_Kernel)(int i, int j, void* kernelData),void* input_KernelData);
  ~kernelMatrix();
  Eigen::MatrixXd block(const int min_i,const int min_j,const int numRows,const int numCols) const;
  Eigen::MatrixXd operator*(const Eigen::MatrixXd & X) const;
  double operator()(const int nRow, const int nCol) const;
  int rows() const;
  int cols() const;





private:
  int numRows;
  int numCols;
  void* kernelData;
  double (*kernel)(int i, int j , void* kernelData);





};

#endif
