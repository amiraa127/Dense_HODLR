#ifndef LOW_RANK_HPP
#define LOW_RANK_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
#include <assert.h>
#include <iostream>
#include <algorithm>
#include <matrixIO.hpp>
#include <HODLR_Tree.hpp>

double fullPivACA_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int minPivot = 0);

double partialPivACA_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int minPivot = 0);

void PS_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank, const std::string pointChoosingMethod = "Chebyshev",const int minRank = -1);
  
void PS_LowRankApprox_Sp(const Eigen::SparseMatrix<double> & matrixData_Sp,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j,const int numRows, const int numCols, const double tolerance, int &calculatedRank);

void SVD_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank, const int minRank = -1);

int  SVD_LowRankApprox(const Eigen::MatrixXd & matrixData, const double accuracy, Eigen::MatrixXd* Wptr = NULL, Eigen::MatrixXd* Vptr = NULL, Eigen::MatrixXd* Kptr = NULL, int minRank = -1); 

void PS_Boundary_LowRankApprox(const Eigen::MatrixXd & matrixData,const Eigen::SparseMatrix<double> graphData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols,const double tolerance,int & calculatedRank, const int maxDepth = 2,const std::string savePath = "none");  

void extractRowsCols(const Eigen::MatrixXd & matrixData, int min_i,int min_j,int numRows,int numCols,Eigen::MatrixXd & W, Eigen::MatrixXd & K, Eigen::MatrixXd & V,const std::vector<int> & rowIndex,const std::vector<int> & colIndex,const double tolerance, int & calculatedRank,const std::string mode = "fullPivLU");

int  getBoundaryRowColIdx(const Eigen::SparseMatrix<double>  & graphData,const int min_i, const int min_j,const int numRows,const int numCols,const int depth,std::vector<int> & rowIdx,std::vector<int> & colIdx);

int add_LR(Eigen::MatrixXd & result_U,Eigen::MatrixXd & result_K,Eigen::MatrixXd & result_V,const Eigen::MatrixXd & U1, const Eigen::MatrixXd & V1, const Eigen::MatrixXd & U2, const Eigen::MatrixXd & V2,double tol,std::string mode);

int PS_PseudoInverse(Eigen::MatrixXd & colMatrix,Eigen::MatrixXd & rowMatrix, Eigen::MatrixXd & U, Eigen::MatrixXd & V,Eigen::MatrixXd & K,std::vector<int> rowIdxVec,const  double tol,const std::string mode);

#endif
