#ifndef LOW_RANK_HPP
#define LOW_RANK_HPP


//Standard C++
#include <assert.h>
#include <algorithm>
#include <iostream>
#include <set>

//External Libraries
#include <Eigen/Dense>
#include <Eigen/Sparse>

//Custom Dependancies
#include "helperFunctions.hpp"
#include "HODLR_Tree.hpp"
#include "kernel.hpp"
#include "matrixIO.hpp"

template <typename T>
double fullPivACA_LowRankApprox(const T & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int minPivot = 0);
 
template <typename T>
double partialPivACA_LowRankApprox(const T & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int maxRank = -1,const int minPivot = 0);

void PS_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank, const std::string pointChoosingMethod = "Chebyshev",const int minRank = -1);
  
void PS_LowRankApprox_Sp(const Eigen::SparseMatrix<double> & matrixData_Sp,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j,const int numRows, const int numCols, const double tolerance, int &calculatedRank);

template <typename T>
void SVD_LowRankApprox(const T & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank, const int minRank = -1);

int  SVD_LowRankApprox(const Eigen::MatrixXd & matrixData, const double accuracy, Eigen::MatrixXd* Wptr = NULL, Eigen::MatrixXd* Vptr = NULL, Eigen::MatrixXd* Kptr = NULL, int minRank = -1); 

template <typename T>
void PS_Boundary_LowRankApprox(const T & matrixData,const Eigen::SparseMatrix<double> graphData,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j, const int numRows, const int numCols,const double tolerance,int & calculatedRank, const int maxDepth = 2,const std::string savePath = "none",int numSel = 2);  

template <typename T>
void extractRowsCols(const T & matrixData, int min_i,int min_j,int numRows,int numCols,Eigen::MatrixXd & W,Eigen::MatrixXd & V,const std::vector<int> & rowIndex,const std::vector<int> & colIndex,const double tolerance, int & calculatedRank,const std::string mode = "fullPivLU");

int  getBoundaryRowColIdx(const Eigen::SparseMatrix<double>  & graphData,const int min_i, const int min_j,const int numRows,const int numCols,const int depth,std::vector<int> & rowIdx,std::vector<int> & colIdx,int numSel = 2);

int add_LR(Eigen::MatrixXd & result_U,Eigen::MatrixXd & result_V,const Eigen::MatrixXd & U1, const Eigen::MatrixXd & V1, const Eigen::MatrixXd & U2, const Eigen::MatrixXd & V2,double tol,std::string mode);

int PS_PseudoInverse(Eigen::MatrixXd & colMatrix,Eigen::MatrixXd & rowMatrix, Eigen::MatrixXd & U, Eigen::MatrixXd & V,std::vector<int> rowIdxVec,const  double tol,const std::string mode);

#endif
