#ifndef LOW_RANK_HPP
#define LOW_RANK_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
#include <assert.h>
#include <iostream>

double fullPivACA_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int minPivot = 0);

double partialPivACA_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int minPivot = 0);

  
void PS_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank, const std::string pointChoosingMethod = "Chebyshev",const int minRank = -1);
  
void PS_LowRankApprox_Sp(const Eigen::SparseMatrix<double> & matrixData_Sp,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j,const int numRows, const int numCols, const double tolerance, int &calculatedRank);

void SVD_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numCols, const int numRows, const double tolerance, int & calculatedRank, const int minRank = -1);

int SVD_LowRankApprox(const Eigen::MatrixXd & inputMatrix, const double accuracy, Eigen::MatrixXd* Wptr = NULL, Eigen::MatrixXd* Vptr = NULL, Eigen::MatrixXd* Kptr = NULL, int minRank = -1);  

void identifyBoundary(const Eigen::SparseMatrix<double> & inputGraph,const std::set<int> &rowSet,const std::set<int> &colSet,std::map<int,std::vector<int> > & rowPos,std::map<int,std::vector<int> > & colPos);


#endif
