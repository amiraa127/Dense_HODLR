#ifndef HELPERFUNCTIONS_HODLR_SOLVER_HPP
#define HELPERFUNCTIONS_HODLR_SOLVER_HPP

//C++ Dependencies
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

//External Dependencies
#include <Eigen/Dense>
#include <Eigen/Sparse>

//Custom Dependencies
#include "HODLR_Matrix.hpp"
#include "user_IndexTree.hpp"

/******************Pre-programmed radial basis function kernels:**********************/

//1+r^2
double quadraticKernel(const double r);

//sqrt(1+r^2)
double multiQuadraticKernel(const double r);

//1/(1+r^2)
double inverseQuadraticKernel(const double r);

//1/sqrt(1+r^2)
double inverseMultiQuadraticKernel(const double r);

//exp(-r^2)
double gaussianKernel(const double r);

//exp(-r)
double exponentialKernel(const double r);

//log(1+r)
double logarithmicKernel(const double r);

//1/r
double oneOverRKernel(const double r);

//1/sqrt(r)
double oneOverSqRKernel(const double r);

//log(r)
double logRKernel(const double r); 

/**************************************************************************************/


std::vector<int> createUniqueRndIdx(const int min, const int max,const int n);

std::vector<int> createSequentialVec(const int min,const int size);

/* Function: makeMatrixFrom1DInterval
 * ----------------------------------
 * This function creates a dense interaction matrix given a set of row and column points given a kernel.
 * The points must lie on a 1D manifold.
 * rowPts : Coordinates of the row points.
 * colPts : Coordinates of the column points.
 * diagValue : The diagonal entry value of the dense matrix.
 * kernel : Pointer to the kernel function.
 */
Eigen::MatrixXd makeMatrixFrom1DInterval(const std::vector<double> rowPts,const  std::vector<double> colPts,const  double diagValue ,double (*kernel)(const double));


/* Function: makeMatrix1DUniformPts
 * --------------------------------
 * This function creates a dense interaction matrix arising from the interaction of uniform points on a 1D interval.
 * This function is mostly a wrapper for makeMatrixFrom1DInterval.
 * minRowPt : Lower bound of the 1D interval for the row points.
 * maxRowPt : Upper bound for the 1D interval for the row points.
 * minColPt : Lower bound of the 1D interval for the column points.
 * maxColpt : Upper bound of the 1D interval for the row points.
 * numRows : Number of rows (row points).
 * numCols : Number of columns (column points).
 * diagValue : The diagonal entry value of the dense matrix.   
`* kernel : Pointer to the kernel function.  
 */
Eigen::MatrixXd makeMatrix1DUniformPts(const int minRowPt, const int maxRowPt, const int minColPt, const int maxColPt, const int numRows, const int numCols, const double diagValue, double (*kernel)(const double)); 


/* Function: testACASolverConv1DUniformPts
 * ---------------------------------------
 * This function creates a convergence plot of solver relative error vs ACA tolerance for a dense self interaction matrix.
 * The test dense matrix, is an interaction matrix arising from the self interaction of uniform points on a 1D interval.
 * intervalMin : Lower bound of the 1D interval.
 * intervalMax : Upper bound of the 1D interval.
 * numPts : Number of interacting points (matrix size).
 * diagValue : The diagonal entry value of the dense matrix.
 * exactSoln : The test right hand side of the linear system.
 * outputFileName : Path of the output file.
 * kernel : Pointer to the kernel function. 
 * solverType : Type of HODLR solver. Default is recLU.
 */
void testACASolverConv1DUnifromPts(const double intervalMin,const double intervalMax, const int numPts, const double diagValue, Eigen::VectorXd exactSoln, std::string outputFileName, double (*kernel)(const double r), std::string solverType = "recLU");


/* Function: testACASolverSpeed1DUniformPts
 * ---------------------------------------
 * This function creates a plot of solver CPU Time vs matrix size for dense self interaction matrice with various sizes.
 * The test dense matrices, are interaction matrix arising from the self interaction of uniform points on a 1D interval.
 * intervalMin : Lower bound of the 1D interval.
 * intervalMax : Upper bound of the 1D interval.
 * diagValue : The diagonal entry value of the dense matrix.
 * ACATolerance : The ACA tolerance of the HODLR solver.
 * outputFileName : Path of the output file.
 * kernel : Pointer to the kernel function. 
 * solverType : Type of HODLR solver. Default is recLU.
 */
void testACASolverSpeed1DUniformPts(const double intervalMin, const double intervalMax,const double diagValue,const double ACATolerance, std::string outputFileName, double (*kernel)(const double), std::string solverType = "recLU");

 /* Function: testACASolverSpeed1DUniformPts_FixedSize
 * ---------------------------------------
 * This function ccalculates the cpu time of the various stages in a solve process for a fixed-size matrix.
 * The test dense matrices, are interaction matrix arising from the self interaction of uniform points on a 1D interval.
 * intervalMin : Lower bound of the 1D interval.
 * intervalMax : Upper bound of the 1D interval.
 * diagValue : The diagonal entry value of the dense matrix.
 * LR_Tolerance : The ACA tolerance of the HODLR solver.
 * outputFileName : Path of the output file.
 * kernel : Pointer to the kernel function. 
 * matrixSize : Size of the matrix.
 * solverType : Type of HODLR solver. Default is recLU.
 */
void testACASolverSpeed1DUniformPts_FixedSize(const double intervalMin, const double intervalMax,const double diagValue,const double LR_Tolerance, std::string outputFileName, double (*kernel)(const double), const int matrixSize, std::string solverType);

double eigenPartialPivLUSpeed(const Eigen::MatrixXd & inputMatrix);

void testSolverSpeed(const std::string inputFilePath,const std::string outputFilePath,const int sizeThreshold,std::string solverType,user_IndexTree & usrTree);

void testBoundaryLRSolver(const std::string inputMatrixFileName,const std::string inputGraphFileName,const std::string outputFileName,const double iterInitTol,const int sizeThreshold,const int depth);

void testBoundaryLRSolver(const std::string inputMatrixFileName,const std::string inputGraphFileName,const std::string outputFileName,const double iterInitTol,const int sizeThreshold,const int depth,user_IndexTree &usrTree);

void analyzeRank(const std::string inputMatrixFileName,const std::string inputGraphFileName,const std::string outputFileName,const int input_Min_i,const int input_Min_j,const int input_NumRows, const int input_NumCols,std::string mode = "default");

#endif
