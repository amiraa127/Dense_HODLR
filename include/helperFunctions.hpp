#ifndef HELPERFUNCTIONS_HODLR_SOLVER_HPP
#define HELPERFUNCTIONS_HODLR_SOLVER_HPP

#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>

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

/* Function: readTxtIntoMatrix
 * ---------------------------
 * This function reads a text file and outputs a dense matrix (Eigen's MatrixXd) on return.
 * The text file format is as follows:
 * Fisrt line : nRows nCols
 * All other lines: i j value(i,j)
 * This function is very slow for large matrices. Use the binary format instead.
 * inputFileName : Path of the input text file.
 */
Eigen::MatrixXd readTxtIntoMatrix(const std::string inputFileName);


/* Function: saveMatrixXdToBinary
 * ------------------------------
 * This function saves a dense matrix (Eigen's MatrixXd) as a binary file (SVD_F_DB) file format.
 * inputMatrix : The dense matrix being saved.
 * outputFileName : Path of the output file.
 */
void saveMatrixXdToBinary(const Eigen::MatrixXd& inputMatrix, const std::string outputFileName);


/* Function: readBinaryIntoMatrixXd
 * --------------------------------
 * This function reads a dense matrix binary file (SVD_F_DB) and outputs on return, a dense matrix (Eigen's MatrixXd).
 * inputFileName : Path of the input file.
 */
Eigen::MatrixXd readBinaryIntoMatrixXd(const std::string inputFileName);


/* Function : readMtxIntoSparseMatrix
 *-------------------------------------
 * This function reads a sparse matrix market format (*.mmx) file and returns an Eigen sparse matrix object.
 * Currently it only supports matrix object type with coordinate format. Only real or double data types are acceptable at this time.
 * The symmetricity can only be general or symmetric.
 * inputFileName : The path of the input matrix market file.
 */
Eigen::SparseMatrix<double> readMtxIntoSparseMatrix(const std::string inputFileName);


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
#endif
