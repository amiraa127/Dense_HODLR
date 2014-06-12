#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "HODLR_Matrix.hpp"
#include <cmath>
#include <ctime>

double quadraticKernel(const double r){
  return (1 + pow(r,2));
}

double multiQuadraticKernel(const double r){
  return sqrt(1 + pow(r,2));
}

double inverseQuadraticKernel(const double r){
  return 1/(1+pow(r,2));
}

double inverseMultiQuadraticKernel(const double r){
  return 1/sqrt(1+pow(r,2));
}

double gaussianKernel(const double r){
  return exp(-pow(r,2));
}

double exponentialKernel(const double r){
  return exp(-r);
}

double logarithmicKernel(const double r){
  return log(1+r);
}

double oneOverRKernel(const double r){
  return 1/r;
}

double oneOverSqRKernel(const double r){
  return 1/sqrt(r);
}

double logRKernel(const double r){
  return log(r);
}


/* Function: readTxtIntoMatrix
 * ---------------------------
 * This function reads a text file and outputs a dense matrix (Eigen's MatrixXd) on return.
 * The text file format is as follows:
 * Fisrt line : nRows nCols
 * All other lines: i j value(i,j)
 * This function is very slow for large matrices. Use the binary format instead.
 * inputFileName : Path of the input text file.
 */
Eigen::MatrixXd readTxtIntoMatrix(const std::string inputFileName){
  std::ifstream inputFile;
  inputFile.open(inputFileName.c_str());
  if (!inputFile.fail()){
    std::string currLine;
    int numRows,numCols;
    double value;
    int currRow = 0;
    int currCol = 0;
    int error;
    getline(inputFile,currLine);
    error = sscanf(currLine.c_str(),"%u %u",&numRows,&numCols);
    //check format
    if (error != 2){
      std::cout<<"Error! Bad format."<<std::endl;
      exit(EXIT_FAILURE);
    }
    int numEntries = 0;
    int maxEntries = numRows * numCols;
    Eigen::MatrixXd result(numRows,numCols);
    while ((!inputFile.eof()) && (numEntries < maxEntries)){
      getline(inputFile,currLine);
      error = sscanf(currLine.c_str(),"%u %u %lf",&currRow,&currCol,&value);
      if ((currRow == (numRows - 2)) && (currCol == (numCols - 2)))
	  std::cout<<value<<std::endl;
      //check format
      if (error != 3){
	std::cout<<"Error! Bad format."<<std::endl;
	exit(EXIT_FAILURE);
      }
      result(currRow,currCol) = value;
      numEntries ++;
    }
    inputFile.close();
    return result;
  }else{
    std::cout<<"Error! File "<<inputFileName<<" could not be opened"<<std::endl;
    exit(EXIT_FAILURE); 
  }
}


/* Function: saveMatrixXdToBinary
 * ------------------------------
 * This function saves a dense matrix (Eigen's MatrixXd) as a binary file (SVD_F_DB) file format.
 * inputMatrix : The dense matrix being saved.
 * outputFileName : Path of the output file.
 */
void saveMatrixXdToBinary(const Eigen::MatrixXd& inputMatrix, const std::string outputFileName){
  std::ofstream outputFile;
  int nRows = inputMatrix.rows();
  int nCols = inputMatrix.cols();
  outputFile.open(outputFileName.c_str(),std::ios::binary);
  if (outputFile.is_open()){
    outputFile.write((char*)&nRows,sizeof(int));
    outputFile.write((char*)&nCols,sizeof(int));
    for (int i = 0; i < nRows ;i++)
      for (int j = 0; j< nCols ;j++){
	double currValue = inputMatrix(i,j);
	outputFile.write((char*)&currValue,sizeof(double));
      }
  }
  outputFile.close();
}


/* Function: readBinaryIntoMatrixXd
 * -------------------------------
 * This function reads a dense matrix binary file (SVD_F_DB) and outputs on return, a dense matrix (Eigen's MatrixXd).
 * inputFileName : Path of the input file.
 */
Eigen::MatrixXd readBinaryIntoMatrixXd(const std::string inputFileName){
  std::ifstream inputFile;
  inputFile.open(inputFileName.c_str());
  if (inputFile.is_open()){
    int nRows,nCols;
    inputFile.read((char*)&nRows,sizeof(int));
    inputFile.read((char*)&nCols,sizeof(int));
    Eigen::MatrixXd result(nRows,nCols);
    for (int i = 0; i < nRows ;i++)
      for (int j = 0; j< nCols ;j++){
	double currValue;
	inputFile.read((char*)&currValue,sizeof(double));
	result(i,j) = currValue;
      }
    inputFile.close();
    return result;
  }else{
    std::cout<<"Error! File "<<inputFileName<<" could not be opened"<<std::endl;
    exit(EXIT_FAILURE);
  } 
}



/* Function : readMtxIntoSparseMatrix
 *-------------------------------------
 * This function reads a sparse matrix market format (*.mmx) file and returns an Eigen sparse matrix object.
 * Currently it only supports matrix object type with coordinate format. Only real or double data types are acceptable at this time.
 * The symmetricity can only be general or symmetric.
 * inputFileName : The path of the input matrix market file.
 */
Eigen::SparseMatrix<double> readMtxIntoSparseMatrix(const std::string inputFileName){
  //open the file
  std::ifstream inputFile;
  inputFile.open(inputFileName.c_str());
  if (!inputFile.fail()){
    std::string currLine;
    int numRows,numCols,nnz;
    double value;
    int currRow,currCol;
    int error;
    bool isSymmetric;
    //check header
    char str1[20],str2[20],str3[20],str4[20],str5[20];
    getline(inputFile,currLine);
    error = sscanf(currLine.c_str(),"%s %s %s %s %s",str1,str2,str3,str4,str5);
    if ((error != 5) || (strcmp(str1,"%%MatrixMarket") != 0)){
      std::cout<<"Error! Incorrect file header."<<std::endl;
      exit(EXIT_FAILURE);
    }
    if (strcmp(str2,"matrix") != 0){
      std::cout<<"Error! Only matrix object type is acceptable at this time."<<std::endl;
      exit(EXIT_FAILURE);
    }
    if (strcmp(str3,"coordinate") != 0){
      std::cout<<"Error! Only coordinate format is acceptable at this time."<<std::endl;
      exit(EXIT_FAILURE);
    }
    if ((strcmp(str4,"real") != 0) && (strcmp(str4,"double") != 0)){
      std::cout<<"Error! Only real or double data types are acceptable at this time."<<std::endl;
      exit(EXIT_FAILURE);
    }
    if ((strcmp(str5,"general") == 0))
      isSymmetric = false;
    else if ((strcmp(str5,"symmetric") == 0))
      isSymmetric = true;
    else{
      std::cout<<"Error! Only general or symmetric symmetry types are acceptable at this time."<<std::endl;
      exit(EXIT_FAILURE);
    } 
      
    //start filling the matrix
    while (inputFile.peek() == '%')
      inputFile.ignore(2048,'\n');
    getline(inputFile,currLine);
    error = sscanf(currLine.c_str(),"%u %u %u",&numRows,&numCols,&nnz);
    //check format correctness
    if (error != 3){
      std::cout<<"Error! Bad format."<<std::endl;
      exit(EXIT_FAILURE);
    }
    Eigen::SparseMatrix<double> result(numRows,numCols);
    std::vector<Eigen::Triplet<double,int> > tripletVector;
    int numEntries = 0;
    while ((!inputFile.eof()) && (numEntries < nnz)){
      getline(inputFile,currLine);
      error = sscanf(currLine.c_str(),"%u %u %lf",&currRow,&currCol,&value);
      //check format correctness
      if (error != 3){
	std::cout<<"Error! Bad format."<<std::endl;
	exit(EXIT_FAILURE);
      }
      Eigen::Triplet<double,int> currTriplet(currRow-1,currCol-1,value);
      tripletVector.push_back(currTriplet);
      // push back adjoint value into the matrix
      if (isSymmetric){
	Eigen::Triplet<double,int> adjTriplet(currCol-1,currRow-1,value);
	tripletVector.push_back(adjTriplet);
      }
      numEntries++;
    }
    inputFile.close();
    result.setFromTriplets(tripletVector.begin(),tripletVector.end());
    return result;
  }else{
    std::cout<<"Error! File "<<inputFileName<<" could not be opened."<<std::endl;
    exit(EXIT_FAILURE);
  }
}

/* Function: makeMatrixFrom1DInterval
 * ----------------------------------
 * This function creates a dense interaction matrix given a set of row and column points given a kernel.
 * The points must lie on a 1D manifold.
 * rowPts : Coordinates of the row points.
 * colPts : Coordinates of the column points.
 * diagValue : The diagonal entry value of the dense matrix.
 * kernel : Pointer to the kernel function
 */
Eigen::MatrixXd makeMatrixFrom1DInterval(const std::vector<double> rowPts,const std::vector<double> colPts,const double diagValue ,double (*kernel)(const double)){
  int numRows = rowPts.size();
  int numCols = colPts.size();
  assert((numRows > 0) && (numCols > 0));
  Eigen::MatrixXd result(numRows,numCols);
  for (int i = 0; i < numRows; i++)
    for (int j = 0; j< numCols; j++){
      double r = std::abs(rowPts[i] - colPts[j]);
      if (i == j)
	result(i,i) = diagValue;
      else
	result(i,j) = kernel(r);
    }
  return result;
}


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
Eigen::MatrixXd makeMatrix1DUniformPts(const int minRowPt, const int maxRowPt, const int minColPt, const int maxColPt, const int numRows, const int numCols, const double diagValue, double (*kernel)(const double)){
  assert((maxRowPt > minRowPt) && (maxColPt > minColPt));
  assert((numRows > 0) && (numCols > 0));
  std::vector<double> rowPts;
  std::vector<double> colPts;
  // Fill row points
  double rowStride,colStride;
  if (numRows != 1)
    rowStride = (maxRowPt - minRowPt + 0.0)/(numRows -1);
  else
    rowStride = 0;
  for (int i = 0; i < numRows; i++)
    rowPts.push_back(minRowPt + rowStride * i);
  // Fill col points
  if (numCols != 1)
    colStride = (maxColPt - minColPt + 0.0)/(numCols -1);
  else
    colStride = 0;
  for (int i = 0; i < numCols; i++)
    colPts.push_back(minColPt + colStride * i);

  // Make matrix
  return makeMatrixFrom1DInterval(rowPts,colPts,diagValue,kernel);
}


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
void testACASolverConv1DUnifromPts(const double intervalMin,const double intervalMax, const int numPts, const double diagValue, Eigen::VectorXd exactSoln, std::string outputFileName, double (*kernel)(const double r),std::string solverType){
  assert(intervalMax > intervalMin);
  assert(numPts > 0);
  assert(exactSoln.rows() == numPts);
  int minTol = -5;
  int maxTol = -10;
  int sizeThreshold = 30;
  
  Eigen::MatrixXd denseMatrix = makeMatrix1DUniformPts (intervalMin, intervalMax, intervalMin, intervalMax, numPts, numPts, diagValue, kernel);
  Eigen::VectorXd RHS = denseMatrix * exactSoln;
  HODLR_Matrix denseHODLR(denseMatrix, sizeThreshold);
  std::ofstream outputFile;
  outputFile.open(outputFileName.c_str());
  for (int i = minTol; i >= maxTol; i--){
    double tol = pow(10,i);
    denseHODLR.set_LRTolerance(tol);
    Eigen::VectorXd solverSoln;
    if (solverType == "recLU")
      solverSoln = denseHODLR.recLU_Solve(RHS);
    if (solverType == "extendedSp")
      solverSoln = denseHODLR.extendedSp_Solve(RHS);
    Eigen::VectorXd residual = solverSoln-exactSoln;
    double relError = residual.norm()/exactSoln.norm();
    outputFile<<tol<<"       "<<relError<<std::endl;
  }
  outputFile.close();
}


/* Function: testACASolverSpeed1DUniformPts
 * ---------------------------------------
 * This function creates a plot of solver CPU Time vs matrix size for dense self interaction matrice with various sizes.
 * The test dense matrices, are interaction matrix arising from the self interaction of uniform points on a 1D interval.
 * intervalMin : Lower bound of the 1D interval.
 * intervalMax : Upper bound of the 1D interval.
 * diagValue : The diagonal entry value of the dense matrix.
 * LR_Tolerance : The ACA tolerance of the HODLR solver.
 * outputFileName : Path of the output file.
 * kernel : Pointer to the kernel function. 
 * solverType : Type of HODLR solver. Default is recLU.
 */
void testACASolverSpeed1DUniformPts(const double intervalMin, const double intervalMax,const double diagValue,const double LR_Tolerance, std::string outputFileName, double (*kernel)(const double), std::string solverType){
  assert(intervalMax > intervalMin);
  int sizeThreshold = 30;
  int numIterations = 5;
  int minSize = 0;
  int maxSize = 4;

  std::ofstream outputFile_Total;
  std::ofstream outputFile_Factorization;
  std::ofstream outputFile_LR;
  std::ofstream outputFile_Solve;
  std::ofstream outputFile_Assembly;

  outputFile_Total.open(outputFileName.c_str());
  outputFile_Factorization.open((outputFileName + "_Factorization").c_str());
  outputFile_LR.open((outputFileName + "_LR").c_str());
  outputFile_Solve.open((outputFileName + "_Solve").c_str());
  if (solverType == "extendedSp")
    outputFile_Assembly.open((outputFileName + "_Assembly").c_str());

  for (int i = minSize; i <= maxSize; i++){
    int matrixSize = 1000 * pow(2,i);
    Eigen::MatrixXd denseMatrix = makeMatrix1DUniformPts (intervalMin, intervalMax, intervalMin, intervalMax, matrixSize, matrixSize, diagValue, kernel);
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd RHS = denseMatrix * exactSoln;
   
    double sum_Factorization = 0;
    double sum_LR = 0;
    double sum_Solve = 0;
    double sum_Assembly = 0;
    double sum_Total = 0;
    for (int i = 0; i < numIterations; i++){
      HODLR_Matrix denseHODLR(denseMatrix , sizeThreshold);
      denseHODLR.set_LRTolerance(LR_Tolerance);
      Eigen::VectorXd solverSoln;
      
      if (solverType == "recLU"){
	solverSoln = denseHODLR.recLU_Solve(RHS);
	sum_Factorization += denseHODLR.get_recLU_FactorizationTime();
	sum_LR += denseHODLR.get_LR_ComputationTime();
	sum_Solve += denseHODLR.get_recLU_SolveTime();
	sum_Total += denseHODLR.get_recLU_FactorizationTime() + denseHODLR.get_LR_ComputationTime() + denseHODLR.get_recLU_SolveTime();
      }
      if (solverType == "extendedSp"){
	solverSoln = denseHODLR.extendedSp_Solve(RHS);
	sum_Factorization += denseHODLR.get_extendedSp_FactorizationTime();
	sum_LR += denseHODLR.get_LR_ComputationTime();
	sum_Solve += denseHODLR.get_extendedSp_SolveTime();
	sum_Assembly += denseHODLR.get_extendedSp_AssemblyTime();
	sum_Total += denseHODLR.get_extendedSp_FactorizationTime() + denseHODLR.get_LR_ComputationTime() + denseHODLR.get_extendedSp_AssemblyTime() + denseHODLR.get_extendedSp_SolveTime();
      }

    }
    
  outputFile_Factorization<<matrixSize<<"       "<<sum_Factorization/numIterations<<std::endl;
  outputFile_LR<<matrixSize<<"       "<<sum_LR/numIterations<<std::endl;
  outputFile_Solve<<matrixSize<<"       "<<sum_Solve/numIterations<<std::endl;
  outputFile_Total<<matrixSize<<"       "<<sum_Total/numIterations<<std::endl;
  if (solverType == "extendedSp")
    outputFile_Assembly<<matrixSize<<"       "<<sum_Assembly/numIterations<<std::endl; 
  }

  outputFile_Factorization.close();
  outputFile_LR.close();
  outputFile_Solve.close();
  outputFile_Total.close();
  if (solverType == "extendedSp")
    outputFile_Assembly.close(); 
}

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
void testACASolverSpeed1DUniformPts_FixedSize(const double intervalMin, const double intervalMax,const double diagValue,const double LR_Tolerance, std::string outputFileName, double (*kernel)(const double), const int matrixSize, std::string solverType){
  assert(intervalMax > intervalMin);
  int sizeThreshold = 30;
  int numIterations = 5;

  std::ofstream outputFile;
  outputFile.open(outputFileName.c_str());
 
  Eigen::MatrixXd denseMatrix = makeMatrix1DUniformPts (intervalMin, intervalMax, intervalMin, intervalMax, matrixSize, matrixSize, diagValue, kernel);
  Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
  Eigen::VectorXd RHS = denseMatrix * exactSoln;
   
  double sum_Factorization = 0;
  double sum_LR = 0;
  double sum_Solve = 0;
  double sum_Assembly = 0;
  double sum_Total = 0;
  for (int i = 0; i < numIterations; i++){
    HODLR_Matrix denseHODLR(denseMatrix , sizeThreshold);
    denseHODLR.set_LRTolerance(LR_Tolerance);
    Eigen::VectorXd solverSoln;
    
    if (solverType == "recLU"){
      solverSoln = denseHODLR.recLU_Solve(RHS);
      sum_Factorization += denseHODLR.get_recLU_FactorizationTime();
      sum_LR += denseHODLR.get_LR_ComputationTime();
      sum_Solve += denseHODLR.get_recLU_SolveTime();
      sum_Total += denseHODLR.get_recLU_FactorizationTime() + denseHODLR.get_LR_ComputationTime() + denseHODLR.get_recLU_SolveTime();
    }
    if (solverType == "extendedSp"){
      solverSoln = denseHODLR.extendedSp_Solve(RHS);
      sum_Factorization += denseHODLR.get_extendedSp_FactorizationTime();
      sum_LR += denseHODLR.get_LR_ComputationTime();
      sum_Solve += denseHODLR.get_extendedSp_SolveTime();
      sum_Assembly += denseHODLR.get_extendedSp_AssemblyTime();
      sum_Total += denseHODLR.get_extendedSp_FactorizationTime() + denseHODLR.get_LR_ComputationTime() + denseHODLR.get_extendedSp_AssemblyTime() + denseHODLR.get_extendedSp_SolveTime();
    }
      
  }
  outputFile<<"Factorization   "<<"       "<<sum_Factorization/numIterations<<std::endl;
  outputFile<<"LR-Approximation"<<"       "<<sum_LR/numIterations<<std::endl;
  outputFile<<"Solve           "<<"       "<<sum_Solve/numIterations<<std::endl;
  outputFile<<"Total           "<<"       "<<sum_Total/numIterations<<std::endl;
  if (solverType == "extendedSp")
    outputFile<<"Assembly        "<<"       "<<sum_Assembly/numIterations<<std::endl; 
  
  outputFile.close();
}
