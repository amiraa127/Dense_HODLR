#include "helperFunctions.hpp"

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


void testBoundaryLRSolver(const std::string inputMatrixFileName,const std::string inputGraphFileName,const std::string outputFileName,const double iterInitTol,const int sizeThreshold,const int depth){
  Eigen::MatrixXd inputMatrix = readBinaryIntoMatrixXd(inputMatrixFileName);
  Eigen::SparseMatrix<double> inputGraph = readMtxIntoSparseMatrix(inputGraphFileName);
  HODLR_Matrix testHODLR(inputMatrix,inputGraph,sizeThreshold);
  int matrixSize = inputMatrix.rows();
  Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
  Eigen::VectorXd inputF    = inputMatrix * exactSoln;
  testHODLR.set_BoundaryDepth(depth);
  testHODLR.printResultInfo = true;
  Eigen::VectorXd solverSoln = testHODLR.iterative_Solve(inputF,100,1e-10,iterInitTol,"PS_Boundary","recLU");
  Eigen::VectorXd difference = solverSoln - exactSoln;
  //double relError = difference.norm()/exactSoln.norm();
  //std::cout<<relError<<std::endl;
  testHODLR.saveSolverInfo(outputFileName);
  double startTime = clock();
  Eigen::PartialPivLU<Eigen::MatrixXd> LU (inputMatrix);
  solverSoln = LU.solve(inputF);
  double endTime = clock();
  double LU_SolveTime = (endTime - startTime)/CLOCKS_PER_SEC;
  std::cout<<"LU Solve Time = "<<LU_SolveTime<<" seconds"<<std::endl;
  std::cout<<"Matrix Size   = "<<inputMatrix.rows()<<std::endl;
}

void testBoundaryLRSolver(const std::string inputMatrixFileName,const std::string inputGraphFileName,const std::string outputFileName,const double iterInitTol,const int sizeThreshold,const int depth,user_IndexTree &usrTree){
  Eigen::MatrixXd inputMatrix = readBinaryIntoMatrixXd(inputMatrixFileName);
  Eigen::SparseMatrix<double> inputGraph = readMtxIntoSparseMatrix(inputGraphFileName);
  HODLR_Matrix testHODLR(inputMatrix,inputGraph,sizeThreshold,usrTree);
  int matrixSize = inputMatrix.rows();
  Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
  Eigen::VectorXd inputF    = inputMatrix * exactSoln;
  testHODLR.set_BoundaryDepth(depth);
  testHODLR.printResultInfo = true;
  Eigen::VectorXd solverSoln = testHODLR.iterative_Solve(inputF,100,1e-10,iterInitTol,"PS_Boundary","recLU");
  Eigen::VectorXd difference = solverSoln - exactSoln;
  //double relError = difference.norm()/exactSoln.norm();
  //std::cout<<relError<<std::endl;
  testHODLR.saveSolverInfo(outputFileName);
  double startTime = clock();
  Eigen::PartialPivLU<Eigen::MatrixXd> LU (inputMatrix);
  solverSoln = LU.solve(inputF);
  double endTime = clock();
  double LU_SolveTime = (endTime - startTime)/CLOCKS_PER_SEC;
  std::cout<<"LU Solve Time = "<<LU_SolveTime<<" seconds"<<std::endl;
  std::cout<<"Matrix Size   = "<<inputMatrix.rows()<<std::endl;
}

void analyzeRank(const std::string inputMatrixFileName,const std::string inputGraphFileName,const std::string outputFileName,const int input_Min_i,const int input_Min_j, const int input_NumRows,const int input_NumCols,std::string mode){
  Eigen::MatrixXd inputMatrix = readBinaryIntoMatrixXd(inputMatrixFileName);
  Eigen::SparseMatrix<double> inputGraph = readMtxIntoSparseMatrix(inputGraphFileName);
  int min_i,min_j,numRows,numCols;
  int matrixSize = inputMatrix.rows();
  if (mode == "topOffDiag"){  
    int split = matrixSize/2; 
    min_i = 0;
    min_j = split + 1;
    numRows = split + 1;
    numCols = matrixSize - split - 1;
  } else if (mode == "bottOffDiag"){  
    int split = matrixSize/2; 
    min_i = split + 1;
    min_j = 0;
    numRows = matrixSize - split - 1;
    numCols = split + 1; 
  }else{
    min_i   = input_Min_i;
    min_j   = input_Min_j;
    numRows = input_NumRows;
    numCols = input_NumCols;
  }
  int currRank  = 0;
  int currRankComp = 0;
  int nextRank  = -1;
  int depth     = 0;
  Eigen::MatrixXd currBlock = inputMatrix.block(min_i,min_j,numRows,numCols);
  Eigen::MatrixXd U,V,K;
  std::vector<double> numPoints,boundaryError,boundaryErrorComp,SVDError,singularValues;
  while (currRank != nextRank){
    nextRank = currRank;
    PS_Boundary_LowRankApprox(inputMatrix,inputGraph,U,V,K,min_i,min_j,numRows,numCols,1e-15,currRank,depth);
    numPoints.push_back(currRank);
    double relError = (U * K * V.transpose() - currBlock).norm()/currBlock.norm();
    boundaryError.push_back(relError);
    PS_Boundary_LowRankApprox(inputMatrix,inputGraph,U,V,K,min_i,min_j,numRows,numCols,1e-1,currRankComp,depth);
    relError = (U * K * V.transpose() - currBlock).norm()/currBlock.norm();
    boundaryErrorComp.push_back(relError);
    depth ++;
  }
  saveVectorAsText(outputFileName + "numPointsVsBoundaryDistance",numPoints);
  saveVectorAsText(outputFileName + "boundaryErrorVsBoundaryDistance",boundaryError);
  saveVectorAsText(outputFileName + "boundaryErrorCompVsBoundaryDistance",boundaryErrorComp);
  // Use svd to calculate the optimal low-rank approximation error
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(currBlock,Eigen::ComputeThinU|Eigen::ComputeThinV);
  Eigen::VectorXd s = svd.singularValues();
  for (unsigned int i = 0; i < s.rows(); i++)
    singularValues.push_back(s(i));

  for (unsigned int i = 0; i < numPoints.size(); i++){
    int rank = numPoints[i];
    U = svd.matrixU().leftCols(rank);
    V = svd.matrixV().leftCols(rank);
    K = Eigen::MatrixXd::Zero(rank,rank);
    for (int j = 0; j < rank; j++)
      K(j,j) = singularValues[j];
    double relError = (U * K * V.transpose() - currBlock).norm()/currBlock.norm();
    SVDError.push_back(relError);
  }
  saveVectorAsText(outputFileName + "SVDErrorVsBoundaryDistance",SVDError);
  saveVectorAsText(outputFileName + "SingularValueDecay",singularValues);

  // Calculate distance from boundary vs pivot size
  PS_Boundary_LowRankApprox(inputMatrix,inputGraph,U,V,K,min_i,min_j,numRows,numCols,1e-15,currRank,matrixSize,outputFileName + "distanceFromBoundaryVsPivotSize");
}
