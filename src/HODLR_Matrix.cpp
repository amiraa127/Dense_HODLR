#include "HODLR_Matrix.hpp"

void HODLR_Matrix::setDefaultValues(){
  
  LR_Tolerance = 1e-6;
  minPivot  = 0;
  
  recLU_FactorizationTime      = 0;
  recLU_SolveTime              = 0;
  recLU_TotalTime              = 0;
  recSM_FactorizationTime      = 0;
  recSM_SolveTime              = 0;
  recSM_TotalTime              = 0;
  LR_ComputationTime           = 0;
  extendedSp_Size              = 0;
  extendedSp_AssemblyTime      = 0;
  extendedSp_FactorizationTime = 0;
  extendedSp_SolveTime         = 0;
  extendedSp_TotalTime         = 0;
  totalIter_SolveTime          = 0;
  matrixSize                   = 0; 
  matrixNumRows                = 0;
  matrixNumCols                = 0; 
  constLeafSize                = 0;
  boundaryDepth                = 1;

  LRStoredInTree         = false;
  recLU_Factorized       = false;
  recSM_Factorized       = false;
  assembled_ExtendedSp   = false;
  saveExtendedSp_Matrix  = false;
  freeMatrixMemory       = false;
  freeMatrixMemory_Sp    = false;
  freeGraphMemmory       = false;
  matrixDataAvail        = false;
  matrixDataAvail_Sp     = false;
  graphDataAvail         = false;
  kernelDataAvail        = false;
  isSquareMatrix         = false;
  calculatedDet          = false;

  isLeafConst            = false;
  constLeafSet           = false;
  constLeafFactorized    = false;

  printLevelRankInfo     = false;
  printLevelAccuracy     = false;
  printLevelInfo         = false;
  printResultInfo        = false;  
  extendedSp_SavePath    = "/extendedSp_Data";
}


HODLR_Matrix::HODLR_Matrix(){
  setDefaultValues();
}


HODLR_Matrix::HODLR_Matrix(Eigen::MatrixXd &inputMatrix){
  setDefaultValues();
  matrixData = inputMatrix;
  sizeThreshold = 30;
  matrixDataAvail = true;
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  matrixSize     = inputMatrix.rows();
  matrixNumRows  = inputMatrix.rows();
  matrixNumCols  = inputMatrix.cols();
}

HODLR_Matrix::HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,std::string LR_Method){
  setDefaultValues();
  matrixData_Sp = inputMatrix;
  sizeThreshold = 30;
  matrixDataAvail_Sp = true;
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  indexTree.set_LRMethod(LR_Method);
  matrixSize    = inputMatrix.rows();
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
}

HODLR_Matrix::HODLR_Matrix(Eigen::MatrixXd &inputMatrix,Eigen::SparseMatrix<double> &inputGraph){
  setDefaultValues();
  assert(inputMatrix.cols() == inputGraph.cols());
  assert(inputMatrix.rows() == inputGraph.rows());
  graphData     = inputGraph;
  matrixData    = inputMatrix;
  sizeThreshold = 30;
  graphDataAvail      = true;
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  matrixSize     = inputMatrix.rows();
  matrixNumRows  = inputMatrix.rows();
  matrixNumCols  = inputMatrix.cols();
  indexTree.set_LRMethod("PS_Boundary");
}

HODLR_Matrix::HODLR_Matrix(Eigen::MatrixXd &inputMatrix,int inputSizeThreshold){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true); // Currently unable to build trees for non squared matrices
  matrixData    = inputMatrix;
  matrixSize    = inputMatrix.rows();
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.createDefaultTree(matrixSize);
  initializeInfoVecotrs(indexTree.get_numLevels());
  matrixDataAvail = true;  
}

HODLR_Matrix::HODLR_Matrix(int numRows, int numCols,double (*inputKernel)(int i,int j,void* kernelData),void* inputKernelData,int inputSizeThreshold){
  setDefaultValues();
  isSquareMatrix = (numRows == numCols);
  assert(isSquareMatrix == true); // Currently unable to build trees for non squared matrices
  kernelMatrixData = kernelMatrix(numRows,numCols,inputKernel,inputKernelData);
  matrixSize    = numRows;
  matrixNumRows = numRows;
  matrixNumCols = numCols;
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.createDefaultTree(matrixSize);
  initializeInfoVecotrs(indexTree.get_numLevels());
  kernelDataAvail = true;  
}

HODLR_Matrix::HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,int inputSizeThreshold,std::string LR_Method){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true); // Currently unable to build trees for non squared matrices
  matrixData_Sp = inputMatrix;
  matrixSize    = inputMatrix.rows();
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.set_LRMethod(LR_Method);
  indexTree.createDefaultTree(matrixSize);
  initializeInfoVecotrs(indexTree.get_numLevels());
  matrixDataAvail_Sp = true;
}


HODLR_Matrix::HODLR_Matrix(Eigen::MatrixXd &inputMatrix,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true); // Currently unable to build trees for non squared matrices
  assert(inputMatrix.cols() == inputGraph.cols());
  assert(inputMatrix.rows() == inputGraph.rows());
  graphData     = inputGraph;
  matrixData    = inputMatrix;
  matrixSize    = inputMatrix.rows();
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.set_LRMethod("PS_Boundary");
  indexTree.createDefaultTree(matrixSize);
  initializeInfoVecotrs(indexTree.get_numLevels());
  matrixDataAvail     = true;
  graphDataAvail      = true;
}

HODLR_Matrix::HODLR_Matrix(int numRows, int numCols,double (*inputKernel)(int i,int j,void* inputKernelData),void* inputKernelData,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold){
  setDefaultValues();
  isSquareMatrix = (numRows == numCols);
  assert(isSquareMatrix == true); // Currently unable to build trees for non squared matrices
  assert(numCols == inputGraph.cols());
  assert(numRows == inputGraph.rows());
  graphData           = inputGraph;
  kernelMatrixData    = kernelMatrix(numRows,numCols,inputKernel,inputKernelData);;
  matrixSize          = numRows;
  matrixNumRows       = numRows;
  matrixNumCols       = numCols;
  sizeThreshold       = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.set_LRMethod("PS_Boundary");
  indexTree.createDefaultTree(matrixSize);
  initializeInfoVecotrs(indexTree.get_numLevels());
  kernelDataAvail     = true;
  graphDataAvail      = true;
}
 
HODLR_Matrix::HODLR_Matrix(Eigen::MatrixXd &inputMatrix, int inputSizeThreshold, user_IndexTree &input_IndexTree){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true); // Currently unable to build trees for non squared matrices
  matrixData    = inputMatrix;
  matrixSize    = inputMatrix.rows();
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.createFromUsrTree(matrixSize,input_IndexTree);
  initializeInfoVecotrs(indexTree.get_numLevels());
  matrixDataAvail = true;
}

HODLR_Matrix::HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix, int inputSizeThreshold, user_IndexTree &input_IndexTree,std::string LR_Method){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true);  // Currently unable to build trees for non squared matrices
  matrixData_Sp = inputMatrix;
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  matrixSize    = inputMatrix.rows();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.set_LRMethod(LR_Method);
  indexTree.createFromUsrTree(matrixSize,input_IndexTree);
  initializeInfoVecotrs(indexTree.get_numLevels());
  matrixDataAvail_Sp = true;
}

HODLR_Matrix::HODLR_Matrix(int numRows, int numCols,double (*inputKernel)(int i,int j,void* inputKernelData),void* inputKernelData,int inputSizeThreshold,user_IndexTree &input_IndexTree){
  setDefaultValues();
  isSquareMatrix = (numRows == numCols);
  assert(isSquareMatrix == true); // Currently unable to build trees for non squared matrices
  kernelMatrixData    = kernelMatrix(numRows,numCols,inputKernel,inputKernelData);;
  matrixSize          = numRows;
  matrixNumRows       = numRows;
  matrixNumCols       = numCols;
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.createFromUsrTree(matrixSize,input_IndexTree);
  initializeInfoVecotrs(indexTree.get_numLevels());
  kernelDataAvail = true;
}

HODLR_Matrix::HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix, Eigen::SparseMatrix<double> & inputGraph,int inputSizeThreshold, user_IndexTree &input_IndexTree,std::string LR_Method){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true);  // Currently unable to build trees for non squared matrices
  matrixData_Sp = inputMatrix;
  graphData     = inputGraph;
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  matrixSize    = inputMatrix.rows();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.set_LRMethod(LR_Method);
  indexTree.createFromUsrTree(matrixSize,input_IndexTree);
  initializeInfoVecotrs(indexTree.get_numLevels());
  matrixDataAvail_Sp = true;
  graphDataAvail     = true;
}

HODLR_Matrix::HODLR_Matrix(Eigen::MatrixXd & inputMatrix,Eigen::SparseMatrix<double> &inputGraph, int inputSizeThreshold, user_IndexTree &input_IndexTree){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true);  // Currently unable to build trees for non squared matrices
  assert(inputMatrix.cols() == inputGraph.cols());
  assert(inputMatrix.rows() == inputGraph.rows());
  matrixData    = inputMatrix;
  //matrixData_Sp = inputGraph;
  graphData     = inputGraph;
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  matrixSize    = inputMatrix.rows();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.set_LRMethod("PS_Boundary");
  indexTree.createFromUsrTree(matrixSize,input_IndexTree);
  initializeInfoVecotrs(indexTree.get_numLevels());
  matrixDataAvail = true;
  graphDataAvail  = true;
  //freeMatrixMemory_Sp = true;
  //freeGraphMemmory = true;
}



HODLR_Matrix::~HODLR_Matrix(){

}

HODLR_Matrix:: HODLR_Matrix(const HODLR_Matrix & rhs){
    
  //public attributes
  printLevelRankInfo = rhs.printLevelRankInfo;
  printLevelAccuracy = rhs.printLevelAccuracy;
  printLevelInfo     = rhs.printLevelInfo;
  printResultInfo    = rhs.printResultInfo;
  
  //private attributes
  sizeThreshold      = rhs.sizeThreshold;
  extendedSp_Size    = rhs.extendedSp_Size;
  matrixSize         = rhs.matrixSize;
  matrixNumRows      = rhs.matrixNumRows;
  matrixNumCols      = rhs.matrixNumCols;
  boundaryDepth      = rhs.boundaryDepth;

  recLU_FactorizationTime      = rhs.recLU_FactorizationTime;
  recLU_SolveTime              = rhs.recLU_SolveTime;
  recLU_TotalTime              = rhs.recLU_TotalTime;
  LR_ComputationTime           = rhs.LR_ComputationTime;
  extendedSp_AssemblyTime      = rhs.extendedSp_AssemblyTime;
  extendedSp_FactorizationTime = rhs.extendedSp_FactorizationTime;
  extendedSp_SolveTime         = rhs.extendedSp_SolveTime;
  extendedSp_TotalTime         = rhs.extendedSp_TotalTime;
  totalIter_SolveTime          = rhs.totalIter_SolveTime;

  recLU_FactorLevelTimeVec     = rhs.recLU_FactorLevelTimeVec;
  recLU_SolveLevelTimeVec      = rhs.recLU_FactorLevelTimeVec;
  LR_ComputationLevelTimeVec   = rhs.LR_ComputationLevelTimeVec;
  iter_IterTimeVec             = rhs.iter_IterTimeVec;
  iter_IterAccuracyVec         = rhs.iter_IterAccuracyVec;
  levelRankAverageVec          = rhs.levelRankAverageVec; 
  levelRankAverageVecCnt       = rhs.levelRankAverageVecCnt;

  LRStoredInTree         = rhs.LRStoredInTree;
  recLU_Factorized       = rhs.recLU_Factorized;
  assembled_ExtendedSp   = rhs.assembled_ExtendedSp;
  saveExtendedSp_Matrix  = rhs.saveExtendedSp_Matrix;
  freeMatrixMemory       = rhs.freeMatrixMemory;
  matrixDataAvail        = rhs.matrixDataAvail;    
  matrixDataAvail_Sp     = rhs.matrixDataAvail_Sp;
  isSquareMatrix         = rhs.isSquareMatrix;

  LR_Tolerance           = rhs.LR_Tolerance;
  minPivot               = rhs.minPivot;

  matrixData             = rhs.matrixData;
  matrixData_Sp          = rhs.matrixData_Sp;
  extendedSp_Solver      = rhs.extendedSp_Solver; 
  extendedSp_SavePath    = rhs.extendedSp_SavePath;
  indexTree              = rhs.indexTree;
  //recLUfactorTree needs to be copied :)) TODO

}

void HODLR_Matrix::initInfoVectors(){
  initializeInfoVecotrs(indexTree.get_numLevels());
}

void HODLR_Matrix::initializeInfoVecotrs(int numLevels){
  recLU_FactorLevelTimeVec   = std::vector<double>(numLevels,0.0);
  recLU_SolveLevelTimeVec    = std::vector<double>(numLevels,0.0);
  LR_ComputationLevelTimeVec = std::vector<double>(numLevels,0.0);
  levelRankAverageVec        = std::vector<double>(numLevels,0.0);
  levelRankAverageVecCnt     = std::vector<double>(numLevels,0.0);
}


void HODLR_Matrix::storeLRinTree(){
  assert(indexTree.rootNode != NULL);
  if (LRStoredInTree == false){
    double startTime = clock();
    assert((matrixDataAvail == true) || (matrixDataAvail_Sp == true) || (kernelDataAvail == true));
    storeLRinTree(indexTree.rootNode);
    double endTime = clock();
    LR_ComputationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    LRStoredInTree = true;
  }
  if (freeMatrixMemory == true)
    freeDenseMatMem();
  if (freeMatrixMemory_Sp == true)
    freeSparseMatMem();
  return;
}

void HODLR_Matrix::storeLRinTree(HODLR_Tree::node* HODLR_Root){

  // Base cases;
  if (HODLR_Root == NULL)
    return;
  if (HODLR_Root->isLeaf == true){
    double startTime = clock();
    int numRows = HODLR_Root->max_i - HODLR_Root->min_i + 1;
    int numCols = HODLR_Root->max_j - HODLR_Root->min_j + 1;
    if (matrixDataAvail_Sp == true){
      if (isLeafConst == true){
	if (constLeafSet == false){
	  constLeaf = Eigen::MatrixXd(matrixData_Sp.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols));
	  constLeafSet = true;
	  constLeafSize = numRows;
	}
	assert(numRows == constLeafSize);
      }else
	HODLR_Root->leafMatrix = Eigen::MatrixXd(matrixData_Sp.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols));
    }else if (matrixDataAvail == true){
      if (isLeafConst == true){
	if (constLeafSet == false){
	  constLeaf = matrixData.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols);
	  constLeafSet = true;
	  constLeafSize = numRows;
	}
	assert(numRows == constLeafSize);
      }else
	HODLR_Root->leafMatrix = matrixData.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols);
    }else if (kernelDataAvail == true){
      if (isLeafConst == true){
	if (constLeafSet == false){
	  constLeaf = kernelMatrixData.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols);
	  constLeafSet = true;
	  constLeafSize = numRows;
	}
	assert(numRows == constLeafSize);
      }else
	HODLR_Root->leafMatrix = kernelMatrixData.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols);
    }else{
      std::cout<<"Error! No matrix data available."<<std::endl;
      exit(EXIT_FAILURE);
    }
    double endTime = clock();
    LR_ComputationLevelTimeVec[HODLR_Root->currLevel] += (endTime - startTime)/CLOCKS_PER_SEC;
    return;
  }
  double startTime = clock();
  int numRows_TopOffDiag  = HODLR_Root->splitIndex_i - HODLR_Root->min_i + 1; 
  int numRows_BottOffDiag = HODLR_Root->max_i - HODLR_Root->splitIndex_i;
  int numCols_TopOffDiag  = HODLR_Root->max_j - HODLR_Root->splitIndex_j;
  int numCols_BottOffDiag = HODLR_Root->splitIndex_j - HODLR_Root->min_j + 1;

  // Calculate the LR factorizations
  if (HODLR_Root->LR_Method == "partialPiv_ACA"){
    
    if (matrixDataAvail == true){
      partialPivACA_LowRankApprox(matrixData,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank,HODLR_Root->topOffDiag_maxRank,minPivot);
      partialPivACA_LowRankApprox(matrixData,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank,HODLR_Root->bottOffDiag_maxRank,minPivot);
    }else if (kernelDataAvail == true){    
      partialPivACA_LowRankApprox(kernelMatrixData,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank,HODLR_Root->topOffDiag_maxRank,minPivot);
      partialPivACA_LowRankApprox(kernelMatrixData,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank,HODLR_Root->bottOffDiag_maxRank,minPivot);
    }else{
      std::cout<<"Error! No matrix data available."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }else if (HODLR_Root->LR_Method == "fullPiv_ACA"){

    if (matrixDataAvail == true){
      fullPivACA_LowRankApprox(matrixData,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank,minPivot);
      fullPivACA_LowRankApprox(matrixData,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank,minPivot);
    }else if (kernelDataAvail == true){
      fullPivACA_LowRankApprox(kernelMatrixData,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank,minPivot);
      fullPivACA_LowRankApprox(kernelMatrixData,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank,minPivot);
    }else{
      std::cout<<"Error! No matrix data available."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }else if (HODLR_Root->LR_Method == "PS_Cheby"){
    
    assert(matrixDataAvail ==  true);   
    PS_LowRankApprox(matrixData,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,"Chebyshev",HODLR_Root->topOffDiag_minRank);
    PS_LowRankApprox(matrixData,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,"Chebyshev",HODLR_Root->bottOffDiag_minRank);
      
  }else if (HODLR_Root->LR_Method == "SVD"){ 
    
    //assert(matrixDataAvail == true);
    //SVD_LowRankApprox(matrixData,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->topOffDiagK,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank);
    //SVD_LowRankApprox(matrixData,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->bottOffDiagK,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank);
    Eigen::MatrixXd tempU,tempV,tempK;
    if (matrixDataAvail == true){
      SVD_LowRankApprox(matrixData,tempU,tempV,tempK,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank);
      HODLR_Root->topOffDiagU  = tempU * tempK;
      HODLR_Root->topOffDiagV  = tempV;
      SVD_LowRankApprox(matrixData,tempU,tempV,tempK,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank);
      HODLR_Root->bottOffDiagU = tempU * tempK;
      HODLR_Root->bottOffDiagV = tempV;
    } else if (kernelDataAvail == true){
      SVD_LowRankApprox(kernelMatrixData,tempU,tempV,tempK,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank);
      HODLR_Root->topOffDiagU  = tempU * tempK;
      HODLR_Root->topOffDiagV  = tempV;
      SVD_LowRankApprox(kernelMatrixData,tempU,tempV,tempK,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank);
      HODLR_Root->bottOffDiagU = tempU * tempK;
      HODLR_Root->bottOffDiagV = tempV;
    } else{
      std::cout<<"Error! No matrix data available"<<std::endl;
      exit(EXIT_FAILURE);
    }
  }else if (HODLR_Root->LR_Method == "identifyBoundary"){
    if (graphDataAvail == true){ 
      getBoundaryRowColIdx(graphData,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,boundaryDepth,HODLR_Root->topOffDiagRowIdx,HODLR_Root->topOffDiagColIdx);
      getBoundaryRowColIdx(graphData,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,boundaryDepth,HODLR_Root->bottOffDiagRowIdx,HODLR_Root->bottOffDiagColIdx);
     
    } else if (matrixDataAvail_Sp == true){
      getBoundaryRowColIdx(matrixData_Sp,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,boundaryDepth,HODLR_Root->topOffDiagRowIdx,HODLR_Root->topOffDiagColIdx);
      getBoundaryRowColIdx(matrixData_Sp,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,boundaryDepth,HODLR_Root->bottOffDiagRowIdx,HODLR_Root->bottOffDiagColIdx);
    } else {
      std::cout<<"Error! No Graph data available."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }else if (HODLR_Root->LR_Method == "PS_Sparse"){
    assert(matrixDataAvail_Sp == true);  
    PS_LowRankApprox_Sp(matrixData_Sp,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank);
    PS_LowRankApprox_Sp(matrixData_Sp,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank);
    if (graphDataAvail == true){ 
      getBoundaryRowColIdx(graphData,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,boundaryDepth,HODLR_Root->topOffDiagRowIdx,HODLR_Root->topOffDiagColIdx);
      getBoundaryRowColIdx(graphData,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,boundaryDepth,HODLR_Root->bottOffDiagRowIdx,HODLR_Root->bottOffDiagColIdx);
    }else if (matrixDataAvail_Sp == true){
      getBoundaryRowColIdx(matrixData_Sp,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,boundaryDepth,HODLR_Root->topOffDiagRowIdx,HODLR_Root->topOffDiagColIdx);
      getBoundaryRowColIdx(matrixData_Sp,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,boundaryDepth,HODLR_Root->bottOffDiagRowIdx,HODLR_Root->bottOffDiagColIdx);
    }else{
      std::cout<<"Error! No Graph data available."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }else if (HODLR_Root->LR_Method == "PS_Boundary"){
      //assert(matrixDataAvail == true);  
    //PS_Boundary_LowRankApprox(matrixData,matrixData_Sp,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,boundaryDepth);
    //PS_Boundary_LowRankApprox(matrixData,matrixData_Sp,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,boundaryDepth);
    assert(graphDataAvail == true);  
    if (matrixDataAvail == true){
      PS_Boundary_LowRankApprox(matrixData,graphData,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,boundaryDepth);
      PS_Boundary_LowRankApprox(matrixData,graphData,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,boundaryDepth);
    } else if (kernelDataAvail == true){
      PS_Boundary_LowRankApprox(kernelMatrixData,graphData,HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag,LR_Tolerance,HODLR_Root->topOffDiagRank,boundaryDepth);
      PS_Boundary_LowRankApprox(kernelMatrixData,graphData,HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag,LR_Tolerance,HODLR_Root->bottOffDiagRank,boundaryDepth);
    }else{
      std::cout<<"Error! No matrix data available."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }else{
    std::cout<<"Error!. Invalid low-rank approximation scheme ( "<<HODLR_Root->LR_Method<<")."<<std::endl;
    exit(EXIT_FAILURE);
  }
  
  double endTime = clock();
  levelRankAverageVec[HODLR_Root->currLevel] = levelRankAverageVec[HODLR_Root->currLevel] * levelRankAverageVecCnt[HODLR_Root->currLevel] + HODLR_Root->topOffDiagRank + HODLR_Root->bottOffDiagRank;
  levelRankAverageVecCnt[HODLR_Root->currLevel] += 2;
  levelRankAverageVec[HODLR_Root->currLevel] /= levelRankAverageVecCnt[HODLR_Root->currLevel];
  LR_ComputationLevelTimeVec[HODLR_Root->currLevel] += (endTime - startTime)/CLOCKS_PER_SEC;
  storeLRinTree(HODLR_Root->left);
  storeLRinTree(HODLR_Root->right);
}

Eigen::MatrixXd HODLR_Matrix::createExactHODLR(const int rank,int input_MatrixSize,const int input_SizeThreshold){
  assert(indexTree.rootNode == NULL);
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(input_MatrixSize,input_MatrixSize);
  isSquareMatrix = true;  // Currently unable to build trees for non squared matrices
  matrixSize    = input_MatrixSize;
  matrixNumRows = input_MatrixSize;
  matrixNumCols = input_MatrixSize;
  sizeThreshold = input_SizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.createDefaultTree(matrixSize);
  initializeInfoVecotrs(indexTree.get_numLevels());
  matrixDataAvail = false;
  createExactHODLR(indexTree.rootNode,rank,result);
  LRStoredInTree = true;
  matrixData = result;
  matrixDataAvail = true;
  return result;

}

void HODLR_Matrix::createExactHODLR(HODLR_Tree::node* HODLR_Root,const int rank,Eigen::MatrixXd & result){
  if (HODLR_Root->isLeaf == true){
    int numRows = HODLR_Root->max_i - HODLR_Root->min_i + 1;
    int numCols = HODLR_Root->max_j - HODLR_Root->min_j + 1;
    HODLR_Root->leafMatrix = Eigen::MatrixXd::Random(numRows,numCols);
    result.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols) =  HODLR_Root->leafMatrix;
    return;
  }
  
  
  int numRows_TopOffDiag  = HODLR_Root->splitIndex_i - HODLR_Root->min_i + 1; 
  int numRows_BottOffDiag = HODLR_Root->max_i - HODLR_Root->splitIndex_i;
  int numCols_TopOffDiag  = HODLR_Root->max_j - HODLR_Root->splitIndex_j;
  int numCols_BottOffDiag = HODLR_Root->splitIndex_j - HODLR_Root->min_j + 1; 

  if (rank > 0){
    HODLR_Root->topOffDiagRank   = rank;
    HODLR_Root->bottOffDiagRank  = rank;
  
    HODLR_Root->topOffDiagU      = Eigen::MatrixXd::Random(numRows_TopOffDiag,rank);
    HODLR_Root->topOffDiagV      = Eigen::MatrixXd::Random(numCols_TopOffDiag,rank);
    HODLR_Root->bottOffDiagU     = Eigen::MatrixXd::Random(numRows_BottOffDiag,rank);
    HODLR_Root->bottOffDiagV     = Eigen::MatrixXd::Random(numCols_BottOffDiag,rank);
    
    HODLR_Root->topOffDiagRowIdx     = createSequentialVec(0,numRows_TopOffDiag);
    HODLR_Root->topOffDiagColIdx     = createSequentialVec(0,numCols_TopOffDiag);
    HODLR_Root->bottOffDiagRowIdx    = createSequentialVec(0,numRows_BottOffDiag);
    HODLR_Root->bottOffDiagColIdx    = createSequentialVec(0,numCols_BottOffDiag);

    result.block(HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1,numRows_TopOffDiag,numCols_TopOffDiag) = HODLR_Root->topOffDiagU * HODLR_Root->topOffDiagV.transpose(); 
    result.block(HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,numRows_BottOffDiag,numCols_BottOffDiag) = HODLR_Root->bottOffDiagU *  HODLR_Root->bottOffDiagV.transpose(); 
  }else{
    HODLR_Root->topOffDiagRank   = 1;
    HODLR_Root->bottOffDiagRank  = 1;
    
    HODLR_Root->topOffDiagU      = Eigen::MatrixXd::Zero(numRows_TopOffDiag,1);
    HODLR_Root->topOffDiagV      = Eigen::MatrixXd::Zero(numCols_TopOffDiag,1);
    HODLR_Root->bottOffDiagU     = Eigen::MatrixXd::Zero(numRows_BottOffDiag,1);
    HODLR_Root->bottOffDiagV     = Eigen::MatrixXd::Zero(numCols_BottOffDiag,1);
    
    HODLR_Root->topOffDiagRowIdx     = createSequentialVec(0,numRows_TopOffDiag);
    HODLR_Root->topOffDiagColIdx     = createSequentialVec(0,numCols_TopOffDiag);
    HODLR_Root->bottOffDiagRowIdx    = createSequentialVec(0,numRows_BottOffDiag);
    HODLR_Root->bottOffDiagColIdx    = createSequentialVec(0,numCols_BottOffDiag);
 
  }
  
  createExactHODLR(HODLR_Root->left ,rank,result);
  createExactHODLR(HODLR_Root->right,rank,result);
 
}


void HODLR_Matrix::recLU_Factorize(){
   if (recLU_Factorized == false){
    double startTime = clock();
    recLUfactorTree.rootNode = new recLU_FactorTree::node;
    (recLUfactorTree.rootNode)->isLeaf = false; 
    Eigen::VectorXd dummyF = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2); 
    Eigen::MatrixXd dummyX = recLU_Factorize(dummyF, indexTree.rootNode, recLUfactorTree.rootNode);
    double endTime = clock();
    recLU_FactorizationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    recLU_Factorized = true;
  }

}

void HODLR_Matrix::recSM_Factorize(){
  if (recSM_Factorized == false){
    double startTime = clock();
    int tree_numLevels = indexTree.get_numLevels();
    for (int i = tree_numLevels - 1; i > 0; i--){
      std::vector<HODLR_Tree::node*> left,right; 
      recSM_Factorize(indexTree.rootNode,left,right,i);
    }
    HODLR_Tree::node* HODLR_Root = indexTree.rootNode;
    HODLR_Root->nodePerturbI = perturbI(&(HODLR_Root->topOffDiagU_SM),&(HODLR_Root->topOffDiagV),&(HODLR_Root->bottOffDiagU_SM),&(HODLR_Root->bottOffDiagV));
    double endTime = clock();
    recSM_FactorizationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    recSM_Factorized = true;
  }
}


void HODLR_Matrix::recSM_Factorize(HODLR_Tree::node* HODLR_Root,std::vector<HODLR_Tree::node*> &leftChildren, std::vector<HODLR_Tree::node*> &rightChildren,int desLevel){

  if (HODLR_Root->left->currLevel == desLevel){
    if (HODLR_Root->left->isLeaf)
      HODLR_Root->left->leafLU = Eigen::PartialPivLU<Eigen::MatrixXd>(HODLR_Root->left->leafMatrix);
    leftChildren.push_back(HODLR_Root->left);
  }else if (!HODLR_Root->left->isLeaf){
    std::vector<HODLR_Tree::node*> left,right;
    recSM_Factorize(HODLR_Root->left,left,right,desLevel);
    left.insert(left.end(),right.begin(),right.end());
    leftChildren = left;
  }
  
  if (HODLR_Root->right->currLevel == desLevel){
    if (HODLR_Root->right->isLeaf)
      HODLR_Root->right->leafLU = Eigen::PartialPivLU<Eigen::MatrixXd>(HODLR_Root->right->leafMatrix);
    rightChildren.push_back(HODLR_Root->right);
  }else if (!HODLR_Root->right->isLeaf){
    std::vector<HODLR_Tree::node*> left,right; 
    recSM_Factorize(HODLR_Root->right,left,right,desLevel);
    right.insert(right.end(),left.begin(),left.end());
    rightChildren = right;
  }
  
  if (HODLR_Root->topOffDiagU_SM.rows() == 0)
    HODLR_Root->topOffDiagU_SM = HODLR_Root->topOffDiagU;
  
  if (HODLR_Root->bottOffDiagU_SM.rows() == 0)
    HODLR_Root->bottOffDiagU_SM = HODLR_Root->bottOffDiagU;
  
  for (unsigned int i = 0; i < leftChildren.size(); i++){
    HODLR_Tree::node* currChild = leftChildren[i];
    int currChildSize   = currChild->max_i - currChild->min_i + 1;
    int minIdx = currChild->min_i - HODLR_Root->min_i;
    Eigen::MatrixXd RHS = HODLR_Root->topOffDiagU_SM.block(minIdx,0,currChildSize,HODLR_Root->topOffDiagU_SM.cols());
    if (currChild->isLeaf){
      HODLR_Root->topOffDiagU_SM.block(minIdx,0,currChildSize,HODLR_Root->topOffDiagU_SM.cols()) = currChild->leafLU.solve(RHS);
    }else{
      currChild->nodePerturbI = perturbI(&(currChild->topOffDiagU_SM),&(currChild->topOffDiagV),&(currChild->bottOffDiagU_SM),&(currChild->bottOffDiagV));
      HODLR_Root->topOffDiagU_SM.block(minIdx,0,currChildSize,HODLR_Root->topOffDiagU_SM.cols()) = currChild->nodePerturbI.solve(RHS);
    }
  }
  
  for (unsigned int i = 0; i < rightChildren.size(); i++){
    HODLR_Tree::node* currChild = rightChildren[i];
    int currChildSize   = currChild->max_i - currChild->min_i + 1;
    int minIdx = currChild->min_i - HODLR_Root->splitIndex_i - 1;
    Eigen::MatrixXd RHS = HODLR_Root->bottOffDiagU_SM.block(minIdx,0,currChildSize,HODLR_Root->bottOffDiagU_SM.cols());
    if (currChild->isLeaf){
      HODLR_Root->bottOffDiagU_SM.block(minIdx,0,currChildSize,HODLR_Root->bottOffDiagU_SM.cols()) = currChild->leafLU.solve(RHS);
    }else{
      currChild->nodePerturbI = perturbI(&(currChild->topOffDiagU_SM),&(currChild->topOffDiagV),&(currChild->bottOffDiagU_SM),&(currChild->bottOffDiagV));
      HODLR_Root->bottOffDiagU_SM.block(minIdx,0,currChildSize,HODLR_Root->bottOffDiagU_SM.cols()) = currChild->nodePerturbI.solve(RHS);
    }
  }
}

void HODLR_Matrix::recSM_Solve(HODLR_Tree::node* HODLR_Root,Eigen::MatrixXd &RHS){
  
  int currBlockSize =  HODLR_Root->max_i - HODLR_Root->min_i + 1;
  
  if (HODLR_Root->isLeaf == true){
    Eigen::MatrixXd currRHS = RHS.block(HODLR_Root->min_i,0,currBlockSize,RHS.cols());
    RHS.block(HODLR_Root->min_i,0,currBlockSize,RHS.cols()) = HODLR_Root->leafLU.solve(currRHS);
    return;
  }
  
  recSM_Solve(HODLR_Root->left,RHS);
  recSM_Solve(HODLR_Root->right,RHS);
  
  Eigen::MatrixXd currRHS = RHS.block(HODLR_Root->min_i,0,currBlockSize,RHS.cols());
  RHS.block(HODLR_Root->min_i,0,currBlockSize,RHS.cols()) = HODLR_Root->nodePerturbI.solve(currRHS);
  
}





Eigen::MatrixXd HODLR_Matrix::recLU_Factorize(const Eigen::MatrixXd & input_RHS,const HODLR_Tree::node* HODLR_Root, recLU_FactorTree::node* factorRoot){

  // Base case
  if (HODLR_Root->isLeaf == true){
    double startTime = clock();
    factorRoot->isLeaf = true;
    factorRoot->left   = NULL;
    factorRoot->right  = NULL;  
    Eigen::MatrixXd result;
    if (isLeafConst == true){
      if (constLeafFactorized == false){
	constLeafLU = Eigen::PartialPivLU<Eigen::MatrixXd>(constLeaf);
	constLeafFactorized = true;
      }
      result = constLeafLU.solve(input_RHS);
    }else{
      factorRoot->LU = Eigen::PartialPivLU<Eigen::MatrixXd>(HODLR_Root->leafMatrix);
      result = factorRoot->LU.solve(input_RHS); 
    }
    double endTime = clock();
    recLU_FactorLevelTimeVec[HODLR_Root->currLevel] += (endTime - startTime)/CLOCKS_PER_SEC;
    return result;
  }

  int topDiagSize = HODLR_Root->splitIndex_i - HODLR_Root->min_i + 1;
  int bottDiagSize = HODLR_Root->max_i - HODLR_Root->splitIndex_i;
  int parentRHS_Cols = input_RHS.cols();
  
  // Factorize and Solve for top diagonal matrix
  int topOffDiag_LR_Cols = HODLR_Root->topOffDiagU.cols();
  int topDiagRHS_Cols = parentRHS_Cols + topOffDiag_LR_Cols;
  Eigen::MatrixXd topDiagRHS = Eigen::MatrixXd::Zero(topDiagSize,topDiagRHS_Cols);
  topDiagRHS.leftCols(parentRHS_Cols) = input_RHS.topRows(topDiagSize);
  topDiagRHS.rightCols(topOffDiag_LR_Cols) = HODLR_Root->topOffDiagU;  
  
  recLU_FactorTree::node* leftFactorRoot = new recLU_FactorTree::node;
  leftFactorRoot->isLeaf = false;
  factorRoot->left = leftFactorRoot;
  Eigen::MatrixXd topDiagSoln = recLU_Factorize(topDiagRHS,HODLR_Root->left,leftFactorRoot);
  Eigen::MatrixXd topDiagSoln_pRHS = topDiagSoln.leftCols(parentRHS_Cols);
  factorRoot->topDiagSoln_LR = topDiagSoln.rightCols(topOffDiag_LR_Cols);
  
  
  // Factorize and Solve for bottom diagonal matrix
  int bottOffDiag_LR_Cols = HODLR_Root->bottOffDiagU.cols();
  int bottDiagRHS_Cols = parentRHS_Cols + bottOffDiag_LR_Cols;
  Eigen::MatrixXd bottDiagRHS = Eigen::MatrixXd::Zero(bottDiagSize,bottDiagRHS_Cols);
  bottDiagRHS.leftCols(parentRHS_Cols) = input_RHS.bottomRows(bottDiagSize);
  bottDiagRHS.rightCols(bottOffDiag_LR_Cols) = HODLR_Root->bottOffDiagU;
 
  recLU_FactorTree::node* rightFactorRoot = new recLU_FactorTree::node;
  rightFactorRoot->isLeaf = false;
  factorRoot->right = rightFactorRoot;
  Eigen::MatrixXd bottDiagSoln = recLU_Factorize(bottDiagRHS,HODLR_Root->right,rightFactorRoot);
  Eigen::MatrixXd bottDiagSoln_pRHS = bottDiagSoln.leftCols(parentRHS_Cols);
  factorRoot->bottDiagSoln_LR = bottDiagSoln.rightCols(bottOffDiag_LR_Cols);
 

  // Update the remaining of the matrix(generate Schur complement (S matrix));
  double startTime = clock();
  int Sdim = HODLR_Root->topOffDiagRank + HODLR_Root->bottOffDiagRank;  
  Eigen::MatrixXd S  = Eigen::MatrixXd::Zero(Sdim,Sdim);

  Eigen::MatrixXd IC = Eigen::MatrixXd::Identity(HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiagRank);
  Eigen::MatrixXd IB = Eigen::MatrixXd::Identity(HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiagRank);
  
  S.topLeftCorner(IC.rows(),IC.cols()) = IC;
  S.bottomRightCorner(IB.rows(),IB.cols()) = IB;
  
  S.topRightCorner(HODLR_Root->bottOffDiagRank,HODLR_Root->topOffDiagRank)   = (HODLR_Root->bottOffDiagV.transpose() * factorRoot->topDiagSoln_LR);
  S.bottomLeftCorner(HODLR_Root->topOffDiagRank,HODLR_Root->bottOffDiagRank) = (HODLR_Root->topOffDiagV.transpose()  * factorRoot->bottDiagSoln_LR);

  factorRoot->LU = Eigen::PartialPivLU<Eigen::MatrixXd>(S);
 
  Eigen::MatrixXd schurRHS = Eigen::MatrixXd::Zero(Sdim,parentRHS_Cols);

  schurRHS.topRows(HODLR_Root->bottOffDiagRank)    = HODLR_Root->bottOffDiagV.transpose() * topDiagSoln_pRHS;  
  schurRHS.bottomRows(HODLR_Root->topOffDiagRank)  = HODLR_Root->topOffDiagV.transpose()  * bottDiagSoln_pRHS;  

  Eigen::MatrixXd y = factorRoot->LU.solve(schurRHS);  
  Eigen::MatrixXd result(topDiagSize + bottDiagSize,input_RHS.cols());
  
  result.topRows(topDiagSize)     = topDiagSoln_pRHS  - factorRoot->topDiagSoln_LR  * y.bottomRows(HODLR_Root->topOffDiagRank);
  result.bottomRows(bottDiagSize) = bottDiagSoln_pRHS - factorRoot->bottDiagSoln_LR * y.topRows(HODLR_Root->bottOffDiagRank);


  double endTime = clock();
  recLU_FactorLevelTimeVec[HODLR_Root->currLevel] += (endTime - startTime)/CLOCKS_PER_SEC;
  return result;
}

Eigen::MatrixXd HODLR_Matrix::recLU_Solve(const Eigen::MatrixXd & input_RHS,const HODLR_Tree::node* HODLR_Root,const recLU_FactorTree::node* factorRoot){

  // Base case
  if (HODLR_Root->isLeaf == true){
      double startTime = clock();
      Eigen::MatrixXd result;
      if (isLeafConst == true)
	result = constLeafLU.solve(input_RHS);
      else
	result = factorRoot->LU.solve(input_RHS);
      double endTime = clock();
      recLU_SolveLevelTimeVec[HODLR_Root->currLevel] += (endTime - startTime)/CLOCKS_PER_SEC;
      return result;
  }
  // Low Rank Approximation
  int calculatedRankB,calculatedRankC;
    
  calculatedRankB = HODLR_Root->topOffDiagRank;
  calculatedRankC = HODLR_Root->bottOffDiagRank;
  		
  // Factorize and Solve for top diagonal matrix
  int topDiagSize = HODLR_Root->splitIndex_i - HODLR_Root->min_i + 1;
  int parentRHS_Cols = input_RHS.cols();
  int topDiagRHS_Cols = parentRHS_Cols;
  Eigen::MatrixXd topDiagRHS = Eigen::MatrixXd::Zero(topDiagSize,topDiagRHS_Cols);
  topDiagRHS.leftCols(parentRHS_Cols) = input_RHS.topRows(topDiagSize);
  Eigen::MatrixXd topDiagSoln = recLU_Solve(topDiagRHS,HODLR_Root->left,factorRoot->left);
  Eigen::MatrixXd topDiagSoln_pRHS = topDiagSoln.leftCols(parentRHS_Cols);

  // Factorize and Solve for bottom diagonal matrix
  int bottDiagSize = HODLR_Root->max_i - HODLR_Root->splitIndex_i;
  int bottDiagRHS_Cols = parentRHS_Cols;
  Eigen::MatrixXd bottDiagRHS = Eigen::MatrixXd::Zero(bottDiagSize,bottDiagRHS_Cols);
  bottDiagRHS.leftCols(parentRHS_Cols) = input_RHS.bottomRows(bottDiagSize);
  Eigen::MatrixXd bottDiagSoln = recLU_Solve(bottDiagRHS,HODLR_Root->right,factorRoot->right);
  Eigen::MatrixXd bottDiagSoln_pRHS = bottDiagSoln.leftCols(parentRHS_Cols);

  // Print Information
  if (printLevelAccuracy){
    Eigen::MatrixXd WB,WC;
    Eigen::MatrixXd VB,VC;
    WB = HODLR_Root->topOffDiagU;
    VB = HODLR_Root->topOffDiagV;
    VC = HODLR_Root->bottOffDiagV;
    WC = HODLR_Root->bottOffDiagU;;
  
    if (matrixDataAvail){
      Eigen::MatrixXd topDiag = matrixData.block(HODLR_Root->min_i,HODLR_Root->min_j,topDiagSize,topDiagSize);
      Eigen::MatrixXd topOffDiag = matrixData.block(HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1 ,topDiagSize,bottDiagSize);
      Eigen::MatrixXd bottOffDiag = matrixData.block(HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,bottDiagSize,topDiagSize);
      Eigen::MatrixXd bottDiag = matrixData.block(HODLR_Root->splitIndex_i + 1,HODLR_Root->splitIndex_j + 1,bottDiagSize,bottDiagSize);
      std::cout<<"**********************************************************************"<<std::endl;
      std::cout<<"**********************************************************************"<<std::endl;
      std::cout<<"Current Recursion Level                     = "<<HODLR_Root->currLevel<<std::endl;
      std::cout<<"Top Diagonal Matrix min_i                   = "<<HODLR_Root->min_i<<" | "<<"Top Diagonal Matrix max_i = "<<HODLR_Root->splitIndex_i<<std::endl;
      std::cout<<"Top Diagonal Matrix min_j                   = "<<HODLR_Root->min_j<<" | "<<"Top Diagonal Matrix max_j = "<<HODLR_Root->splitIndex_j<<std::endl;
      std::cout<<"Top Diagonal Matrix Solve Rel Error         = "<<(topDiag * topDiagSoln-topDiagRHS).norm()/(topDiagRHS.norm())<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"Top Off-Diagonal Approximation Rel Error    = "<<(topOffDiag - WB * VB.transpose()).norm()/topOffDiag.norm()<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"Bottom Diagonal Matrix min_i                = "<<HODLR_Root->splitIndex_i + 1<<" | "<<"Bottom Diagonal Matrix max_i = "<<HODLR_Root->max_i<<std::endl;
      std::cout<<"Bottom Diagonal Matrix min_j                = "<<HODLR_Root->splitIndex_j + 1<<" | "<<"Bottom Diagonal Matrix max_j = "<<HODLR_Root->max_j<<std::endl;
      std::cout<<"Bottom Diagonal Matrix Solve Rel Error      = "<<(bottDiag*bottDiagSoln-bottDiagRHS).norm()/(bottDiagRHS.norm())<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"Bottom Off-Diagonal Approximation Rel Error = "<<(bottOffDiag- WC * VC.transpose()).norm()/bottOffDiag.norm()<<std::endl;
    }else{
      std::cout<<"Error! Matrix data has been deleted from memory."<<std::endl;
    }
  }

  if (printLevelRankInfo){
    if (matrixDataAvail){
      Eigen::MatrixXd topOffDiag = matrixData.block(HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1 ,topDiagSize,bottDiagSize);
      Eigen::MatrixXd bottOffDiag = matrixData.block(HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,bottDiagSize,topDiagSize);
      int actualRankTop = ::SVD_LowRankApprox(topOffDiag, LR_Tolerance);
      int actualRankBott =  ::SVD_LowRankApprox(bottOffDiag, LR_Tolerance);
      std::cout<<"*****************************************************"<<std::endl;
      std::cout<<"Current Recursion Level                     = "<<HODLR_Root->currLevel<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"Top Off-Diagonal Matrix min_i               = "<<HODLR_Root->min_i<<" | "<<" Top Off-Diagonal Matrix max_i = "<<HODLR_Root->splitIndex_i<<std::endl;
      std::cout<<"Top Off-Diagonal Matrix min_j               = "<<HODLR_Root->splitIndex_j + 1 <<" | "<<" Top Off-Diagonal Matrix max_j = "<<HODLR_Root->max_j<<std::endl;
      std::cout<<"Top off-Diagonal Matrix Calculated Rank     = "<<calculatedRankB<<std::endl;
      std::cout<<"Top off-Diagonal Matrix Actual Rank         = "<<actualRankTop<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"Bottom Off-Diagonal Matrix min_i            = "<<HODLR_Root->splitIndex_i + 1 <<" | "<<" Bottom Off-Diagonal Matrix max_i = "<<HODLR_Root->max_i<<std::endl;
      std::cout<<"Bottom Off-Diagonal Matrix min_j            = "<<HODLR_Root->min_j<<" | "<<" Bottom Off-Diagonal Matrix max_j = "<<HODLR_Root->splitIndex_j<<std::endl;
      std::cout<<"Bottom Off-Diagonal Matrix Calculated Rank  = "<<calculatedRankC<<std::endl; 
      std::cout<<"Bottom off-Diagonal Matrix Actual Rank      = "<<actualRankBott<<std::endl;
    }else{
       std::cout<<"Error! Matrix data has been deleted from memory."<<std::endl;
    }
  }
  
  double startTime = clock();
  
  int Sdim = calculatedRankB + calculatedRankC;  

  Eigen::MatrixXd schurRHS = Eigen::MatrixXd::Zero(Sdim,parentRHS_Cols);
  schurRHS.topRows(HODLR_Root->bottOffDiagRank)   = HODLR_Root->bottOffDiagV.transpose() * topDiagSoln_pRHS;  
  schurRHS.bottomRows(HODLR_Root->topOffDiagRank) = HODLR_Root->topOffDiagV.transpose()  * bottDiagSoln_pRHS;  
  
  Eigen::MatrixXd y = factorRoot->LU.solve(schurRHS);  
  Eigen::MatrixXd result(topDiagSize + bottDiagSize,input_RHS.cols());
  
  result.topRows(topDiagSize)     = topDiagSoln_pRHS  - factorRoot->topDiagSoln_LR  * y.bottomRows(HODLR_Root->topOffDiagRank);
  result.bottomRows(bottDiagSize) = bottDiagSoln_pRHS - factorRoot->bottDiagSoln_LR * y.topRows(HODLR_Root->bottOffDiagRank);
  
  double endTime = clock();
  recLU_SolveLevelTimeVec[HODLR_Root->currLevel] += (endTime - startTime)/CLOCKS_PER_SEC;
  return result;
}

Eigen::MatrixXd HODLR_Matrix::recLU_Solve(const Eigen::MatrixXd & input_RHS){
  
  assert(isSquareMatrix == true);
  assert(input_RHS.rows() == matrixSize);

  if (indexTree.rootNode == NULL){
    indexTree.set_sizeThreshold(sizeThreshold);
    indexTree.createDefaultTree(matrixSize);
    initializeInfoVecotrs(indexTree.get_numLevels());
  }

  storeLRinTree();
  recLU_Factorize();
  
  double startTime = clock();
  Eigen::MatrixXd solution = recLU_Solve(input_RHS,indexTree.rootNode,recLUfactorTree.rootNode);
  double endTime = clock();
  recLU_SolveTime = (endTime-startTime)/CLOCKS_PER_SEC;
  recLU_TotalTime = recLU_FactorizationTime + recLU_SolveTime + LR_ComputationTime;
  if (printResultInfo){
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"Solver Type                      = recLU"<<std::endl;
    std::cout<<"Low-Rank Computation Time        = "<<LR_ComputationTime<<" seconds"<<std::endl;
    std::cout<<"Factorization Time               = "<<recLU_FactorizationTime<<" seconds"<<std::endl;
    std::cout<<"Solve Time                       = "<<recLU_SolveTime<<" seconds"<<std::endl; 
    std::cout<<"Total Solve Time                 = "<<recLU_TotalTime<<" seconds"<<std::endl;
    std::cout<<"LR Tolerance                     = "<<LR_Tolerance<<std::endl;
    if (matrixDataAvail == true)
      std::cout<<"Residual l2 Relative Error       = "<<((matrixData * solution) - input_RHS).norm()/input_RHS.norm()<<std::endl;
    else if (kernelDataAvail == true)
      std::cout<<"Residual l2 Relative Error       = "<<((kernelMatrixData * solution) - input_RHS).norm()/input_RHS.norm()<<std::endl;   
  } 

  return solution;

}

void HODLR_Matrix::recLU_Compute(){
  if (indexTree.rootNode == NULL){
    indexTree.set_sizeThreshold(sizeThreshold);
    indexTree.createDefaultTree(matrixSize);
    initializeInfoVecotrs(indexTree.get_numLevels());
  }
  storeLRinTree();
  recLU_Factorize();
 
}

Eigen::MatrixXd HODLR_Matrix::recSM_Solve(const Eigen::MatrixXd & input_RHS){
  assert(isSquareMatrix == true);
  assert(input_RHS.rows() == matrixSize);
  
  if (indexTree.rootNode == NULL){
    indexTree.set_sizeThreshold(sizeThreshold);
    indexTree.createDefaultTree(matrixSize);
    initializeInfoVecotrs(indexTree.get_numLevels());
  }

  storeLRinTree();
  recSM_Factorize();
  
  Eigen::MatrixXd solution = input_RHS;
  double startTime = clock();
  recSM_Solve(indexTree.rootNode,solution);
  double endTime = clock();
  recSM_SolveTime = (endTime-startTime)/CLOCKS_PER_SEC;
  recSM_TotalTime = recSM_FactorizationTime + recSM_SolveTime + LR_ComputationTime;
  
  if (printResultInfo){
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"Solver Type                      = recSM"<<std::endl;
    std::cout<<"Low-Rank Computation Time        = "<<LR_ComputationTime<<" seconds"<<std::endl;
    std::cout<<"Factorization Time               = "<<recSM_FactorizationTime<<" seconds"<<std::endl;
    std::cout<<"Solve Time                       = "<<recSM_SolveTime<<" seconds"<<std::endl; 
    std::cout<<"Total Solve Time                 = "<<recSM_TotalTime<<" seconds"<<std::endl;
    std::cout<<"LR Tolerance                     = "<<LR_Tolerance<<std::endl;
    if (matrixDataAvail == true)
      std::cout<<"Residual l2 Relative Error       = "<<((matrixData * solution) - input_RHS).norm()/input_RHS.norm()<<std::endl;
    else if (kernelDataAvail == true)
      std::cout<<"Residual l2 Relative Error       = "<<((kernelMatrixData * solution) - input_RHS).norm()/input_RHS.norm()<<std::endl;   
  } 
  
  return solution;

}

void HODLR_Matrix::findNodesAtLevel(HODLR_Tree::node* HODLR_Root, const int level,std::vector<HODLR_Tree::node*> & outputVector){

  // Base case
  
  if (HODLR_Root->isLeaf == true){
  return;
  }

  if (HODLR_Root->currLevel == level){
    outputVector.push_back(HODLR_Root);
    return;
  }

  findNodesAtLevel(HODLR_Root->left,level,outputVector);
  findNodesAtLevel(HODLR_Root->right,level,outputVector);
}

void HODLR_Matrix::findLeafNodes(HODLR_Tree::node* HODLR_Root, std::vector<HODLR_Tree::node*>& outputVector){
  // Base case
  
  if (HODLR_Root->isLeaf == true){
    outputVector.push_back(HODLR_Root);
    return;
  }
  
    findLeafNodes(HODLR_Root->left,outputVector);
    findLeafNodes(HODLR_Root->right,outputVector);  
}

int HODLR_Matrix::sumRanks(HODLR_Tree::node* HODLR_Root) {
  // Base case
  
  if (HODLR_Root->isLeaf == true){
    return 0;
  }

  int currOffDiagRank = HODLR_Root->topOffDiagRank + HODLR_Root->bottOffDiagRank;
  return currOffDiagRank + sumRanks(HODLR_Root->left) + sumRanks(HODLR_Root->right);
  
}


void HODLR_Matrix::insertIdentityIntoSpMatrix(std::vector<Eigen::Triplet<double,int> > & Sp_TripletVec, const int start_i,const  int start_j, const int input_MatrixSize, const int constant){
  
  assert(input_MatrixSize > 0);
  for (int i = 0; i < input_MatrixSize; i++){
    Eigen::Triplet<double,int> entryTriplet(start_i + i, start_j + i, constant);
    Sp_TripletVec.push_back(entryTriplet);
  }
  
}

void HODLR_Matrix::insertDenseBlockIntoSpMatrix(std::vector<Eigen::Triplet<double,int> > & Sp_TripletVec, const Eigen::MatrixXd & denseMatrix, const int start_i, const int start_j){
  
  int numRows = denseMatrix.rows();
  int numCols = denseMatrix.cols();
  for (int i = 0; i < numRows; i++)
    for (int j = 0; j < numCols; j++){
      Eigen::Triplet<double,int> entryTriplet(start_i + i, start_j + j, denseMatrix(i,j));
      Sp_TripletVec.push_back(entryTriplet);
    }
}

Eigen::SparseMatrix<double>  HODLR_Matrix::assembleExtendedSPMatrix(){
  
  int rankSum = sumRanks(indexTree.rootNode);
  int spMatrixSize = matrixData.rows() + 2 * rankSum;
  
  Eigen::SparseMatrix<double> extendedSp_Matrix(spMatrixSize, spMatrixSize);
  std::vector<Eigen::Triplet<double,int> > extendedSp_TripletVec;
  int currLevelStartIndex = spMatrixSize - 1;
 
  int tree_numLevels = indexTree.get_numLevels();
  for (int i = 0; i < tree_numLevels; i++){
    
    std::vector<HODLR_Tree::node*> currLevelNodesVec;
    findNodesAtLevel(indexTree.rootNode,i,currLevelNodesVec);
    //currLevelNodesVec = indexTree.nodeLevelVec[i];
    // Add ranks at current level
    int levelRankSum = 0;
    
    for (unsigned int j = 0; j < currLevelNodesVec.size(); j++)
      levelRankSum += currLevelNodesVec[j]->topOffDiagRank + currLevelNodesVec[j]->bottOffDiagRank;
    
    int block_StartIndex = currLevelStartIndex - 2 * levelRankSum + 1;
    
    // Assemble low rank factors to sparse matrix
    
    // Construct K, U and V subblocks
   
    int local_IstartIndex_i = levelRankSum;
    int local_IstartIndex_j = 0;
    int local_KstartIndex   = levelRankSum;
    int local_UstartIndex_i = 0;
    int local_UstartIndex_j = 0;
    int local_VstartIndex_i = 0;
    int local_VstartIndex_j = 0;
    
    for (unsigned int j = 0; j < currLevelNodesVec.size(); j++){
    
      
      int offset;
      if (j == 0)
	offset = currLevelNodesVec[j]->min_i;
      else
	offset = currLevelNodesVec[j]->min_i - currLevelNodesVec[j-1]->max_i - 1;
   
      local_UstartIndex_i += offset;
      local_VstartIndex_j += offset;
     
      int topOffDiagRank    = currLevelNodesVec[j]->topOffDiagRank;
      int bottOffDiagRank   = currLevelNodesVec[j]->bottOffDiagRank;
            
      Eigen::MatrixXd topOffDiagK  = Eigen::MatrixXd::Identity(topOffDiagRank,topOffDiagRank);
      Eigen::MatrixXd bottOffDiagK = Eigen::MatrixXd::Identity(bottOffDiagRank,bottOffDiagRank);
    
      Eigen::MatrixXd topOffDiagU  = currLevelNodesVec[j]->topOffDiagU;
      Eigen::MatrixXd bottOffDiagU = currLevelNodesVec[j]->bottOffDiagU;
      Eigen::MatrixXd topOffDiagV  = currLevelNodesVec[j]->topOffDiagV;
      Eigen::MatrixXd bottOffDiagV = currLevelNodesVec[j]->bottOffDiagV;
  
      // Insert Is
    
      int global_IstartIndex_i = block_StartIndex + local_IstartIndex_i;
      int global_IstartIndex_j = block_StartIndex + local_IstartIndex_j;

      insertIdentityIntoSpMatrix(extendedSp_TripletVec,global_IstartIndex_i,global_IstartIndex_j,topOffDiagRank + bottOffDiagRank,-1);
      insertIdentityIntoSpMatrix(extendedSp_TripletVec,global_IstartIndex_j,global_IstartIndex_i,topOffDiagRank + bottOffDiagRank,-1);

      local_IstartIndex_i += topOffDiagRank + bottOffDiagRank;
      local_IstartIndex_j += topOffDiagRank + bottOffDiagRank;
      
      // Insert Ks
      int global_KstartIndex = block_StartIndex + local_KstartIndex;
      
      insertDenseBlockIntoSpMatrix(extendedSp_TripletVec,topOffDiagK,global_KstartIndex,global_KstartIndex + bottOffDiagRank);
      insertDenseBlockIntoSpMatrix(extendedSp_TripletVec,bottOffDiagK,global_KstartIndex + topOffDiagRank,global_KstartIndex);

      local_KstartIndex += topOffDiagRank + bottOffDiagRank;
      
      // Insert Us
      int global_UstartIndex_i = local_UstartIndex_i;
      int global_UstartIndex_j = block_StartIndex + local_UstartIndex_j;
      
       insertDenseBlockIntoSpMatrix(extendedSp_TripletVec,topOffDiagU,global_UstartIndex_i,global_UstartIndex_j);
       
      local_UstartIndex_i += topOffDiagU.rows();
      local_UstartIndex_j += topOffDiagRank;
      
      global_UstartIndex_i = local_UstartIndex_i;
      global_UstartIndex_j = block_StartIndex + local_UstartIndex_j;
   
      insertDenseBlockIntoSpMatrix(extendedSp_TripletVec,bottOffDiagU,global_UstartIndex_i,global_UstartIndex_j);

      local_UstartIndex_i += bottOffDiagU.rows();
      local_UstartIndex_j += bottOffDiagRank;
     
      // Insert Vs
      int global_VstartIndex_i = block_StartIndex + local_VstartIndex_i;
      int global_VstartIndex_j = local_VstartIndex_j;
      
      insertDenseBlockIntoSpMatrix(extendedSp_TripletVec,bottOffDiagV.transpose(),global_VstartIndex_i,global_VstartIndex_j);
      
      local_VstartIndex_i += bottOffDiagRank;
      local_VstartIndex_j += bottOffDiagV.rows();
      global_VstartIndex_i = block_StartIndex + local_VstartIndex_i;
      global_VstartIndex_j = local_VstartIndex_j; 
     
      insertDenseBlockIntoSpMatrix(extendedSp_TripletVec,topOffDiagV.transpose(),global_VstartIndex_i,global_VstartIndex_j);

      local_VstartIndex_i += topOffDiagRank;
      local_VstartIndex_j += topOffDiagV.rows();
      
    }

    currLevelStartIndex -= 2*levelRankSum;
			      
  }

  // Insert leaf nodes
  int leaf_assemblyIndex = 0;
  std::vector<HODLR_Tree::node*> leafNodesVec;
  findLeafNodes(indexTree.rootNode,leafNodesVec);
  //leafNodesVec = indexTree.leafNodesVec;
  for (unsigned int i = 0; i < leafNodesVec.size(); i++){
    int leafSize = leafNodesVec[i]->max_i - leafNodesVec[i]->min_i + 1;
    //Eigen::MatrixXd leafMatrix = matrixData.block(leafNodesVec[i]->min_i,leafNodesVec[i]->min_j,leafSize,leafSize);
    insertDenseBlockIntoSpMatrix(extendedSp_TripletVec,leafNodesVec[i]->leafMatrix,leaf_assemblyIndex,leaf_assemblyIndex);
    leaf_assemblyIndex += leafSize;
  }
  extendedSp_Matrix.setFromTriplets(extendedSp_TripletVec.begin(),extendedSp_TripletVec.end());   
  return extendedSp_Matrix;

}

Eigen::MatrixXd HODLR_Matrix::extendedSp_Solve(const Eigen::MatrixXd & input_RHS){

  assert(isSquareMatrix);
  assert(input_RHS.rows() == matrixSize);
  if (indexTree.rootNode == NULL){
    indexTree.set_sizeThreshold(sizeThreshold);
    indexTree.createDefaultTree(matrixSize);
    initializeInfoVecotrs(indexTree.get_numLevels());
  }
  
  storeLRinTree();
    
  if (assembled_ExtendedSp == false){
    double startTime = clock();
    Eigen::SparseMatrix<double> extendedSp_Matrix = assembleExtendedSPMatrix();
    extendedSp_Matrix.makeCompressed();
    double endTime = clock();
    if (saveExtendedSp_Matrix == true){
      //Eigen::MatrixXd sp(extendedSp_Matrix);
      //saveMatrixXdToBinary(sp,extendedSp_SavePath);
      // To Do
    }
    
    extendedSp_AssemblyTime = (endTime-startTime)/CLOCKS_PER_SEC;
    extendedSp_Size = extendedSp_Matrix.rows();
    startTime = clock();
    extendedSp_Solver.compute(extendedSp_Matrix);
    endTime = clock();
    if (extendedSp_Solver.info() != Eigen::Success){
      std::cout<<"Extended sparse matrix factorization failed."<<std::endl;
      exit(EXIT_FAILURE);
    }
    extendedSp_FactorizationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    assembled_ExtendedSp = true;
  }
  

  // Construct Extended RHS
  Eigen::MatrixXd sp_Solution;
  Eigen::MatrixXd extendedSp_RHS = Eigen::MatrixXd::Zero(extendedSp_Size,input_RHS.cols());
  extendedSp_RHS.topLeftCorner(input_RHS.rows(),input_RHS.cols()) = input_RHS;
  double startTime = clock();
  sp_Solution = extendedSp_Solver.solve(extendedSp_RHS);
  double endTime = clock();
  extendedSp_SolveTime = (endTime-startTime)/CLOCKS_PER_SEC;
  if (extendedSp_Solver.info() != Eigen::Success){
    std::cout<<"Extended sparse matrix solve failed."<<std::endl;
    exit(EXIT_FAILURE);
  }

  // Extract Dense Solution
  Eigen::MatrixXd solution = sp_Solution.topLeftCorner(input_RHS.rows(),input_RHS.cols());
  extendedSp_TotalTime = LR_ComputationTime + extendedSp_AssemblyTime + extendedSp_FactorizationTime + extendedSp_SolveTime;
  if (printResultInfo){
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"Solver Type                      = extendedSp"<<std::endl;
    std::cout<<"Low-Rank Computation Time        = "<<LR_ComputationTime<<" seconds"<<std::endl;
    std::cout<<"Assembly Time                    = "<<extendedSp_AssemblyTime<<" seconds"<<std::endl;
    std::cout<<"Factorization Time               = "<<extendedSp_FactorizationTime<<" seconds"<<std::
endl;
    std::cout<<"Solve Time                       = "<<extendedSp_SolveTime<<" seconds"<<std::endl; 
    std::cout<<"Total Solve Time                 = "<<extendedSp_TotalTime<<" seconds"<<std::endl;
    std::cout<<"LR Tolerance                     = "<<LR_Tolerance<<std::endl;
    std::cout<<"Residual l2 Relative Error       = "<<((matrixData * solution)-input_RHS).norm()/input_RHS.norm()<<std::endl;
  }
  
  return solution;

}

Eigen::MatrixXd HODLR_Matrix::oneStep_Iterate(const Eigen::MatrixXd & prevStep_result,const Eigen::MatrixXd & RHS, const Eigen::MatrixXd & initSolveGuess, Eigen::MatrixXd & prevStep_Product,const std::string directSolve_Method){

  if (matrixDataAvail == true)
    prevStep_Product = matrixData * prevStep_result;
  else if (kernelDataAvail == true)
    prevStep_Product = kernelMatrixData * prevStep_result;
  else{
    std::cout<<"Error! No matrix data available."<<std::endl;
    exit(EXIT_FAILURE);
  }

  if (directSolve_Method == "recLU"){
    Eigen::MatrixXd update = prevStep_result - recLU_Solve(prevStep_Product);
    return initSolveGuess + update;
  }

  if (directSolve_Method == "extendedSp"){
    Eigen::MatrixXd update = prevStep_result - extendedSp_Solve(prevStep_Product);
    return initSolveGuess + update;
  }

  std::cout<<"Error! No such direct solver type."<<std::endl;
  exit(EXIT_FAILURE);
  
}

Eigen::MatrixXd HODLR_Matrix::iterative_Solve(const Eigen::MatrixXd & input_RHS, const int maxIterations, const double stop_tolerance,const double init_LRTolerance,const std::string input_LRMethod, const std::string directSolve_Method = "recLU"){
    
  assert(input_RHS.rows() == matrixSize);
  // double prev_LRTolerance = LR_Tolerance;
  set_LRTolerance(init_LRTolerance);

  std::string prev_LRMethod = indexTree.get_def_LRMethod(); 
  set_LRMethod(input_LRMethod);

  bool save_printResultInfo = printResultInfo;
  printResultInfo = false;

  if (indexTree.rootNode == NULL){
    indexTree.set_sizeThreshold(sizeThreshold);
    indexTree.createDefaultTree(matrixSize);
    initializeInfoVecotrs(indexTree.get_numLevels());
  }
  
  storeLRinTree();

  double startTime = clock();
  Eigen::MatrixXd init_Guess;
  if (directSolve_Method == "recLU")
    init_Guess = recLU_Solve(input_RHS);
  else if (directSolve_Method == "extendedSp")
    init_Guess = extendedSp_Solve(input_RHS);
  else{
    std::cout<<"Error! Unknown direct solver type!"<<std::endl;
    exit(EXIT_FAILURE);
  }
  Eigen::MatrixXd currStep_Soln = init_Guess;
  Eigen::MatrixXd nextStep_Soln;
  Eigen::MatrixXd currStep_Product;
  int num_Iter = 1;
  double tolerance = 1;
  double iterSolveTime = recLU_SolveTime;
  while (tolerance > stop_tolerance){
    iter_IterTimeVec.push_back(iterSolveTime);
    double iterStartTime = clock();
    nextStep_Soln = oneStep_Iterate(currStep_Soln,input_RHS,init_Guess,currStep_Product,directSolve_Method);
    //tolerance = (nextStep_Soln - currStep_Soln).norm()/currStep_Soln.norm();
    tolerance = (currStep_Product - input_RHS).norm()/input_RHS.norm();
    double iterEndTime = clock();
    iterSolveTime = (iterEndTime - iterStartTime)/CLOCKS_PER_SEC;
    iter_IterAccuracyVec.push_back(tolerance);
    currStep_Soln = nextStep_Soln;      
    std::cout<<num_Iter<<" "<<tolerance<<std::endl;
    num_Iter ++;
    if (num_Iter > maxIterations)
      break;
  }
  Eigen::MatrixXd solution = currStep_Soln;
  double endTime = clock();
  totalIter_SolveTime = (endTime-startTime)/CLOCKS_PER_SEC;
  printResultInfo = save_printResultInfo;
  if (printResultInfo){
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"Solver Type                        = iterative"<<std::endl;
    std::cout<<"Low-Rank Computation Time          = "<<LR_ComputationTime<<" seconds"<<std::endl;
    std::cout<<"HODLR Factorization Time           = "<<recLU_FactorizationTime<<" seconds"<<std::endl;
    std::cout<<"HODLR Direct Solver Total Time     = "<<recLU_FactorizationTime + LR_ComputationTime + recLU_SolveTime<<" seconds"<<std::endl; 
    std::cout<<"Total Iteration Time               = "<<totalIter_SolveTime<<" seconds"<<std::endl;
    std::cout<<"Total Solve Time                   = "<<totalIter_SolveTime + LR_ComputationTime<<" seconds"<<std::endl;
    std::cout<<"LR Tolerance                       = "<<LR_Tolerance<<std::endl;
    std::cout<<"Number of Iterations               = "<<num_Iter<<std::endl;
    if (matrixDataAvail == true)
      std::cout<<"Residual l2 Relative Error         = "<<((matrixData * solution) - input_RHS).norm()/input_RHS.norm()<<std::endl;
    else if (kernelDataAvail == true)
      std::cout<<"Residual l2 Relative Error         = "<<((kernelMatrixData * solution) - input_RHS).norm()/input_RHS.norm()<<std::endl;  
  }
  
  // restore previous state;
  //set_LRTolerance(prev_LRTolerance);
  //set_LRMethod(prev_LRMethod);
  return solution;
}
  
void HODLR_Matrix::saveSolverInfo(const std::string outputFileName){
  std::string LR_LevelTiming          = outputFileName + "LR_Timing";
  std::string recLU_FactorLevelTiming = outputFileName + "recLU_FactorTiming";
  std::string recLU_SolveLevelTiming  = outputFileName + "recLU_SolveTiming";
  std::string iter_IterTiming         = outputFileName + "iter_IterTiming";
  std::string iter_Accuracy           = outputFileName + "iter_Accuracy";
  std::string levelRankAverage        = outputFileName + "levelRankAverage";

  saveVectorAsText(LR_LevelTiming,LR_ComputationLevelTimeVec);
  saveVectorAsText(recLU_FactorLevelTiming,recLU_FactorLevelTimeVec);
  saveVectorAsText(recLU_SolveLevelTiming,recLU_SolveLevelTimeVec);
  saveVectorAsText(iter_IterTiming,iter_IterTimeVec);
  saveVectorAsText(iter_Accuracy,iter_IterAccuracyVec);
  saveVectorAsText(levelRankAverage,levelRankAverageVec);

}

void HODLR_Matrix::reset_attributes(){
  LRStoredInTree       = false;
  recLU_Factorized     = false;
  assembled_ExtendedSp = false;

  recLU_FactorizationTime = 0;
  recLU_SolveTime = 0;
  LR_ComputationTime = 0;
  extendedSp_AssemblyTime = 0;
  extendedSp_FactorizationTime = 0;
  extendedSp_SolveTime = 0;
  extendedSp_Size = 0;
  return;
}

void HODLR_Matrix::set_LRTolerance(double input_tolerance){
  if ((matrixDataAvail == false) && (matrixDataAvail_Sp == false) && (kernelDataAvail == false)){
    std::cout<<"Error! Matrix data has been deleted from memory!"<<std::endl;
    exit(EXIT_FAILURE);
  }
  if (LR_Tolerance != input_tolerance){
    LR_Tolerance = input_tolerance;
    reset_attributes();
  }
}  

void HODLR_Matrix::set_minPivot(double input_minPivot){
  if ((matrixDataAvail == false) && (matrixDataAvail_Sp == false) && (kernelDataAvail == false)){
    std::cout<<"Error! Matrix data has been deleted from memory!"<<std::endl;
    exit(EXIT_FAILURE);
  }
  if (minPivot != input_minPivot){
    minPivot = input_minPivot;
    reset_attributes();
  }
}

void HODLR_Matrix::set_LRMethod(std::string input_LRMethod){
  if ((matrixDataAvail == false) && (matrixDataAvail_Sp == false) && (kernelDataAvail == false)){
    std::cout<<"Error! Matrix data has been deleted from memory!"<<std::endl;
    exit(EXIT_FAILURE);
  }
  if (indexTree.get_def_LRMethod() != input_LRMethod){
    indexTree.set_LRMethod(input_LRMethod);
    reset_attributes();
  }
}

void HODLR_Matrix::set_FreeMatrixMemory(bool inputVal){
  freeMatrixMemory = inputVal;
  freeMatrixMemory_Sp = inputVal;
}

void HODLR_Matrix::set_BoundaryDepth(int inputBoundaryDepth){
  boundaryDepth = inputBoundaryDepth;
}

void HODLR_Matrix::set_LeafConst(){
  isLeafConst = true;
}

void HODLR_Matrix::saveExtendedSp(std::string savePath){
  saveExtendedSp_Matrix = true;
  extendedSp_SavePath = savePath;
}

void HODLR_Matrix::freeDenseMatMem(){
  matrixData.resize(0,0);
  matrixDataAvail = false;
}

void HODLR_Matrix::freeSparseMatMem(){
  matrixData_Sp.resize(0,0);
  matrixDataAvail_Sp = false;
}

void HODLR_Matrix::freeMatrixData(){
  freeDenseMatMem();
  freeSparseMatMem();
}

void HODLR_Matrix::destroyAllData(){
  freeMatrixData();
  if (indexTree.rootNode != NULL)
    indexTree.freeTree(indexTree.rootNode);
  indexTree.rootNode = NULL;
}
void HODLR_Matrix::recalculateSize(){
  if (indexTree.rootNode != NULL){
    matrixNumRows = indexTree.rootNode->max_i - indexTree.rootNode->min_i + 1;
    matrixNumCols = indexTree.rootNode->max_j - indexTree.rootNode->min_j + 1;
    matrixSize    = matrixNumRows;
  }
}

void HODLR_Matrix::correctIndices(){
  if (indexTree.rootNode != NULL)
    indexTree.correctIndices();
}

double HODLR_Matrix::get_recLU_FactorizationTime() const{
  return recLU_FactorizationTime;
}

double HODLR_Matrix::get_recLU_SolveTime() const{
  return recLU_SolveTime;
}

double HODLR_Matrix::get_recLU_TotalTime() const{
  return recLU_TotalTime;
}

double HODLR_Matrix::get_extendedSp_AssemblyTime() const{
  return extendedSp_AssemblyTime;
}

double HODLR_Matrix::get_extendedSp_FactorizationTime() const{
  return extendedSp_FactorizationTime;
}

double HODLR_Matrix::get_extendedSp_SolveTime() const{
  return extendedSp_SolveTime;
}

double HODLR_Matrix::get_extendedSp_TotalTime() const{
  return extendedSp_TotalTime;
}

double HODLR_Matrix::get_LR_ComputationTime() const{
  return LR_ComputationTime;
}

double HODLR_Matrix::get_totalIter_SolveTime() const {
  return totalIter_SolveTime;
}
double HODLR_Matrix::get_MatrixSize() const{
  return matrixSize;
}

int HODLR_Matrix::rows() const {
  return matrixNumRows;
}

int HODLR_Matrix::cols() const {
  return matrixNumCols;
}

double HODLR_Matrix:: norm(){
  return block(0,0,rows(),cols()).norm();
}
/*
HODLR_Tree::node* HODLR_Matrix::get_TreeRootNode(){
  return indexTree.rootNode;
}
*/
/************************************ Acessing HODLR Entries *******************************/
Eigen::MatrixXd HODLR_Matrix::row(int row){
  return block(row,0,1,matrixNumCols);
}

Eigen::MatrixXd HODLR_Matrix::col(int col){
  return block(0,col,matrixNumRows,1);
}

Eigen::MatrixXd HODLR_Matrix::block(int min_i,int min_j,int numRows,int numCols){
  
  if (matrixDataAvail == true){
    return matrixData.block(min_i,min_j,numRows,numCols);
  }else  if(matrixDataAvail_Sp == true){
    return Eigen::MatrixXd(matrixData_Sp.block(min_i,min_j,numRows,numCols));    
  }else if (kernelDataAvail == true){
    return kernelMatrixData.block(min_i,min_j,numRows,numCols);
  }else if (LRStoredInTree == true){
    int max_i = min_i + numRows - 1;
    int max_j = min_j + numCols - 1;
    Eigen::MatrixXd blkMatrix = Eigen::MatrixXd::Zero(numRows,numCols);
    fill_Block(blkMatrix,indexTree.rootNode,min_i,min_j,max_i,max_j);
    return blkMatrix;
  }else{
    std::cout<<"Error! Cannot retrieve matrix information."<<std::endl; 
    exit(EXIT_FAILURE);
  }
}

void HODLR_Matrix::fill_Block(Eigen::MatrixXd & blkMatrix,HODLR_Tree::node* root,int min_i,int min_j,int max_i,int max_j){
  
  if ((max_i < root->min_i) || (max_j < root->min_j) || (min_i > root->max_i) || (min_j > root->max_j ))
    return;
  // Find parts corresponding to off-diagonal blocks
  if (root->isLeaf == true){
    int leaf_Min_i  = std::max(root->min_i,min_i);
    int leaf_Min_j  = std::max(root->min_j,min_j);
    int leaf_Max_i  = std::min(root->max_i,max_i);
    int leaf_Max_j  = std::min(root->max_j,max_j);
    int blk_Min_i   = leaf_Min_i - min_i;
    int blk_Min_j   = leaf_Min_j - min_j;
    int blk_numRows = leaf_Max_i - leaf_Min_i + 1;
    int blk_numCols = leaf_Max_j - leaf_Min_j + 1;
    blkMatrix.block(blk_Min_i,blk_Min_j,blk_numRows,blk_numCols) = root->leafMatrix.block(leaf_Min_i - root->min_i,leaf_Min_j - root->min_j,blk_numRows,blk_numCols);
    return;
  }
  
  //Top off diagonal
  if ((max_j >= root->splitIndex_j + 1) && (min_i <= root->splitIndex_i)){
    int topOffDiag_Min_i = std::max(root->min_i,min_i);
    int topOffDiag_Min_j = std::max(root->splitIndex_j + 1,min_j);
    int topOffDiag_Max_i = std::min(root->splitIndex_i,max_i);
    int topOffDiag_Max_j = std::min(root->max_j,max_j);
    int LR_Min_i   = topOffDiag_Min_i - root->min_i;
    int LR_Min_j   = topOffDiag_Min_j - (root->splitIndex_j + 1);
    int LR_numRows = topOffDiag_Max_i - topOffDiag_Min_i + 1;
    int LR_numCols = topOffDiag_Max_j - topOffDiag_Min_j + 1;
    int blk_Min_i  = topOffDiag_Min_i - min_i;
    int blk_Min_j  = topOffDiag_Min_j - min_j;
    fill_BlockWithLRProduct(blkMatrix,LR_Min_i,LR_Min_j,LR_numRows,LR_numCols,root->topOffDiagU,root->topOffDiagV,blk_Min_i,blk_Min_j);
    
  }
  //Bottom off diagonal
  if ((min_j <= root->splitIndex_j) && (max_i >= root->splitIndex_i + 1)){
    int bottOffDiag_Min_i = std::max(root->splitIndex_i + 1,min_i);
    int bottOffDiag_Min_j = std::max(root->min_j,min_j);
    int bottOffDiag_Max_i = std::min(root->max_i,max_i);
    int bottOffDiag_Max_j = std::min(root->splitIndex_j,max_j);
    int LR_Min_i   = bottOffDiag_Min_i - (root->splitIndex_i + 1);
    int LR_Min_j   = bottOffDiag_Min_j - root->min_j;
    int LR_numRows = bottOffDiag_Max_i - bottOffDiag_Min_i + 1;
    int LR_numCols = bottOffDiag_Max_j - bottOffDiag_Min_j + 1;
    int blk_Min_i  = bottOffDiag_Min_i - min_i;
    int blk_Min_j  = bottOffDiag_Min_j - min_j;
    fill_BlockWithLRProduct(blkMatrix,LR_Min_i,LR_Min_j,LR_numRows,LR_numCols,root->bottOffDiagU,root->bottOffDiagV,blk_Min_i,blk_Min_j);  
  }
  
  // Find Parts corresponding to diagonal blocks
  fill_Block(blkMatrix,root->right,min_i,min_j,max_i,max_j);
  fill_Block(blkMatrix,root->left,min_i,min_j,max_i,max_j);
  
}

void HODLR_Matrix::fill_BlockWithLRProduct(Eigen::MatrixXd & blkMatrix,int LR_Min_i,int LR_Min_j, int LR_numRows, int LR_numCols,Eigen::MatrixXd & LR_U,Eigen::MatrixXd & LR_V,int blk_Min_i,int blk_Min_j){
  blkMatrix.block(blk_Min_i,blk_Min_j,LR_numRows,LR_numCols) = LR_U.block(LR_Min_i,0,LR_numRows,LR_U.cols()) * LR_V.block(LR_Min_j,0,LR_numCols,LR_V.cols()).transpose();
}


Eigen::MatrixXd& HODLR_Matrix::returnTopOffDiagU(){
  assert(indexTree.rootNode != NULL);
  assert(LRStoredInTree     == true);
  return indexTree.rootNode->topOffDiagU;
}

Eigen::MatrixXd& HODLR_Matrix::returnTopOffDiagV(){
  assert(indexTree.rootNode != NULL);
  assert(LRStoredInTree     == true);
  return indexTree.rootNode->topOffDiagV;
}

Eigen::MatrixXd& HODLR_Matrix::returnBottOffDiagU(){
  assert(indexTree.rootNode != NULL);
  assert(LRStoredInTree     == true);
  return indexTree.rootNode->bottOffDiagU;
}

Eigen::MatrixXd& HODLR_Matrix::returnBottOffDiagV(){
  assert(indexTree.rootNode != NULL);
  assert(LRStoredInTree     == true);
  return indexTree.rootNode->bottOffDiagV;
}

void splitAtTop(HODLR_Matrix& self,HODLR_Matrix& topHODLR, HODLR_Matrix& bottHODLR){
  assert(self.isSquareMatrix == true);
  self.storeLRinTree();
  if (self.indexTree.rootNode->isLeaf == true)
    return;
  topHODLR.indexTree.rootNode  = self.indexTree.rootNode->left;
  bottHODLR.indexTree.rootNode = self.indexTree.rootNode->right; 
  int topDiagRows  = self.indexTree.rootNode->splitIndex_i - self.indexTree.rootNode->min_i + 1;  
  int topDiagCols  = self.indexTree.rootNode->splitIndex_j - self.indexTree.rootNode->min_j + 1;
  int bottDiagRows = self.indexTree.rootNode->max_i - self.indexTree.rootNode->splitIndex_i;
  int bottDiagCols = self.indexTree.rootNode->max_j - self.indexTree.rootNode->splitIndex_j;
  if (self.matrixDataAvail == true){
    Eigen::MatrixXd topBlk  = self.matrixData.topLeftCorner(topDiagRows,topDiagCols);
    Eigen::MatrixXd bottBlk = self.matrixData.bottomRightCorner(bottDiagRows,bottDiagCols); 
    topHODLR.matrixData  = topBlk;
    bottHODLR.matrixData = bottBlk;
  }  
  if (self.matrixDataAvail_Sp == true){
    Eigen::SparseMatrix<double> topBlk  = self.matrixData_Sp.topLeftCorner(topDiagCols,topDiagRows);
    Eigen::SparseMatrix<double> bottBlk = self.matrixData_Sp.bottomRightCorner(bottDiagRows,bottDiagCols); 
    topHODLR.matrixData_Sp = topBlk;
    bottHODLR.matrixData_Sp = bottBlk;
  }
  
  topHODLR.correctIndices();
  bottHODLR.correctIndices();
  topHODLR.recalculateSize();
  
  bottHODLR.recalculateSize();
  topHODLR.initInfoVectors();
  bottHODLR.initInfoVectors();
  
  topHODLR.isSquareMatrix  = (topDiagRows == topDiagCols);
  bottHODLR.isSquareMatrix = (bottDiagRows == bottDiagCols);
  topHODLR.LRStoredInTree = true;
  bottHODLR.LRStoredInTree = true;
  
  delete self.indexTree.rootNode;
  self.indexTree.rootNode = NULL;
}

HODLR_Matrix  HODLR_Matrix::topDiag(){
  storeLRinTree();
  std::shared_ptr<HODLR_Matrix> topDiagMatrixPtr(new HODLR_Matrix);
  *(topDiagMatrixPtr) = *this;
  topDiagMatrixPtr->keepTopDiag();
  return *topDiagMatrixPtr;
  
}

HODLR_Matrix HODLR_Matrix::bottDiag(){
  storeLRinTree();
  std::shared_ptr<HODLR_Matrix> bottDiagMatrixPtr(new HODLR_Matrix);
  *(bottDiagMatrixPtr) = *this;
  bottDiagMatrixPtr->keepBottDiag();
  return *bottDiagMatrixPtr;
}

void HODLR_Matrix::keepTopDiag(){
  assert(indexTree.rootNode != NULL);
  if (indexTree.rootNode->isLeaf == true)
    return;
  int topDiagSize = indexTree.rootNode->splitIndex_i - indexTree.rootNode->min_i + 1;
  matrixSize = topDiagSize;
  if (matrixDataAvail == true)
    matrixData = matrixData.topLeftCorner(topDiagSize,topDiagSize);
  if (matrixDataAvail_Sp == true)
    matrixData_Sp = matrixData_Sp.topLeftCorner(topDiagSize,topDiagSize);
  HODLR_Tree::node* topDiagPtr = indexTree.rootNode->left;
  indexTree.freeTree(indexTree.rootNode->right);
  delete indexTree.rootNode;
  indexTree.rootNode = topDiagPtr;
  recalculateSize();
}

void HODLR_Matrix::keepBottDiag(){
  assert(indexTree.rootNode != NULL);
  if (indexTree.rootNode->isLeaf == true){
    freeDenseMatMem();
    freeSparseMatMem();
    indexTree.rootNode = NULL;
    std::cout<<"Warning! Matrix too small for splitting!"<<std::endl;
    return;
  }
  int bottDiagSize = indexTree.rootNode->max_i - indexTree.rootNode->splitIndex_i;
  matrixSize = bottDiagSize;
  if (matrixDataAvail == true)
    matrixData = matrixData.bottomRightCorner(bottDiagSize,bottDiagSize);
  if (matrixDataAvail_Sp == true)
    matrixData_Sp = matrixData_Sp.bottomRightCorner(bottDiagSize,bottDiagSize);
  HODLR_Tree::node* bottDiagPtr = indexTree.rootNode->right;
  indexTree.freeTree(indexTree.rootNode->left);
  delete indexTree.rootNode;
  indexTree.rootNode = bottDiagPtr;
  indexTree.correctIndices();
  recalculateSize();
}

double HODLR_Matrix::calcAbsDiff(){

  assert(LRStoredInTree == true);
  if (matrixDataAvail){
    matrixDataAvail = false;
    double result = (matrixData - block(0,0,rows(),cols())).norm();
    matrixDataAvail = true;
    return result;
  }else if (matrixDataAvail_Sp){
    matrixDataAvail_Sp = false;
    double result = (Eigen::MatrixXd(matrixData_Sp) - block(0,0,rows(),cols())).norm();
    matrixDataAvail_Sp = true;
    return result;
  }else if (kernelDataAvail){
    kernelDataAvail = false;
    double result =  (kernelMatrixData.block(0,0,rows(),cols()) - block(0,0,rows(),cols())).norm();
    kernelDataAvail = true;
    return result;
  }else{
    std::cout<<"Error! No matrix data avilable."<<std::endl;
    exit(EXIT_FAILURE);
  }
}
void HODLR_Matrix::check_Structure(){
  assert(indexTree.rootNode != NULL);
  assert(LRStoredInTree    == true);
  check_Structure(indexTree.rootNode);
}

void HODLR_Matrix::check_Structure(HODLR_Tree::node* HODLR_Root){
    assert(HODLR_Root->currLevel >= 0 );
    assert(HODLR_Root->min_i >= 0);
    assert(HODLR_Root->min_j >= 0);
    assert(HODLR_Root->min_i < HODLR_Root->max_i);
    assert(HODLR_Root->min_j < HODLR_Root->max_j);
    if (HODLR_Root->isLeaf == true){
      assert(std::isfinite(HODLR_Root->leafMatrix.norm()));
      return;
    }
    assert(HODLR_Root->splitIndex_i > HODLR_Root->min_i); 
    assert(HODLR_Root->splitIndex_i < HODLR_Root->max_i);
    assert(HODLR_Root->splitIndex_j > HODLR_Root->min_j); 
    assert(HODLR_Root->splitIndex_j < HODLR_Root->max_j);
    assert(std::isfinite(HODLR_Root->topOffDiagU.norm()));
    assert(std::isfinite(HODLR_Root->topOffDiagV.norm()));
    assert(std::isfinite(HODLR_Root->bottOffDiagU.norm()));
    assert(std::isfinite(HODLR_Root->bottOffDiagV.norm()));
    assert(HODLR_Root->topOffDiagRank > 0);
    assert(HODLR_Root->bottOffDiagRank > 0);
    check_Structure(HODLR_Root->left);
    check_Structure(HODLR_Root->right);
}



double HODLR_Matrix::determinant(){
  calcDeterminant();
  return determinant_;
}


double HODLR_Matrix::logAbsDeterminant(){
  calcDeterminant();
  return logAbsDeterminant_;
}

void HODLR_Matrix::calcDeterminant(){
  storeLRinTree();
  if (calculatedDet == false){
    determinant_ = 1;
    logAbsDeterminant_ = 0;
    recSM_Factorize();
    calcDeterminant(indexTree.rootNode);
  }
}

void HODLR_Matrix::calcDeterminant(HODLR_Tree::node* HODLR_Root){
  
  if (HODLR_Root->isLeaf == true){
    Eigen::PartialPivLU<Eigen::MatrixXd> luFact(HODLR_Root->leafMatrix);
    Eigen::MatrixXd luMatrix = luFact.matrixLU();
    for (int i = 0; i < luMatrix.rows(); i++){
      determinant_ *= luMatrix(i,i);
      logAbsDeterminant_ += log(fabs(luMatrix(i,i)));
    }
    return;
  }
  
  int numRows_TopOffDiag  = HODLR_Root->splitIndex_i - HODLR_Root->min_i + 1; 
  int numRows_BottOffDiag = HODLR_Root->max_i - HODLR_Root->splitIndex_i;
  int numCols_TopOffDiag  = HODLR_Root->max_j - HODLR_Root->splitIndex_j;
  int numCols_BottOffDiag = HODLR_Root->splitIndex_j - HODLR_Root->min_j + 1;
  
  int rankTotal = HODLR_Root->topOffDiagRank + HODLR_Root->bottOffDiagRank;
  int blockSize = numRows_TopOffDiag + numRows_BottOffDiag;
  
  Eigen::MatrixXd U  = Eigen::MatrixXd::Zero(blockSize,rankTotal);
  Eigen::MatrixXd VT = Eigen::MatrixXd::Zero(rankTotal,blockSize);
  

  U.topLeftCorner(numRows_TopOffDiag,HODLR_Root->topOffDiagRank) = HODLR_Root->topOffDiagU_SM;
  U.bottomRightCorner(numRows_BottOffDiag,HODLR_Root->bottOffDiagRank) = HODLR_Root->bottOffDiagU_SM;
  VT.topRightCorner(HODLR_Root->topOffDiagRank,numCols_TopOffDiag) = HODLR_Root->topOffDiagV.transpose();
  VT.bottomLeftCorner(HODLR_Root->bottOffDiagRank,numCols_BottOffDiag) = HODLR_Root->bottOffDiagV.transpose();
  
  Eigen::MatrixXd detMatrix = Eigen::MatrixXd::Identity(rankTotal,rankTotal) + VT * U; 
  Eigen::PartialPivLU<Eigen::MatrixXd> luFact(detMatrix);
  Eigen::MatrixXd luMatrix = luFact.matrixLU();
  
  for (int i = 0; i < rankTotal; i++){
    determinant_ *= luMatrix(i,i);
    logAbsDeterminant_ += log(fabs(luMatrix(i,i)));
  }
 
  calcDeterminant(HODLR_Root->left);
  calcDeterminant(HODLR_Root->right);
  
}

