#include "HODLR_Matrix.hpp"

const double pi = 3.14159265359;

void HODLR_Matrix::setDefaultValues(){
  
  LR_Tolerance = 1e-6;
  minValueACA  = 0;
  
  recLU_FactorizationTime      = 0;
  recLU_SolveTime              = 0;
  LR_ComputationTime           = 0;
  extendedSp_Size              = 0;
  extendedSp_AssemblyTime      = 0;
  extendedSp_FactorizationTime = 0;
  extendedSp_SolveTime         = 0;
  iter_SolveTime               = 0;
  matrixSize                   = 0; 
  matrixNumRows                = 0;
  matrixNumCols                = 0; 


  LRStoredInTree         = false;
  createdRecLUfactorTree = false;
  assembled_ExtendedSp   = false;
  saveExtendedSp_Matrix  = false;
  freeMatrixMemory       = false;
  freeMatrixMemory_Sp    = false;
  matrixDataAvail        = false;
  matrixDataAvail_Sp     = false;
  isSquareMatrix         = false;

  printLevelRankInfo     = false;
  printLevelAccuracy     = false;
  printLevelInfo         = false;
  printResultInfo        = false;  
 
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
  matrixSize    = inputMatrix.rows();
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
}

HODLR_Matrix::HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix){
  setDefaultValues();
  matrixData_Sp = inputMatrix;
  sizeThreshold = 30;
  matrixDataAvail_Sp = true;
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  indexTree.set_def_LRMethod("PS_Sparse");
  matrixSize    = inputMatrix.rows();
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
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
  matrixDataAvail = true;  
}

HODLR_Matrix::HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,int inputSizeThreshold){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true); // Currently unable to build trees for non squared matrices
  matrixData_Sp = inputMatrix;
  matrixSize    = inputMatrix.rows();
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.createDefaultTree(matrixSize);
  indexTree.set_def_LRMethod("PS_Sparse");
  matrixDataAvail_Sp = true;
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
  matrixDataAvail = true;
  
}

HODLR_Matrix::HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix, int inputSizeThreshold, user_IndexTree &input_IndexTree){
  setDefaultValues();
  isSquareMatrix = (inputMatrix.rows() == inputMatrix.cols());
  assert(isSquareMatrix == true);  // Currently unable to build trees for non squared matrices
  matrixData_Sp = inputMatrix;
  matrixNumRows = inputMatrix.rows();
  matrixNumCols = inputMatrix.cols();
  matrixSize    = inputMatrix.rows();
  sizeThreshold = inputSizeThreshold;
  indexTree.set_sizeThreshold(sizeThreshold);
  indexTree.createFromUsrTree(matrixSize,input_IndexTree);
  indexTree.set_def_LRMethod("PS_Sparse");
  matrixDataAvail_Sp = true;
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

  recLU_FactorizationTime      = rhs.recLU_FactorizationTime;
  recLU_SolveTime              = rhs.recLU_SolveTime;
  LR_ComputationTime           = rhs.LR_ComputationTime;
  extendedSp_AssemblyTime      = rhs.extendedSp_AssemblyTime;
  extendedSp_FactorizationTime = rhs.extendedSp_FactorizationTime;
  extendedSp_SolveTime         = rhs.extendedSp_SolveTime;
  iter_SolveTime               = rhs.iter_SolveTime;

  LRStoredInTree         = rhs.LRStoredInTree;
  createdRecLUfactorTree = rhs.createdRecLUfactorTree;
  assembled_ExtendedSp   = rhs.assembled_ExtendedSp;
  saveExtendedSp_Matrix  = rhs.saveExtendedSp_Matrix;
  freeMatrixMemory       = rhs.freeMatrixMemory;
  matrixDataAvail        = rhs.matrixDataAvail;    
  matrixDataAvail_Sp     = rhs.matrixDataAvail_Sp;
  isSquareMatrix         = rhs.isSquareMatrix;

  LR_Tolerance           = rhs.LR_Tolerance;
  minValueACA            = rhs.minValueACA; 

 
  matrixData          = rhs.matrixData;
  matrixData_Sp       = rhs.matrixData_Sp;
  extendedSp_Solver   = rhs.extendedSp_Solver; 
  extendedSp_SavePath = rhs.extendedSp_SavePath;
  indexTree           = rhs.indexTree;
  //recLUfactorTree needs to be copied :)) TODO

}

void HODLR_Matrix::storeLRinTree(){
  assert(indexTree.rootNode != NULL);
  assert((matrixDataAvail == true) || (matrixDataAvail_Sp == true));
  storeLRinTree(indexTree.rootNode);
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
    int numRows = HODLR_Root->max_i - HODLR_Root->min_i + 1;
    int numCols = HODLR_Root->max_j - HODLR_Root->min_j + 1;
    if (matrixDataAvail_Sp == true)
      HODLR_Root->leafMatrix = Eigen::MatrixXd(matrixData_Sp.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols));
    else
      HODLR_Root->leafMatrix = matrixData.block(HODLR_Root->min_i,HODLR_Root->min_j,numRows,numCols);
    return;
  }
  // Calculate the LR factorizations
  if (HODLR_Root->LR_Method == "partialPiv_ACA"){
    partialPivACA_LowRankApprox(HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_i,HODLR_Root->splitIndex_j + 1,HODLR_Root->max_j,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank);
    partialPivACA_LowRankApprox(HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->max_i,HODLR_Root->min_j,HODLR_Root->splitIndex_j,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank);
    HODLR_Root->topOffDiagK = Eigen::MatrixXd::Identity(HODLR_Root->topOffDiagRank, HODLR_Root->topOffDiagRank);
    HODLR_Root->bottOffDiagK = Eigen::MatrixXd::Identity(HODLR_Root->bottOffDiagRank, HODLR_Root->bottOffDiagRank);
    
  }else if (HODLR_Root->LR_Method == "fullPiv_ACA"){
    fullPivACA_LowRankApprox(HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->min_i,HODLR_Root->splitIndex_i,HODLR_Root->splitIndex_j + 1,HODLR_Root->max_j,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank);
    fullPivACA_LowRankApprox(HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->splitIndex_i + 1,HODLR_Root->max_i,HODLR_Root->min_j,HODLR_Root->splitIndex_j,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank);
    HODLR_Root->topOffDiagK = Eigen::MatrixXd::Identity(HODLR_Root->topOffDiagRank, HODLR_Root->topOffDiagRank);
    HODLR_Root->bottOffDiagK = Eigen::MatrixXd::Identity(HODLR_Root->bottOffDiagRank, HODLR_Root->bottOffDiagRank);

  }else if (HODLR_Root->LR_Method == "PS_Cheby"){
    PS_LowRankApprox(HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->topOffDiagK,HODLR_Root->min_i,HODLR_Root->splitIndex_i,HODLR_Root->splitIndex_j + 1,HODLR_Root->max_j,LR_Tolerance,HODLR_Root->topOffDiagRank,"Chebyshev",HODLR_Root->topOffDiag_minRank);
    PS_LowRankApprox(HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->bottOffDiagK,HODLR_Root->splitIndex_i + 1,HODLR_Root->max_i,HODLR_Root->min_j,HODLR_Root->splitIndex_j,LR_Tolerance,HODLR_Root->bottOffDiagRank,"Chebyshev",HODLR_Root->bottOffDiag_minRank);

  }else if (HODLR_Root->LR_Method == "SVD"){ 
    SVD_LowRankApprox(HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->topOffDiagK, HODLR_Root->min_i,HODLR_Root->splitIndex_i,HODLR_Root->splitIndex_j + 1,HODLR_Root->max_j,LR_Tolerance,HODLR_Root->topOffDiagRank,HODLR_Root->topOffDiag_minRank);
    SVD_LowRankApprox(HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->bottOffDiagK, HODLR_Root->splitIndex_i + 1,HODLR_Root->max_i,HODLR_Root->min_j,HODLR_Root->splitIndex_j,LR_Tolerance,HODLR_Root->bottOffDiagRank,HODLR_Root->bottOffDiag_minRank);

  }else if(HODLR_Root->LR_Method == "PS_Sparse"){
    PS_LowRankApprox_Sp(HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagV,HODLR_Root->topOffDiagK,HODLR_Root->min_i,HODLR_Root->splitIndex_i,HODLR_Root->splitIndex_j + 1,HODLR_Root->max_j,LR_Tolerance,HODLR_Root->topOffDiagRank);
    PS_LowRankApprox_Sp(HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagV,HODLR_Root->bottOffDiagK,HODLR_Root->splitIndex_i + 1,HODLR_Root->max_i,HODLR_Root->min_j,HODLR_Root->splitIndex_j,LR_Tolerance,HODLR_Root->bottOffDiagRank);

  }else{
    std::cout<<"Error!. Invalid low-rank approximation scheme."<<std::endl;
    exit(EXIT_FAILURE);
  }
  
  storeLRinTree(HODLR_Root->left);
  storeLRinTree(HODLR_Root->right);
  
}


void HODLR_Matrix::recLU_Factorize(){
  recLUfactorTree.rootNode = new recLU_FactorTree::node;
  (recLUfactorTree.rootNode)->isLeaf = false; 
  Eigen::VectorXd dummyF = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2); 
  Eigen::MatrixXd dummyX = recLU_Factorize(dummyF, indexTree.rootNode, recLUfactorTree.rootNode);
}

Eigen::MatrixXd HODLR_Matrix::recLU_Factorize(const Eigen::MatrixXd & input_RHS,const HODLR_Tree::node* HODLR_Root, recLU_FactorTree::node* factorRoot){
  // Base case
  if (HODLR_Root->isLeaf == true){
    factorRoot->isLeaf = true;
    factorRoot->left = NULL;
    factorRoot->right = NULL;
    //int leafMatrixSize = HODLR_Root->max_i - HODLR_Root->min_i + 1;
    //Eigen::MatrixXd leafMatrix = matrixData.block(HODLR_Root->min_i,HODLR_Root->min_j,leafMatrixSize,leafMatrixSize);
    LUDecompose(HODLR_Root->leafMatrix,factorRoot->LU_leaf,factorRoot->P_leaf);
    Eigen::MatrixXd y = (factorRoot->LU_leaf).triangularView<Eigen::UnitLower>().solve(factorRoot->P_leaf * input_RHS);
    return (factorRoot->LU_leaf).triangularView<Eigen::Upper>().solve(y);
  }
  
  // Low Rank Approximation
  Eigen::MatrixXd WB,WC;
  Eigen::MatrixXd VB,VC;
  int calculatedRankB,calculatedRankC;
  int topDiagSize = HODLR_Root->splitIndex_i - HODLR_Root->min_i + 1;
  int bottDiagSize = HODLR_Root->max_i - HODLR_Root->splitIndex_i;
  int parentRHS_Cols = input_RHS.cols();
  WB = HODLR_Root->topOffDiagU * HODLR_Root->topOffDiagK;
  VB = HODLR_Root->topOffDiagV;
  WC = HODLR_Root->bottOffDiagU * HODLR_Root->bottOffDiagK;
  VC = HODLR_Root->bottOffDiagV;
  calculatedRankB = HODLR_Root->topOffDiagRank;
  calculatedRankC = HODLR_Root->bottOffDiagRank;
  
  // Factorize and Solve for top diagonal matrix
  int topOffDiag_LR_Cols = WB.cols();
  int topDiagRHS_Cols = parentRHS_Cols + topOffDiag_LR_Cols;
  Eigen::MatrixXd topDiagRHS = Eigen::MatrixXd::Zero(topDiagSize,topDiagRHS_Cols);
  topDiagRHS.leftCols(parentRHS_Cols) = input_RHS.topRows(topDiagSize);
  topDiagRHS.rightCols(topOffDiag_LR_Cols) = WB;  
  recLU_FactorTree::node* leftFactorRoot = new recLU_FactorTree::node;
  leftFactorRoot->isLeaf = false;
  factorRoot->left = leftFactorRoot;
  Eigen::MatrixXd topDiagSoln = recLU_Factorize(topDiagRHS,HODLR_Root->left,leftFactorRoot);
  Eigen::MatrixXd topDiagSoln_pRHS = topDiagSoln.leftCols(parentRHS_Cols);
  Eigen::MatrixXd topDiagSoln_LR = topDiagSoln.rightCols(topOffDiag_LR_Cols);
  factorRoot->topDiagSoln_LR = topDiagSoln_LR;
  
  // Factorize and Solve for bottom diagonal matrix
  int bottOffDiag_LR_Cols = WC.cols();
  int bottDiagRHS_Cols = parentRHS_Cols + bottOffDiag_LR_Cols;
  Eigen::MatrixXd bottDiagRHS = Eigen::MatrixXd::Zero(bottDiagSize,bottDiagRHS_Cols);
  bottDiagRHS.leftCols(parentRHS_Cols) = input_RHS.bottomRows(bottDiagSize);
  bottDiagRHS.rightCols(bottOffDiag_LR_Cols) = WC;
  recLU_FactorTree::node* rightFactorRoot = new recLU_FactorTree::node;
  rightFactorRoot->isLeaf = false;
  factorRoot->right = rightFactorRoot;
  Eigen::MatrixXd bottDiagSoln = recLU_Factorize(bottDiagRHS,HODLR_Root->right,rightFactorRoot);
  Eigen::MatrixXd bottDiagSoln_pRHS = bottDiagSoln.leftCols(parentRHS_Cols);
  Eigen::MatrixXd bottDiagSoln_LR = bottDiagSoln.rightCols(bottOffDiag_LR_Cols);
  factorRoot->bottDiagSoln_LR = bottDiagSoln_LR;

  // Update the remaining of the matrix(generate Schur complement (S matrix));
  
  int Sdim = calculatedRankB + calculatedRankC;  
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(Sdim,Sdim);
  Eigen::MatrixXd IC = Eigen::MatrixXd::Identity(calculatedRankC,calculatedRankC);
  Eigen::MatrixXd IB = Eigen::MatrixXd::Identity(calculatedRankB,calculatedRankB);
  
  S.topLeftCorner(calculatedRankC,calculatedRankC) = IC;
  S.bottomRightCorner(calculatedRankB,calculatedRankB) = IB;
  
  S.topRightCorner(calculatedRankC,calculatedRankB) = (VC.transpose() * topDiagSoln_LR);
  S.bottomLeftCorner(calculatedRankB,calculatedRankC) = (VB.transpose() * bottDiagSoln_LR);

  LUDecompose(S,factorRoot->LU_S,factorRoot->P_S);	

  Eigen::MatrixXd L_S = factorRoot->LU_S.triangularView<Eigen::UnitLower>();
  Eigen::MatrixXd U_S = factorRoot->LU_S.triangularView<Eigen::Upper>();
  
  Eigen::MatrixXd l1 = L_S.topLeftCorner(calculatedRankC,calculatedRankC);
  Eigen::MatrixXd l2 = L_S.bottomLeftCorner(calculatedRankB,calculatedRankC);
  Eigen::MatrixXd l3 = L_S.bottomRightCorner(calculatedRankB,calculatedRankB);
  
  Eigen::MatrixXd u1 = U_S.topLeftCorner(calculatedRankC,calculatedRankC);
  Eigen::MatrixXd u2 = U_S.topRightCorner(calculatedRankC,calculatedRankB);
  Eigen::MatrixXd u3 = U_S.bottomRightCorner(calculatedRankB,calculatedRankB);
  
  Eigen::MatrixXd equation1 = VC.transpose() * topDiagSoln_pRHS;
  Eigen::MatrixXd equation2 = VB.transpose() * bottDiagSoln_pRHS;

  Eigen::MatrixXd setOfEqns(Sdim,equation1.cols());
  setOfEqns.topRows(calculatedRankC) = equation1;
  setOfEqns.bottomRows(calculatedRankB) = equation2;
  
  setOfEqns = factorRoot->P_S*setOfEqns;
  equation1 = setOfEqns.topRows(calculatedRankC);
  equation2 = setOfEqns.bottomRows(calculatedRankB);
  
  Eigen::MatrixXd zeta1 = l1.triangularView<Eigen::UnitLower>().solve(equation1);
  Eigen::MatrixXd zeta2 = l3.triangularView<Eigen::UnitLower>().solve(equation2- l2 * zeta1);	
  
  // Calculate y1 & y2
  Eigen::MatrixXd y2 = u3.triangularView<Eigen::Upper>().solve(zeta2);
  Eigen::MatrixXd y1 = u1.triangularView<Eigen::Upper>().solve(zeta1 - u2 * y2);

  // Obtain x
  Eigen::MatrixXd x2 = bottDiagSoln_pRHS - bottDiagSoln_LR * y1;
  Eigen::MatrixXd x1 = topDiagSoln_pRHS - topDiagSoln_LR * y2;
  
  Eigen::MatrixXd result(x1.rows() + x2.rows(),x1.cols());
  result.topRows(x1.rows()) = x1;
  result.bottomRows(x2.rows()) = x2;
  
  return result;

}

Eigen::MatrixXd HODLR_Matrix::recLU_Solve(const Eigen::MatrixXd & input_RHS,const HODLR_Tree::node* HODLR_Root,const recLU_FactorTree::node* factorRoot){

  // Base case
  if (HODLR_Root->isLeaf == true){
    Eigen::MatrixXd y = (factorRoot->LU_leaf).triangularView<Eigen::UnitLower>().solve(factorRoot->P_leaf * input_RHS);
    return (factorRoot->LU_leaf).triangularView<Eigen::Upper>().solve(y);
 
  }
  
  // Low Rank Approximation
  Eigen::MatrixXd WB,WC;
  Eigen::MatrixXd VB,VC;
  int calculatedRankB,calculatedRankC;
  WB = HODLR_Root->topOffDiagU * HODLR_Root->topOffDiagK;
  VB = HODLR_Root->topOffDiagV;
  WC = HODLR_Root->bottOffDiagU * HODLR_Root->bottOffDiagK;
  VC = HODLR_Root->bottOffDiagV;
    
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
    if (matrixDataAvail){
      Eigen::MatrixXd topDiag = matrixData.block(HODLR_Root->min_i,HODLR_Root->min_j,topDiagSize,topDiagSize);
      Eigen::MatrixXd topOffDiag = matrixData.block(HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1 ,topDiagSize,bottDiagSize);
      Eigen::MatrixXd bottOffDiag = matrixData.block(HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,bottDiagSize,topDiagSize);
      Eigen::MatrixXd bottDiag = matrixData.block(HODLR_Root->splitIndex_i + 1,HODLR_Root->splitIndex_j + 1,bottDiagSize,bottDiagSize);
      std::cout<<"**********************************************************************"<<std::endl;
      std::cout<<"**********************************************************************"<<std::endl;
      std::cout<<"Current Recursion Level                     = "<<HODLR_Root->currLevel<<std::endl;
      std::cout<<"Top Diagonal Matrix min_i                   = "<<HODLR_Root->min_i<<" | "<<"Top Diagonal Matrix max_i = "<<HODLR_Root->splitIndex_i<<std::endl;
      std::cout<<"Top Diagonal Matrix min_j                   = "<<HODLR_Root->min_j<<" | "<<"Top Diagonal Matrix max_i = "<<HODLR_Root->splitIndex_j<<std::endl;
      std::cout<<"Top Diagonal Matrix Solve Rel Error         = "<<(topDiag * topDiagSoln-topDiagRHS).norm()/(topDiagRHS.norm())<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"Top Off-Diagonal Approximation Rel Error    = "<<(topOffDiag - WB*VB.transpose()).norm()/topOffDiag.norm()<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"Bottom Diagonal Matrix min_i                = "<<HODLR_Root->splitIndex_i + 1<<" | "<<"Bottom Diagonal Matrix max_i = "<<HODLR_Root->max_i<<std::endl;
      std::cout<<"Bottom Diagonal Matrix min_j                = "<<HODLR_Root->splitIndex_j + 1<<" | "<<"Bottom Diagonal Matrix max_i = "<<HODLR_Root->max_j<<std::endl;
      std::cout<<"Bottom Diagonal Matrix Solve Rel Error      = "<<(bottDiag*bottDiagSoln-bottDiagRHS).norm()/(bottDiagRHS.norm())<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"Bottom Off-Diagonal Approximation Rel Error = "<<(bottOffDiag-WC*VC.transpose()).norm()/bottOffDiag.norm()<<std::endl;
    }else{
      std::cout<<"Error! Matrix data has been deleted from memory."<<std::endl;
    }
  }

  if (printLevelRankInfo){
    if (matrixDataAvail){
      Eigen::MatrixXd topOffDiag = matrixData.block(HODLR_Root->min_i,HODLR_Root->splitIndex_j + 1 ,topDiagSize,bottDiagSize);
      Eigen::MatrixXd bottOffDiag = matrixData.block(HODLR_Root->splitIndex_i + 1,HODLR_Root->min_j,bottDiagSize,topDiagSize);
      int actualRankTop = SVD_LowRankApprox(topOffDiag, LR_Tolerance);
      int actualRankBott =  SVD_LowRankApprox(bottOffDiag, LR_Tolerance);
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
  int Sdim = calculatedRankB + calculatedRankC;  
  
  Eigen::MatrixXd L_S = factorRoot->LU_S.triangularView<Eigen::UnitLower>();
  Eigen::MatrixXd U_S = factorRoot->LU_S.triangularView<Eigen::Upper>();

  Eigen::MatrixXd l1 = L_S.topLeftCorner(calculatedRankC,calculatedRankC);
  Eigen::MatrixXd l2 = L_S.bottomLeftCorner(calculatedRankB,calculatedRankC);
  Eigen::MatrixXd l3 = L_S.bottomRightCorner(calculatedRankB,calculatedRankB);
  
  Eigen::MatrixXd u1 = U_S.topLeftCorner(calculatedRankC,calculatedRankC);
  Eigen::MatrixXd u2 = U_S.topRightCorner(calculatedRankC,calculatedRankB);
  Eigen::MatrixXd u3 = U_S.bottomRightCorner(calculatedRankB,calculatedRankB);
  
  Eigen::MatrixXd equation1 = VC.transpose() * topDiagSoln_pRHS;
  Eigen::MatrixXd equation2 = VB.transpose() * bottDiagSoln_pRHS;

  Eigen::MatrixXd setOfEqns(Sdim,equation1.cols());
  setOfEqns.topRows(calculatedRankC) = equation1;
  setOfEqns.bottomRows(calculatedRankB) = equation2;
  
  setOfEqns = factorRoot->P_S * setOfEqns;
  equation1 = setOfEqns.topRows(calculatedRankC);
  equation2 = setOfEqns.bottomRows(calculatedRankB);
  
  Eigen::MatrixXd zeta1 = l1.triangularView<Eigen::UnitLower>().solve(equation1);
  Eigen::MatrixXd zeta2 = l3.triangularView<Eigen::UnitLower>().solve(equation2 -l2 * zeta1);	
  
  // Calculate y1 & y2
  Eigen::MatrixXd y2 = u3.triangularView<Eigen::Upper>().solve(zeta2);
  Eigen::MatrixXd y1 = u1.triangularView<Eigen::Upper>().solve(zeta1-u2*y2);
	
  // Obtain x
  Eigen::MatrixXd x2 = bottDiagSoln_pRHS - factorRoot->bottDiagSoln_LR * y1;
  Eigen::MatrixXd x1 = topDiagSoln_pRHS - factorRoot->topDiagSoln_LR * y2;
  
  Eigen::MatrixXd result(x1.rows() + x2.rows(),x1.cols());
  result.topRows(x1.rows()) = x1;
  result.bottomRows(x2.rows()) = x2;
  
  return result;
}

Eigen::MatrixXd HODLR_Matrix::recLU_Solve(const Eigen::MatrixXd & input_RHS){
  
  assert(isSquareMatrix == true);
  assert(input_RHS.rows() == matrixSize);
  if (indexTree.rootNode == NULL){
    indexTree.set_sizeThreshold(sizeThreshold);
    indexTree.createDefaultTree(matrixSize);
  }
  
  if (LRStoredInTree == false){
    double startTime = clock();
    storeLRinTree();
    double endTime = clock();
    LR_ComputationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    LRStoredInTree = true;
  }

  if (createdRecLUfactorTree == false){
    double startTime = clock();
    recLU_Factorize();
    double endTime = clock();
    recLU_FactorizationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    createdRecLUfactorTree = true;
  }
  
  Eigen::MatrixXd solution;
  double startTime = clock();
  solution = recLU_Solve(input_RHS,indexTree.rootNode,recLUfactorTree.rootNode);
  double endTime = clock();
  recLU_SolveTime = (endTime-startTime)/CLOCKS_PER_SEC;
  if (printResultInfo){
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"Solver Type                      = recLU"<<std::endl;
    std::cout<<"Low-Rank Computation Time        = "<<LR_ComputationTime<<" seconds"<<std::endl;
    std::cout<<"Factorization Time               = "<<recLU_FactorizationTime<<" seconds"<<std::endl;
    std::cout<<"Solve Time                       = "<<recLU_SolveTime<<" seconds"<<std::endl; 
    std::cout<<"LR Tolerance                     = "<<LR_Tolerance<<std::endl;
    std::cout<<"Residual l2 Relative Error       = "<<((matrixData * solution) - input_RHS).norm()/input_RHS.norm()<<std::endl;
  }
  

  return solution;

}

void HODLR_Matrix::recLU_Compute(){
  if (indexTree.rootNode == NULL){
    indexTree.set_sizeThreshold(sizeThreshold);
    indexTree.createDefaultTree(matrixSize);
  }
  if (LRStoredInTree == false){
    double startTime = clock();
    storeLRinTree();
    double endTime = clock();
    LR_ComputationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    LRStoredInTree = true;
  }
  if (createdRecLUfactorTree == false){
    double startTime = clock();
    recLU_Factorize();
    double endTime = clock();
    recLU_FactorizationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    createdRecLUfactorTree = true;
  }
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
    //findNodesAtLevel(indexTree.rootNode,i,currLevelNodesVec);
    currLevelNodesVec = indexTree.nodeLevelVec[i];
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
          
      Eigen::MatrixXd topOffDiagK  = currLevelNodesVec[j]->topOffDiagK;
      Eigen::MatrixXd bottOffDiagK = currLevelNodesVec[j]->bottOffDiagK;
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
  //findLeafNodes(indexTree.rootNode,leafNodesVec);
  leafNodesVec = indexTree.leafNodesVec;
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
  }
  
  if (LRStoredInTree == false){
    double startTime = clock();
    storeLRinTree();
    double endTime = clock();
    LR_ComputationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    LRStoredInTree = true;
  }

    
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
  
  if (printResultInfo){
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"Solver Type                      = extendedSp"<<std::endl;
    std::cout<<"Low-Rank Computation Time        = "<<LR_ComputationTime<<" seconds"<<std::endl;
    std::cout<<"Assembly Time                    = "<<extendedSp_AssemblyTime<<" seconds"<<std::endl;
    std::cout<<"Factorization Time               = "<<extendedSp_FactorizationTime<<" seconds"<<std::
endl;
    std::cout<<"Solve Time                       = "<<extendedSp_SolveTime<<" seconds"<<std::endl; 
    std::cout<<"LR Tolerance                     = "<<LR_Tolerance<<std::endl;
    std::cout<<"Residual l2 Relative Error       = "<<((matrixData * solution)-input_RHS).norm()/input_RHS.norm()<<std::endl;
  }
  
  return solution;

}

Eigen::MatrixXd HODLR_Matrix::oneStep_Iterate(const Eigen::MatrixXd & prevStep_result,const Eigen::MatrixXd & RHS, const Eigen::MatrixXd & initSolveGuess, Eigen::MatrixXd & prevStep_Product,const std::string directSolve_Method){

  prevStep_Product = matrixData * prevStep_result;
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
  set_def_LRMethod(input_LRMethod);

  bool save_printResultInfo = printResultInfo;
  printResultInfo = false;

  if (indexTree.rootNode == NULL){
    indexTree.set_sizeThreshold(sizeThreshold);
    indexTree.createDefaultTree(matrixSize);
  }
  
  if (LRStoredInTree == false){
    double startTime = clock();
    storeLRinTree();
    double endTime = clock();
    LR_ComputationTime = (endTime-startTime)/CLOCKS_PER_SEC;
    LRStoredInTree = true;
  }

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
  while (tolerance > stop_tolerance){
    nextStep_Soln = oneStep_Iterate(currStep_Soln,input_RHS,init_Guess,currStep_Product,directSolve_Method);
    //tolerance = (nextStep_Soln - currStep_Soln).norm()/currStep_Soln.norm();
    tolerance = (currStep_Product - input_RHS).norm()/input_RHS.norm();
    currStep_Soln = nextStep_Soln;      
    //std::cout<<num_Iter<<" "<<tolerance<<std::endl;
    num_Iter ++;
    if (num_Iter > maxIterations)
      break;
  }
  Eigen::MatrixXd solution = currStep_Soln;
  double endTime = clock();
  iter_SolveTime = (endTime-startTime)/CLOCKS_PER_SEC;
  printResultInfo = save_printResultInfo;
  if (printResultInfo){
    std::cout<<"**************************************************"<<std::endl;
    std::cout<<"Solver Type                      = iterative"<<std::endl;
    std::cout<<"Low-Rank Computation Time        = "<<LR_ComputationTime<<" seconds"<<std::endl;
    std::cout<<"Solve Time                       = "<<iter_SolveTime<<" seconds"<<std::endl; 
    std::cout<<"LR Tolerance                     = "<<LR_Tolerance<<std::endl;
    std::cout<<"Number of Iterations             = "<<num_Iter<<std::endl;
    std::cout<<"Residual l2 Relative Error       = "<<((matrixData * solution) - input_RHS).norm()/input_RHS.norm()<<std::endl;
  }
  
  // restore previous state;
  //set_LRTolerance(prev_LRTolerance);
  //set_def_LRMethod(prev_LRMethod);
  return solution;
}


int HODLR_Matrix::SVD_LowRankApprox(const Eigen::MatrixXd & inputMatrix, const double accuracy, Eigen::MatrixXd* Wptr, Eigen::MatrixXd* Vptr, Eigen::MatrixXd* Kptr, int minRank) const{
  
  // Use svd to calculate the low-rank approximation
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(inputMatrix,Eigen::ComputeThinU|Eigen::ComputeThinV);
  
  // Calculate rank
  Eigen::VectorXd singularValues = svd.singularValues();
  int nSingularValues = singularValues.rows();
  int rank = nSingularValues;
  for (int i = 0; i < nSingularValues; i++)
    if (singularValues(i) < accuracy){
      rank = i + 1;
      break;
    } 
  if (rank < minRank)
    rank = minRank;

  if (Wptr != NULL){
    *Wptr = svd.matrixU().leftCols(rank);
  }
  if (Vptr != NULL){
    *Vptr = svd.matrixV().leftCols(rank);
  }
  if (Kptr != NULL){
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(rank,rank);
    for (int i = 0; i < rank; i++)
      K(i,i) = singularValues(i);
    *Kptr = K;
  }
  return rank;
}


void HODLR_Matrix::SVD_LowRankApprox(Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int max_i, const int min_j, const int max_j, const double tolerance, int & calculatedRank, const int minRank ) const{

  int nRows = max_i-min_i+1;
  int nCols = max_j-min_j+1;
  Eigen::MatrixXd lowRankMatrix = matrixData.block(min_i,min_j,nRows,nCols);
  calculatedRank = SVD_LowRankApprox(lowRankMatrix, tolerance, &W, &V, &K, minRank);
   
}


void HODLR_Matrix::PS_LowRankApprox_Sp(Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int max_i,const int min_j, const int max_j, const double tolerance, int &calculatedRank)const{
	
  int nRows = max_i-min_i+1;
  int nCols = max_j-min_j+1;
  //int maxRank = std::min(nRows,nCols);
  //int numPoints;
  std::set<int> rowIdxSet,colIdxSet;
  /*
  if (minRank > 0)
    numPoints = minRank;
  else
    numPoints = maxRank/25 + 1;
  

  double absFrobNormDiff = 1;
  double relFrobNormDiff = 1;
  double approxError = 1;
  */
  
  Eigen::SparseMatrix<double> lowRankMatrix_Sp = matrixData_Sp.block(min_i,min_j,nRows,nCols);
  //find numPoints
  if (lowRankMatrix_Sp.nonZeros() == 0){
    calculatedRank = 1;
    W = Eigen::MatrixXd::Zero(nRows,1);
    K = Eigen::MatrixXd::Zero(1,1);
    V = Eigen::MatrixXd::Zero(nCols,1);
    return;
  }
  for (int k = 0; k < lowRankMatrix_Sp.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(lowRankMatrix_Sp,k); it; ++it){
      rowIdxSet.insert(it.row());
      colIdxSet.insert(it.col());
    }
  std::vector<int> rowIndex(rowIdxSet.begin(),rowIdxSet.end());
  std::vector<int> colIndex(colIdxSet.begin(),colIdxSet.end());
  
  //while ((approxError > tolerance) && (numPoints <= maxRank)){	
    
  //create row and col index using Chebyshev or uniform nodes
  /*
  Eigen::VectorXi rowIndex(numPoints),colIndex(numPoints);
  if (pointChoosingMethod == "Chebyshev")
    for (int i = 0 ; i < numPoints; i++){
      rowIndex(i) = floor((nRows + nRows * cos(pi*(2*i+1)/(2*numPoints)))/2);
      colIndex(i) = floor((nCols + nCols * cos(pi*(2*i+1)/(2*numPoints)))/2);
    }
  if (pointChoosingMethod == "Uniform" && numPoints > 1){
    int rowStride = nRows/(numPoints-1);
    int colStride = nCols/(numPoints-1);
    if (rowStride == 0)
      rowStride = 1;
    if (colStride == 0)
      colStride = 1;
    for (int i=0 ; i<numPoints; i++){
      rowIndex(i) = i * rowStride;
      colIndex(i) = i * colStride;
    }
  }
  */
  
  //choose rows and columns and do the low-rank approximation
  Eigen::MatrixXd dummyW,dummyK,dummyV;
  extractRowsCols(W,K,V,Eigen::MatrixXd(lowRankMatrix_Sp),rowIndex,colIndex);
  calculatedRank = W.cols();
  
  //obtain stopping criterion
  /*
  Eigen::VectorXi rowIndexTest = (rowIndex.head(numPoints-1)+rowIndex.tail(numPoints-1))/2;
  Eigen::VectorXi colIndexTest = (colIndex.head(numPoints-1)+colIndex.tail(numPoints-1))/2;
  
  Eigen::MatrixXd sampleColsTest(nRows,numPoints-1);
  Eigen::MatrixXd approxColsTest(nRows,numPoints-1);
  Eigen::MatrixXd sampleRowsTest(numPoints-1,nCols);
  Eigen::MatrixXd approxRowsTest(numPoints-1,nCols);
  
  //fill KTempApprox
  for (int i = 0; i < numPoints-1;i++){
    sampleRowsTest.row(i) = lowRankMatrix.row(rowIndexTest(i));
    approxRowsTest.row(i) = (W * K).row(rowIndexTest(i))*V.transpose();
    sampleColsTest.col(i) = lowRankMatrix.col(colIndexTest(i));
    approxColsTest.col(i) = (W * K)*(V.row(colIndexTest(i)).transpose());	
  }
  
  Eigen::MatrixXd sampleColsTestBlock = sampleColsTest.block(numPoints-1,0,nRows-numPoints+1,numPoints-1);
  Eigen::MatrixXd approxColsTestBlock = approxColsTest.block(numPoints-1,0,nRows-numPoints+1,numPoints-1);
  absFrobNormDiff = (sampleRowsTest-approxRowsTest).norm()+(sampleColsTestBlock-approxColsTestBlock).norm();
  relFrobNormDiff = absFrobNormDiff/(sampleRowsTest.norm()+sampleColsTestBlock.norm());
  approxError = relFrobNormDiff*(sqrt((nRows*nCols)/((numPoints-1)*(nCols+nRows-numPoints+1))));
  numPoints *= 1.5;
  //}
  calculatedRank = W.cols(); 
  */
}

void HODLR_Matrix::PS_LowRankApprox(Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int max_i,const int min_j, const int max_j, const double tolerance, int &calculatedRank, const std::string pointChoosingMethod,const int minRank)const{
	
  int nRows = max_i-min_i+1;
  int nCols = max_j-min_j+1;
  int maxRank = std::min(nRows,nCols);
  int numPoints;

  if (minRank > 0)
    numPoints = minRank;
  else
    numPoints = maxRank/25 + 1;

  double absFrobNormDiff = 1;
  double relFrobNormDiff = 1;
  double approxError = 1;
  Eigen::MatrixXd lowRankMatrix = matrixData.block(min_i,min_j,nRows,nCols);
  
  while ((approxError > tolerance) && (numPoints <= maxRank)){	
    
    //create row and col index using Chebyshev or uniform nodes
    std::vector<int> rowIndex_Vec(numPoints),colIndex_Vec(numPoints);
    Eigen::VectorXi  rowIndex(numPoints),colIndex(numPoints);
    if (pointChoosingMethod == "Chebyshev")
      for (int i = 0 ; i < numPoints; i++){
	rowIndex(i) = floor((nRows + nRows * cos(pi * (2 * i + 1)/(2 * numPoints))) / 2);
	colIndex(i) = floor((nCols + nCols * cos(pi * (2 * i + 1)/(2 * numPoints))) / 2);
	rowIndex_Vec[i] = rowIndex(i);
	colIndex_Vec[i] = colIndex(i);
      }
    if (pointChoosingMethod == "Uniform" && numPoints > 1){
      int rowStride = nRows / (numPoints - 1);
      int colStride = nCols / (numPoints - 1);
      if (rowStride == 0)
	rowStride = 1;
      if (colStride == 0)
	colStride = 1;
      for (int i=0 ; i<numPoints; i++){
	rowIndex(i) = i * rowStride;
	colIndex(i) = i * colStride;
	rowIndex_Vec[i] = rowIndex(i);
	colIndex_Vec[i] = colIndex(i);
      }
    }

    //choose rows and columns and do the low-rank approximation
    extractRowsCols(W,K,V,lowRankMatrix,rowIndex_Vec,colIndex_Vec);
    calculatedRank = W.cols();

    //obtain stopping criterion
    Eigen::VectorXi rowIndexTest = (rowIndex.head(numPoints-1)+rowIndex.tail(numPoints-1))/2;
    Eigen::VectorXi colIndexTest = (colIndex.head(numPoints-1)+colIndex.tail(numPoints-1))/2;
    
    Eigen::MatrixXd sampleColsTest(nRows,numPoints-1);
    Eigen::MatrixXd approxColsTest(nRows,numPoints-1);
    Eigen::MatrixXd sampleRowsTest(numPoints-1,nCols);
    Eigen::MatrixXd approxRowsTest(numPoints-1,nCols);
    
    //fill KTempApprox
    for (int i = 0; i < numPoints-1;i++){
      sampleRowsTest.row(i) = lowRankMatrix.row(rowIndexTest(i));
      approxRowsTest.row(i) = (W * K).row(rowIndexTest(i))*V.transpose();
      sampleColsTest.col(i) = lowRankMatrix.col(colIndexTest(i));
      approxColsTest.col(i) = (W * K)*(V.row(colIndexTest(i)).transpose());	
    }
    
    Eigen::MatrixXd sampleColsTestBlock = sampleColsTest.block(numPoints-1,0,nRows-numPoints+1,numPoints-1);
    Eigen::MatrixXd approxColsTestBlock = approxColsTest.block(numPoints-1,0,nRows-numPoints+1,numPoints-1);
    absFrobNormDiff = (sampleRowsTest-approxRowsTest).norm()+(sampleColsTestBlock-approxColsTestBlock).norm();
    relFrobNormDiff = absFrobNormDiff/(sampleRowsTest.norm()+sampleColsTestBlock.norm());
    approxError = relFrobNormDiff*(sqrt((nRows*nCols)/((numPoints-1)*(nCols+nRows-numPoints+1))));
    numPoints *= 1.5;
  }
  calculatedRank = W.cols(); 
}


void HODLR_Matrix::extractRowsCols(Eigen::MatrixXd & W, Eigen::MatrixXd & K, Eigen::MatrixXd & V, const Eigen::MatrixXd &inputMatrix,const std::vector<int> & rowIndex,const std::vector<int> & colIndex)const{
	
  int nRowsSelect = rowIndex.size();
  int nColsSelect = colIndex.size();
  
  //double rankTolerance=max(inputMatrix.rows(),inputMatrix.cols())*1e-16;
  double rankTolerance  = 1e-10;
  Eigen::MatrixXd WTemp = Eigen::MatrixXd::Zero(inputMatrix.rows(),nColsSelect);
  Eigen::MatrixXd VTemp = Eigen::MatrixXd::Zero(inputMatrix.cols(),nRowsSelect);
  Eigen::MatrixXd KTemp = Eigen::MatrixXd::Zero(nRowsSelect,nColsSelect);
  
  //fill W
  for (int i = 0; i < nColsSelect; i++)
    WTemp.col(i) = inputMatrix.col(colIndex[i]);
  
  //fill V
  for (int i = 0; i < nRowsSelect; i++)
    VTemp.col(i) = inputMatrix.row(rowIndex[i]).transpose();
  
  
  //fill K
  for (int i = 0; i < nRowsSelect; i++)
    for (int j = 0; j < nColsSelect; j++)
      KTemp(i,j) = inputMatrix(rowIndex[i],colIndex[j]);
  

  Eigen::MatrixXd svdW,svdV,svdK;
  int KTempRank = SVD_LowRankApprox(KTemp, rankTolerance, &svdW, &svdV, &svdK);

  
  //calculate W and V
  W = WTemp * svdV;
  V = VTemp * svdW;
  
  //calculate K
  
  K = Eigen::MatrixXd::Zero(KTempRank,KTempRank);	
  for (int i = 0; i < KTempRank; i++)
    K(i,i) = 1/svdK(i,i);
  
}




double HODLR_Matrix::fullPivACA_LowRankApprox(Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int max_i, const int min_j, const int max_j, const double tolerance, int & calculatedRank,const int minRank)  {
  
  int nRows = max_i - min_i + 1;
  int nCols = max_j - min_j + 1;
  int maxRank = std::min(nRows,nCols);
  int numColsW = 2;
  int numColsV = 2;

  Eigen::MatrixXd tempW(nRows,numColsW);
  Eigen::MatrixXd tempV(nCols,numColsV);
  Eigen::VectorXd colMaxValues(nCols);
  Eigen::VectorXi colMaxIdx(nCols);
  Eigen::MatrixXd residualMatrix = matrixData.block(min_i,min_j,nRows,nCols);
  double origMatrixNorm = residualMatrix.norm();

  double epsilon = 1;
  int k = 0;
  while (((epsilon > tolerance) || (k < minRank)) && (k < maxRank)){
    
    if ( k == numColsW - 1){
      numColsW = 2 * numColsW;
      numColsV = 2 * numColsV;
      tempW.conservativeResize(Eigen::NoChange,numColsW);
      tempV.conservativeResize(Eigen::NoChange,numColsV);
    }

    // Find largest pivot in the residual matrix
   
    for (int i = 0; i < nCols; i++)
      colMaxValues(i) = residualMatrix.col(i).cwiseAbs().maxCoeff(&colMaxIdx(i));
    int currRowIdx,currColIdx;
    double absMaxValue = colMaxValues.maxCoeff(&currColIdx);
    currRowIdx = colMaxIdx(currColIdx);
    double maxValue = residualMatrix(currRowIdx,currColIdx);
    if (absMaxValue <= minValueACA){
      break;
    }
    double currPivot = 1/maxValue;
    // Write to W & V
    tempW.col(k) = currPivot * residualMatrix.col(currColIdx);
    tempV.col(k) = residualMatrix.row(currRowIdx);

    // Update residual matrix
    residualMatrix -= tempW.col(k) * tempV.col(k).transpose();
 
    // Calculate epsilon
    epsilon = residualMatrix.norm()/origMatrixNorm;
    
    // Set Values for next iteration
    k++;
  }
  calculatedRank = k;
  // Return zero for zero matrix
  if ( k == 0){
    W = Eigen::MatrixXd::Zero(nRows,1);
    V = Eigen::MatrixXd::Zero(nCols,1);
    calculatedRank = 1;
    return epsilon;
  }
  
  // Return the original matrix if rank is equal to matrix dimensions
  if (k >= maxRank - 1){
    // Return original matrix
    // Skinny matrix
    if (nCols <= nRows){
      W = matrixData.block(min_i,min_j,nRows,nCols);
      V = Eigen::MatrixXd::Identity(nCols,nCols);
      calculatedRank = nCols;
    }// Fat matrix      
    else {
      W = Eigen::MatrixXd::Identity(nRows,nRows);
      V = matrixData.block(min_i,min_j,nRows,nCols).transpose();
      calculatedRank = nRows;
    } 
    return epsilon;
  }
  
  W = tempW.leftCols(calculatedRank);
  V = tempV.leftCols(calculatedRank);

  return epsilon;
}


double HODLR_Matrix::partialPivACA_LowRankApprox(Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int max_i, const int min_j, const int max_j, const double tolerance, int & calculatedRank,const int minRank)  {
  
  int nRows = max_i - min_i + 1;
  int nCols = max_j - min_j + 1;
  int maxRank = std::min(nRows,nCols);
  int numColsW = 2;
  int numColsV = 2;

  Eigen::MatrixXd tempW(nRows,numColsW);
  Eigen::MatrixXd tempV(nCols,numColsV);

  Eigen::VectorXd residualRow,residualCol;
  std::vector<bool> chosenRows(nRows),chosenCols(nCols);
  for (int i = 0; i < nRows; i++)
    chosenRows[i] = false;
  for (int i = 0; i < nCols; i++)
    chosenCols[i] = false;
  
  double frobNormSq = 0;
  double frobNorm   = 0;   
  double epsilon    = 1;
  int currRowIndex  = 0;
  int currColIndex  = 0;
  int nextRowIndex  = 0;
  int k = 0;
  
  while (((epsilon > tolerance) || (k < minRank)) && (k < maxRank)){
    
    if ( k == numColsW - 1){
      numColsW = 2 * numColsW;
      numColsV = 2 * numColsV;
      tempW.conservativeResize(Eigen::NoChange,numColsW);
      tempV.conservativeResize(Eigen::NoChange,numColsV);
    }

    chosenRows[currRowIndex] = true;
    int globalCurrRowIdx = currRowIndex + min_i;
    Eigen::VectorXd currRow = matrixData.block(globalCurrRowIdx,min_j,1,nCols).transpose();

    // Update row of Residual
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(nCols);
    for (int l = 0; l < k; l++){
      sum += tempW(currRowIndex,l) * tempV.col(l);
    }
    residualRow = (currRow - sum);
    // Find Next Column
    int maxInd;
    Eigen::VectorXd absCurrRow = residualRow.cwiseAbs();
    double maxValue = absCurrRow.maxCoeff(&maxInd);
    if (maxValue <= minValueACA){
      currRowIndex = chooseNNZRowIndex(chosenRows);
      if (currRowIndex == -1)
	break;
      continue;
      absCurrRow = residualRow.cwiseAbs();
      maxValue = absCurrRow.maxCoeff(&maxInd);
      currColIndex = maxInd;
    }
    if (chosenCols[maxInd] == false){
      currColIndex = maxInd;
    }else{
      currColIndex = chooseNextRowCol(chosenCols,residualRow.cwiseAbs());
      if (currColIndex == -1)
	break;
    }
    
    // Update column of Residual
    chosenCols[currColIndex] = true;
    int globalCurrColIdx = currColIndex + min_j;
    double currPivot = 1/residualRow(currColIndex);
    Eigen::VectorXd currColumn = matrixData.block(min_i,globalCurrColIdx,nRows,1);

    sum = Eigen::VectorXd::Zero(nRows);
    for(int l = 0; l < k; l++){
      sum += tempV(currColIndex,l) * tempW.col(l);
    }
    residualCol = currColumn - sum;
    
    // Find Next Row
    Eigen::VectorXd absCurrCol = residualCol.cwiseAbs();
    maxValue = absCurrCol.maxCoeff(&maxInd);
    
    if (chosenRows[maxInd] == false){
      nextRowIndex = maxInd;
    }else{
      nextRowIndex = chooseNextRowCol(chosenRows,residualCol.cwiseAbs());
    }
    if (nextRowIndex == -1)
      break;

    // Write to W & V
    tempW.col(k) = currPivot * residualCol;
    tempV.col(k) = residualRow;
    
    // Update Approximation Matrix
    //approxMatrix += currPivot * residualCol * residualRow.transpose();
    
    // Update Frobenious Norm
    double sumNorm = 0;
    for(int j = 0; j < k; j++){
      sumNorm += ((currPivot * residualCol.transpose()) * tempW.col(j)) * ((tempV.col(j).transpose()) * residualRow);
    }
    frobNormSq += 2 * sumNorm + ((currPivot * residualCol).squaredNorm()) * residualRow.squaredNorm();
    frobNorm = sqrt(frobNormSq);
    // Calculate epsilon
    epsilon = (currPivot * residualCol).norm() * residualRow.norm()/frobNorm;
    
    // Set Values for next iteration
    currRowIndex = nextRowIndex;
    k++;
  }
  calculatedRank = k;
  // Return zero for zero matrix
  if ( k == 0){
    W = Eigen::MatrixXd::Zero(nRows,1);
    V = Eigen::MatrixXd::Zero(nCols,1);
    calculatedRank = 1;
    return epsilon;
  }
  
  // Return the original matrix if rank is equal to matrix dimensions
  if (k >= maxRank - 1){
    // Return original matrix
    // Skinny matrix
    if (nCols <= nRows){
      W = matrixData.block(min_i,min_j,nRows,nCols);
      V = Eigen::MatrixXd::Identity(nCols,nCols);
      calculatedRank = nCols;
    }// Fat matrix      
    else {
      W = Eigen::MatrixXd::Identity(nRows,nRows);
      V = matrixData.block(min_i,min_j,nRows,nCols).transpose();
      calculatedRank = nRows;
    } 
    return epsilon;
  }
  
  W = tempW.leftCols(calculatedRank);
  V = tempV.leftCols(calculatedRank);
  return epsilon;
}

int HODLR_Matrix::chooseNNZRowIndex(const std::vector<bool> &chosenRows) const{
  int n = chosenRows.size();
  for (int i = 0;i < n; i++){
    if (chosenRows[i] == false)
      return i;
  }
  return -1;
}


int HODLR_Matrix::chooseNextRowCol(const std::vector<bool> &chosenRowsCols, const Eigen::VectorXd &currColRow) const{
  int n = currColRow.rows();
  for (int index = 0; index < n; index++){
    if  ((chosenRowsCols[index] == false) && (currColRow(index) > minValueACA))
      return index;
  }
  return -1;
}

void HODLR_Matrix::LUDecompose(const Eigen::MatrixXd &inputMatrix,Eigen::MatrixXd &LU,Eigen::MatrixXd &P) const{
    
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(inputMatrix);
  LU = lu.matrixLU();
  P = lu.permutationP();
  return;	  
}

void HODLR_Matrix::reset_attributes(){
  LRStoredInTree = false;
  createdRecLUfactorTree = false;
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
  if ((matrixDataAvail == false) && (matrixDataAvail_Sp == false)){
    std::cout<<"Error! Matrix data has been deleted from memory!"<<std::endl;
    exit(EXIT_FAILURE);
  }
  if (LR_Tolerance != input_tolerance){
    LR_Tolerance = input_tolerance;
    reset_attributes();
  }
}  

void HODLR_Matrix::set_MinValueACA(double input_minValueACA){
  if (matrixDataAvail == false){
    std::cout<<"Error! Matrix data has been deleted from memory!"<<std::endl;
    exit(EXIT_FAILURE);
  }
  if (minValueACA != input_minValueACA){
    minValueACA = input_minValueACA;
    reset_attributes();
  }
}

void HODLR_Matrix::set_def_LRMethod(std::string input_LRMethod){
  if (matrixDataAvail == false){
    std::cout<<"Error! Matrix data has been deleted from memory!"<<std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (indexTree.get_def_LRMethod() != input_LRMethod){
    indexTree.set_def_LRMethod(input_LRMethod);
    reset_attributes();
  }
}

void HODLR_Matrix::set_FreeMatrixMemory(bool inputVal){
  freeMatrixMemory = inputVal;
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

double HODLR_Matrix::get_recLU_FactorizationTime() const{
  return recLU_FactorizationTime;
}

double HODLR_Matrix::get_recLU_SolveTime() const{
  return recLU_SolveTime;
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

double HODLR_Matrix::get_LR_ComputationTime() const{
  return LR_ComputationTime;
}

double HODLR_Matrix::get_iter_SolveTime() const {
  return iter_SolveTime;
}

Eigen::MatrixXd HODLR_Matrix::get_Block(int min_i,int min_j,int numRows,int numCols){
  if (matrixDataAvail == true){
    return matrixData.block(min_i,min_j,numRows,numCols);
  }else  if(matrixDataAvail_Sp == true){
    Eigen::MatrixXd result(matrixData_Sp.block(min_i,min_j,numRows,numCols));
    return result;
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
    fill_BlockWithLRProduct(blkMatrix,LR_Min_i,LR_Min_j,LR_numRows,LR_numCols,root->topOffDiagU,root->topOffDiagK,root->topOffDiagV,blk_Min_i,blk_Min_j);
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
    fill_BlockWithLRProduct(blkMatrix,LR_Min_i,LR_Min_j,LR_numRows,LR_numCols,root->bottOffDiagU,root->bottOffDiagK,root->bottOffDiagV,blk_Min_i,blk_Min_j);
    }
  
  // Find Parts corresponding to diagonal blocks
  fill_Block(blkMatrix,root->right,min_i,min_j,max_i,max_j);
  fill_Block(blkMatrix,root->left,min_i,min_j,max_i,max_j);
  
}

void HODLR_Matrix::fill_BlockWithLRProduct(Eigen::MatrixXd & blkMatrix,int LR_Min_i,int LR_Min_j, int LR_numRows, int LR_numCols,Eigen::MatrixXd & LR_U,Eigen::MatrixXd & LR_K,Eigen::MatrixXd & LR_V,int blk_Min_i,int blk_Min_j){
  blkMatrix.block(blk_Min_i,blk_Min_j,LR_numRows,LR_numCols) = LR_U.block(LR_Min_i,0,LR_numRows,LR_U.cols()) * LR_K * LR_V.block(LR_Min_j,0,LR_numCols,LR_V.cols()).transpose();
}

/**************************************Extend-Add Functions************************************/

void HODLR_Matrix::extendAddUpdate(std::vector<int> & parentIdxVec,std::vector<Eigen::MatrixXd*> & LR_Update_U_PtrVec,std::vector<Eigen::MatrixXd*> & LR_Update_V_PtrVec,std::vector<std::vector<int>* > &updateIdxPtrVec,int sumChildRanks){
  int numChildren = LR_Update_U_PtrVec.size();
  Eigen::MatrixXd updateExtendU = Eigen::MatrixXd::Zero(matrixNumRows,sumChildRanks);
  Eigen::MatrixXd updateExtendV = Eigen::MatrixXd::Zero(matrixNumCols,sumChildRanks);
  int j_ind = 0;  
  for (int i = 0; i < numChildren; i++){
    //create extended U and V
    int currRank = (LR_Update_U_PtrVec[i])->cols(); 
    int updateMatrixSize = (updateIdxPtrVec[i])->size();
    // Find update matrix extend add indices
    std::vector<int> childUpdateExtendVec(updateMatrixSize);
    for (int j = 0; j < updateMatrixSize; j++){
      std::vector<int>::iterator iter;
      iter = std::lower_bound(parentIdxVec.begin(),parentIdxVec.end(),(*(updateIdxPtrVec[i]))[j]);
      int extendPos = iter - parentIdxVec.begin();
      childUpdateExtendVec[j] = extendPos;
    }
     // Go over all rows and columns in the update matrix
    for (int j = 0; j < updateMatrixSize; j++){
      for (int k = 0; k < currRank; k++){
	int rowIdx = childUpdateExtendVec[j];
	updateExtendU(rowIdx,k + j_ind) = (*(LR_Update_U_PtrVec[i]))(j,k);	
	updateExtendV(rowIdx,k + j_ind) = (*(LR_Update_V_PtrVec[i]))(j,k);	
      }
    }
    j_ind += currRank;
  }
  extendAddLRinTree(indexTree.rootNode,updateExtendU,updateExtendV,sumChildRanks);
  
}

void HODLR_Matrix::extendAddLRinTree(HODLR_Tree::node* HODLR_Root,const Eigen::MatrixXd & updateExtendU,const Eigen::MatrixXd & updateExtendV,int sumChildRanks){
  if (HODLR_Root->isLeaf == true){
    int numRows = HODLR_Root->max_i - HODLR_Root->min_i + 1;
    int numCols = HODLR_Root->max_j - HODLR_Root->min_j + 1;  
    HODLR_Root->leafMatrix += (updateExtendU).block(HODLR_Root->min_i,0,numRows,sumChildRanks) * (updateExtendV).block(HODLR_Root->min_j,0,numCols,sumChildRanks).transpose();       
    return;
  }

  int numRows_TopOffDiag  = HODLR_Root->splitIndex_i - HODLR_Root->min_i + 1; 
  int numRows_BottOffDiag = HODLR_Root->max_i - HODLR_Root->splitIndex_i;
  int numCols_TopOffDiag  = numRows_BottOffDiag;
  int numCols_BottOffDiag = numRows_TopOffDiag; 
  Eigen::MatrixXd U2_TopOffDiag,U2_BottOffDiag;
  Eigen::MatrixXd V2_TopOffDiag,V2_BottOffDiag;

  // Create topDiag U2s
  U2_TopOffDiag  = (updateExtendU).block(HODLR_Root->min_i,0,numRows_TopOffDiag,sumChildRanks);
  U2_BottOffDiag = (updateExtendU).block(HODLR_Root->min_i + numRows_TopOffDiag,0,numRows_BottOffDiag,sumChildRanks);

  // Create V2s
  V2_TopOffDiag  = (updateExtendV).block(HODLR_Root->min_j,0,numCols_TopOffDiag,sumChildRanks);
  V2_BottOffDiag = (updateExtendV).block(HODLR_Root->min_j + numCols_TopOffDiag,0,numCols_BottOffDiag,sumChildRanks);

  // Update current LR
  HODLR_Root->topOffDiagRank  = add_LR(HODLR_Root->topOffDiagU,HODLR_Root->topOffDiagK,HODLR_Root->topOffDiagV,HODLR_Root->topOffDiagU * HODLR_Root->topOffDiagK,HODLR_Root->topOffDiagV,U2_TopOffDiag,V2_TopOffDiag);
  HODLR_Root->bottOffDiagRank = add_LR(HODLR_Root->bottOffDiagU,HODLR_Root->bottOffDiagK,HODLR_Root->bottOffDiagV,HODLR_Root->bottOffDiagU * HODLR_Root->bottOffDiagK,HODLR_Root->bottOffDiagV,U2_BottOffDiag,V2_BottOffDiag);
   
  // Do the same for children
  extendAddLRinTree(HODLR_Root->left ,updateExtendU,updateExtendV,sumChildRanks);
  extendAddLRinTree(HODLR_Root->right,updateExtendU,updateExtendV,sumChildRanks);

}


int HODLR_Matrix::add_LR(Eigen::MatrixXd & result_U,Eigen::MatrixXd & result_K,Eigen::MatrixXd result_V,const Eigen::MatrixXd & U1, const Eigen::MatrixXd & V1, const Eigen::MatrixXd & U2, const Eigen::MatrixXd & V2){
  assert(U1.rows() == U2.rows());
  assert(V1.rows() == V2.rows());
  Eigen::MatrixXd Utot(U1.rows(),U1.cols() + U2.cols());
  Eigen::MatrixXd Vtot(V1.rows(),V1.cols() + V2.cols());
  Utot.topLeftCorner(U1.rows(),U1.cols())  = U1;
  Utot.topRightCorner(U2.rows(),U2.cols()) = U2;
  Vtot.topLeftCorner(V1.rows(),V1.cols())  = V1;
  Vtot.topRightCorner(V2.rows(),V2.cols()) = V2;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Utot);
  Eigen::MatrixXd thinQ;
  thinQ.setIdentity(Utot.rows(),Utot.cols());
  qr.householderQ().applyThisOnTheLeft(thinQ);
  int rank = qr.rank();
  Eigen::MatrixXd Q = thinQ.leftCols(rank);
  Eigen::MatrixXd sigma = (Q.transpose() * Utot) * (Vtot.transpose() * Q);
  Eigen::MatrixXd sigma_W,sigma_V,sigma_K;
  SVD_LowRankApprox(sigma,LR_Tolerance,&sigma_W,&sigma_K,&sigma_V);
  result_U = Q * sigma_W;
  result_K = sigma_K;
  result_V = Q * sigma_V;
  return sigma_K.rows();
}
