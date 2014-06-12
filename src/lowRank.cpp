#include "lowRank.hpp"

const double pi = 3.14159265359;

// Some needed function prototypes
// ACA
int chooseNNZRowIndex(const std::vector<bool> &chosenRows);
int chooseNextRowCol(const std::vector<bool> &chosenRowsCols, const Eigen::VectorXd &currColRow,const int minPivot);

// Boundary Identifier
void identifyBoundary(const Eigen::SparseMatrix<double> & inputGraph,const std::set<int> &rowSet,const std::set<int> &colSet,std::map<int,std::vector<int> > & rowPos,std::map<int,std::vector<int> > & colPos,int depth = -1);
void createIdxFromBoundaryMap(std::map<int,std::vector<int> > & rowPos, std::map<int,std::vector<int> > & colPos, int depth,std::vector<int> &rowIdx,std::vector<int> &colIdx);

// Psuedo Skeleton
void extractRowsCols(Eigen::MatrixXd & W, Eigen::MatrixXd & K, Eigen::MatrixXd & V, const Eigen::MatrixXd &inputMatrix,const std::vector<int> & rowIndex,const std::vector<int> & colIndex);
void extractRowsCols(const Eigen::MatrixXd & matrixData, int min_i,int min_j,int numRows,int numCols,Eigen::MatrixXd & W, Eigen::MatrixXd & K, Eigen::MatrixXd & V,const std::vector<int> & rowIndex,const std::vector<int> & colIndex);


double fullPivACA_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank,const int minPivot){

  int maxRank = std::min(numRows,numCols);
  int numColsW = 2;
  int numColsV = 2;

  Eigen::MatrixXd tempW(numRows,numColsW);
  Eigen::MatrixXd tempV(numCols,numColsV);
  Eigen::VectorXd colMaxValues(numCols);
  Eigen::VectorXi colMaxIdx(numCols);
  Eigen::MatrixXd residualMatrix = matrixData.block(min_i,min_j,numRows,numCols);
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
   
    for (int i = 0; i < numCols; i++)
      colMaxValues(i) = residualMatrix.col(i).cwiseAbs().maxCoeff(&colMaxIdx(i));
    int currRowIdx,currColIdx;
    double absMaxValue = colMaxValues.maxCoeff(&currColIdx);
    currRowIdx = colMaxIdx(currColIdx);
    double maxValue = residualMatrix(currRowIdx,currColIdx);
    if (absMaxValue <= minPivot){
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
    W = Eigen::MatrixXd::Zero(numRows,1);
    V = Eigen::MatrixXd::Zero(numCols,1);
    calculatedRank = 1;
    return epsilon;
  }
  
  // Return the original matrix if rank is equal to matrix dimensions
  if (k >= maxRank - 1){
    // Return original matrix
    // Skinny matrix
    if (numCols <= numRows){
      W = matrixData.block(min_i,min_j,numRows,numCols);
      V = Eigen::MatrixXd::Identity(numCols,numCols);
      calculatedRank = numCols;
    }// Fat matrix      
    else {
      W = Eigen::MatrixXd::Identity(numRows,numRows);
      V = matrixData.block(min_i,min_j,numRows,numCols).transpose();
      calculatedRank = numRows;
    } 
    return epsilon;
  }
  
  W = tempW.leftCols(calculatedRank);
  V = tempV.leftCols(calculatedRank);

  return epsilon;
}

double partialPivACA_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank,const int minPivot)  {
  
  int maxRank = std::min(numRows,numCols);
  int numColsW = 2;
  int numColsV = 2;

  Eigen::MatrixXd tempW(numRows,numColsW);
  Eigen::MatrixXd tempV(numCols,numColsV);

  Eigen::VectorXd residualRow,residualCol;
  std::vector<bool> chosenRows(numRows),chosenCols(numCols);
  for (int i = 0; i < numRows; i++)
    chosenRows[i] = false;
  for (int i = 0; i < numCols; i++)
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
    Eigen::VectorXd currRow = matrixData.block(globalCurrRowIdx,min_j,1,numCols).transpose();

    // Update row of Residual
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(numCols);
    for (int l = 0; l < k; l++){
      sum += tempW(currRowIndex,l) * tempV.col(l);
    }
    residualRow = (currRow - sum);
    // Find Next Column
    int maxInd;
    Eigen::VectorXd absCurrRow = residualRow.cwiseAbs();
    double maxValue = absCurrRow.maxCoeff(&maxInd);
    if (maxValue <= minPivot){
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
      currColIndex = chooseNextRowCol(chosenCols,residualRow.cwiseAbs(),minPivot);
      if (currColIndex == -1)
	break;
    }
    
    // Update column of Residual
    chosenCols[currColIndex] = true;
    int globalCurrColIdx = currColIndex + min_j;
    double currPivot = 1/residualRow(currColIndex);
    Eigen::VectorXd currColumn = matrixData.block(min_i,globalCurrColIdx,numRows,1);

    sum = Eigen::VectorXd::Zero(numRows);
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
      nextRowIndex = chooseNextRowCol(chosenRows,residualCol.cwiseAbs(),minPivot);
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
    W = Eigen::MatrixXd::Zero(numRows,1);
    V = Eigen::MatrixXd::Zero(numCols,1);
    calculatedRank = 1;
    return epsilon;
  }
  
  // Return the original matrix if rank is equal to matrix dimensions
  if (k >= maxRank - 1){
    // Return original matrix
    // Skinny matrix
    if (numCols <= numRows){
      W = matrixData.block(min_i,min_j,numRows,numCols);
      V = Eigen::MatrixXd::Identity(numCols,numCols);
      calculatedRank = numCols;
    }// Fat matrix      
    else {
      W = Eigen::MatrixXd::Identity(numRows,numRows);
      V = matrixData.block(min_i,min_j,numRows,numCols).transpose();
      calculatedRank = numRows;
    } 
    return epsilon;
  }
  
  W = tempW.leftCols(calculatedRank);
  V = tempV.leftCols(calculatedRank);
  return epsilon;
}

int chooseNNZRowIndex(const std::vector<bool> &chosenRows){
  int n = chosenRows.size();
  for (int i = 0;i < n; i++){
    if (chosenRows[i] == false)
      return i;
  }
  return -1;
}

int chooseNextRowCol(const std::vector<bool> &chosenRowsCols, const Eigen::VectorXd &currColRow, const int minPivot) {
  int n = currColRow.rows();
  for (int index = 0; index < n; index++){
    if  ((chosenRowsCols[index] == false) && (currColRow(index) > minPivot))
      return index;
  }
  return -1;
}



void PS_LowRankApprox_Sp(const Eigen::SparseMatrix<double> & matrixData_Sp,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j,const int numRows, const int numCols, const double tolerance, int &calculatedRank){
	

  std::set<int> rowIdxSet,colIdxSet;
 
  Eigen::SparseMatrix<double> lowRankMatrix_Sp = matrixData_Sp.block(min_i,min_j,numRows,numCols);
  //find numPoints
  if (lowRankMatrix_Sp.nonZeros() == 0){
    calculatedRank = 1;
    W = Eigen::MatrixXd::Zero(numRows,1);
    K = Eigen::MatrixXd::Zero(1,1);
    V = Eigen::MatrixXd::Zero(numCols,1);
    return;
  }
  for (int k = 0; k < lowRankMatrix_Sp.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(lowRankMatrix_Sp,k); it; ++it){
      rowIdxSet.insert(it.row());
      colIdxSet.insert(it.col());
    }
  std::vector<int> rowIndex(rowIdxSet.begin(),rowIdxSet.end());
  std::vector<int> colIndex(colIdxSet.begin(),colIdxSet.end());
 
  //choose rows and columns and do the low-rank approximation
  Eigen::MatrixXd dummyW,dummyK,dummyV;
  extractRowsCols(W,K,V,Eigen::MatrixXd(lowRankMatrix_Sp),rowIndex,colIndex);
  //extractRowsCols(Eigen::MatrixXd(lowRankMatrix_Sp),0,0,numRows,numCols,W,K,V,rowIndex,colIndex);
  calculatedRank = W.cols();

}

void PS_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j,const int numRows, const int numCols, const double tolerance, int &calculatedRank, const std::string pointChoosingMethod,const int minRank){
	
  
  int maxRank = std::min(numRows,numCols);
  int numPoints;

  if (minRank > 0)
    numPoints = minRank;
  else
    numPoints = maxRank/25 + 1;

  double absFrobNormDiff = 1;
  double relFrobNormDiff = 1;
  double approxError = 1;
  Eigen::MatrixXd lowRankMatrix = matrixData.block(min_i,min_j,numRows,numCols);
  
  while ((approxError > tolerance) && (numPoints <= maxRank)){	
    
    //create row and col index using Chebyshev or uniform nodes
    std::vector<int> rowIndex_Vec(numPoints),colIndex_Vec(numPoints);
    Eigen::VectorXi  rowIndex(numPoints),colIndex(numPoints);
    if (pointChoosingMethod == "Chebyshev")
      for (int i = 0 ; i < numPoints; i++){
	rowIndex(i) = floor((numRows + numRows * cos(pi * (2 * i + 1)/(2 * numPoints))) / 2);
	colIndex(i) = floor((numCols + numCols * cos(pi * (2 * i + 1)/(2 * numPoints))) / 2);
	rowIndex_Vec[i] = rowIndex(i);
	colIndex_Vec[i] = colIndex(i);
      }
    if (pointChoosingMethod == "Uniform" && numPoints > 1){
      int rowStride = numRows / (numPoints - 1);
      int colStride = numCols / (numPoints - 1);
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
    //extractRowsCols(W,K,V,lowRankMatrix,rowIndex_Vec,colIndex_Vec);
    extractRowsCols(matrixData,min_i,min_j,numRows,numCols,W,K,V,rowIndex_Vec,colIndex_Vec);
    

    calculatedRank = W.cols();

    //obtain stopping criterion
    Eigen::VectorXi rowIndexTest = (rowIndex.head(numPoints-1)+rowIndex.tail(numPoints-1))/2;
    Eigen::VectorXi colIndexTest = (colIndex.head(numPoints-1)+colIndex.tail(numPoints-1))/2;
    
    Eigen::MatrixXd sampleColsTest(numRows,numPoints-1);
    Eigen::MatrixXd approxColsTest(numRows,numPoints-1);
    Eigen::MatrixXd sampleRowsTest(numPoints-1,numCols);
    Eigen::MatrixXd approxRowsTest(numPoints-1,numCols);
    
    //fill KTempApprox
    for (int i = 0; i < numPoints-1;i++){
      sampleRowsTest.row(i) = lowRankMatrix.row(rowIndexTest(i));
      approxRowsTest.row(i) = (W * K).row(rowIndexTest(i))*V.transpose();
      sampleColsTest.col(i) = lowRankMatrix.col(colIndexTest(i));
      approxColsTest.col(i) = (W * K)*(V.row(colIndexTest(i)).transpose());	
    }
    
    Eigen::MatrixXd sampleColsTestBlock = sampleColsTest.block(numPoints-1,0,numRows-numPoints+1,numPoints-1);
    Eigen::MatrixXd approxColsTestBlock = approxColsTest.block(numPoints-1,0,numRows-numPoints+1,numPoints-1);
    absFrobNormDiff = (sampleRowsTest-approxRowsTest).norm()+(sampleColsTestBlock-approxColsTestBlock).norm();
    relFrobNormDiff = absFrobNormDiff/(sampleRowsTest.norm()+sampleColsTestBlock.norm());
    approxError = relFrobNormDiff*(sqrt((numRows*numCols)/((numPoints-1)*(numCols+numRows-numPoints+1))));
    numPoints *= 1.5;
  }
  calculatedRank = W.cols(); 
}


void extractRowsCols(Eigen::MatrixXd & W, Eigen::MatrixXd & K, Eigen::MatrixXd & V, const Eigen::MatrixXd &inputMatrix,const std::vector<int> & rowIndex,const std::vector<int> & colIndex){
	
  int numRowsSelect = rowIndex.size();
  int numColsSelect = colIndex.size();
  int numPoints = std::max(numRowsSelect,numColsSelect);
  
  //double rankTolerance=max(inputMatrix.rows(),inputMatrix.cols())*1e-16;
  
  double rankTolerance  = 1e-10;
  Eigen::MatrixXd WTemp = Eigen::MatrixXd::Zero(inputMatrix.rows(),numColsSelect);
  Eigen::MatrixXd VTemp = Eigen::MatrixXd::Zero(inputMatrix.cols(),numRowsSelect);
  Eigen::MatrixXd KTemp = Eigen::MatrixXd::Zero(numRowsSelect,numColsSelect);
  
  //fill W
  for (int i = 0; i < numColsSelect; i++)
    WTemp.col(i) = inputMatrix.col(colIndex[i]);
 

  //fill V
  for (int i = 0; i < numRowsSelect; i++)
    VTemp.col(i) = inputMatrix.row(rowIndex[i]).transpose();
  
  
  //fill K
  for (int i = 0; i < numRowsSelect; i++)
    for (int j = 0; j < numColsSelect; j++)
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
  
 
  /*
  Eigen::MatrixXd tempW = Eigen::MatrixXd::Zero(inputMatrix.rows(),numPoints);
  Eigen::MatrixXd tempV = Eigen::MatrixXd::Zero(inputMatrix.cols(),numPoints);
  Eigen::MatrixXd tempK = Eigen::MatrixXd::Zero(numPoints,numPoints);
  
  //fill W
  for (int i = 0; i < numPoints; i++){
    if (i < numColsSelect)
      tempW.col(i) = inputMatrix.col(colIndex[i]);
    if (i < numRowsSelect)
      tempV.col(i) = inputMatrix.row(rowIndex[i]).transpose();
  }
 
  //fill K
  for (int i = 0; i < numPoints; i++)
    for (int j = 0; j < numPoints; j++)
      if (i < numRowsSelect && j < numColsSelect)
	tempK(i,j) = inputMatrix(rowIndex[i],colIndex[j]);
  
  Eigen::FullPivLU<Eigen::MatrixXd> lu(tempK);
  lu.setThreshold(1e-10);
  int rank = lu.rank();

  V = ((lu.permutationP() * tempV.transpose()).transpose()).leftCols(rank);
  Eigen::MatrixXd L_Soln = lu.matrixLU().topLeftCorner(rank,rank).triangularView<Eigen::UnitLower>().solve(V.transpose());
  V = lu.matrixLU().topLeftCorner(rank,rank).triangularView<Eigen::Upper>().solve(L_Soln).transpose();
  K = Eigen::MatrixXd::Identity(rank,rank);
  W = (tempW * lu.permutationQ()).leftCols(rank);
  */
}


void extractRowsCols(const Eigen::MatrixXd & matrixData, int min_i,int min_j,int numRows,int numCols,Eigen::MatrixXd & W, Eigen::MatrixXd & K, Eigen::MatrixXd & V,const std::vector<int> & rowIndex,const std::vector<int> & colIndex){
	
  int numRowsSelect = rowIndex.size();
  int numColsSelect = colIndex.size();
  int numPoints = std::max(numRowsSelect,numColsSelect);

  double rankTolerance  = 1e-10;
  Eigen::MatrixXd tempW = Eigen::MatrixXd::Zero(numRows,numPoints);
  Eigen::MatrixXd tempV = Eigen::MatrixXd::Zero(numCols,numPoints);
  Eigen::MatrixXd tempK = Eigen::MatrixXd::Zero(numPoints,numPoints);
  
  //fill W
  for (int i = 0; i < numPoints; i++){
    if (i < numColsSelect)
      tempW.col(i) = matrixData.block(min_i,min_j + colIndex[i],numRows,1);
    if (i < numRowsSelect)
      tempV.col(i) = matrixData.block(min_i + rowIndex[i],min_j,1,numCols).transpose();
  
  }
 
  //fill K
  for (int i = 0; i < numPoints; i++)
    for (int j = 0; j < numPoints; j++)
      if (i < numRowsSelect && j < numColsSelect)
	tempK(i,j) = matrixData(min_i + rowIndex[i],min_j + colIndex[j]);
  
  Eigen::FullPivLU<Eigen::MatrixXd> lu(tempK);
  lu.setThreshold(1e-1);
  int rank = lu.rank();

  V = ((lu.permutationP() * tempV.transpose()).transpose()).leftCols(rank);
  Eigen::MatrixXd L_Soln = lu.matrixLU().topLeftCorner(rank,rank).triangularView<Eigen::UnitLower>().solve(V.transpose());
  V = lu.matrixLU().topLeftCorner(rank,rank).triangularView<Eigen::Upper>().solve(L_Soln).transpose();
  K = Eigen::MatrixXd::Identity(rank,rank);
  W = (tempW * lu.permutationQ()).leftCols(rank);
}



int SVD_LowRankApprox(const Eigen::MatrixXd & inputMatrix, const double accuracy, Eigen::MatrixXd* Wptr, Eigen::MatrixXd* Vptr, Eigen::MatrixXd* Kptr, int minRank){
  
  // Use svd to calculate the low-rank approximation
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(inputMatrix,Eigen::ComputeThinU|Eigen::ComputeThinV);
  
  // Calculate rank
  Eigen::VectorXd singularValues = svd.singularValues();
  int nSingularValues = singularValues.rows();
  int rank = nSingularValues;
  for (int i = 0; i < nSingularValues; i++)
    if (singularValues(i) < accuracy){
      rank = i ;
      break;
    } 
  if (rank < minRank)
    rank = minRank;
  if (rank !=0){
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
  }else{
    if (Wptr != NULL){
      *Wptr = Eigen::MatrixXd::Zero(inputMatrix.rows(),1);
    }
    if (Vptr != NULL){
      *Vptr = Eigen::MatrixXd::Zero(inputMatrix.cols(),1);
    }
    if (Kptr != NULL){
      *Kptr = Eigen::MatrixXd::Identity(1,1);
    }
    return 1;
  }
}


void SVD_LowRankApprox(const Eigen::MatrixXd & matrixData, Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols,const double tolerance, int & calculatedRank, const int minRank ){

  Eigen::MatrixXd lowRankMatrix = matrixData.block(min_i,min_j,numRows,numCols);
  calculatedRank = SVD_LowRankApprox(lowRankMatrix, tolerance, &W, &V, &K, minRank);
   
}

void identifyBoundary(const Eigen::SparseMatrix<double> & inputGraph,const std::set<int> &rowSet,const std::set<int> &colSet,std::map<int,std::vector<int> > & rowPos,std::map<int,std::vector<int> > & colPos,int depth){
  int numRows = rowSet.size();
  int numCols = colSet.size();
  assert(inputGraph.rows() == numRows + numCols);
  assert(inputGraph.cols() == numRows + numCols);
  
  std::vector<int> rowCurrClassVec;
  std::vector<int> colCurrClassVec;
  std::vector<int> rowNextClassVec;
  std::vector<int> colNextClassVec;
  std::map<int,bool> classifiedRows,classifiedCols;
  //initialize
  for (std::set<int>::iterator iter = rowSet.begin(); iter != rowSet.end(); ++ iter)
    classifiedRows[*iter] = false;
  
  for (std::set<int>::iterator iter = colSet.begin(); iter != colSet.end(); ++ iter)
    classifiedCols[*iter] = false;

  int numClassifiedRows = 0;
  int numClassifiedCols = 0;
  //Identify boundary nodes
  for (int k = 0; k < inputGraph.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(inputGraph,k); it; ++it){
      if (rowSet.count(it.row()) == 1 && colSet.count(it.col()) == 1){
	if (classifiedRows[it.row()] == false && classifiedCols[it.col()] == false){
	  rowPos[0].push_back(it.row());
	  colPos[0].push_back(it.col());
	  classifiedRows[it.row()] = true;
	  classifiedCols[it.col()] = true;
	  rowCurrClassVec.push_back(it.row());
	  colCurrClassVec.push_back(it.col());
	  numClassifiedRows ++;
	  numClassifiedCols ++;
	}
      }
    }     

  //Clasify other rows
  int rowCurrClass = 0;
  while (numClassifiedRows < numRows){
    if (rowCurrClass == depth)
      break;			 
    for (unsigned int i = 0; i < rowCurrClassVec.size();i++){
      Eigen::SparseMatrix<double> currNode = inputGraph.block(rowCurrClassVec[i],0,1,numRows+numCols);
      for (int k = 0; k < currNode.outerSize(); ++k)
	for (Eigen::SparseMatrix<double>::InnerIterator it(currNode,k); it; ++it){
	  if (rowSet.count(it.col()) == 1 && classifiedRows[it.col()] == false){
	    rowPos[rowCurrClass + 1].push_back(it.col());
	    classifiedRows[it.col()] = true;
	    rowNextClassVec.push_back(it.col());
	    numClassifiedRows ++;
	  }
	}  
    }
    rowCurrClass ++;
    rowCurrClassVec = rowNextClassVec;
    rowNextClassVec.clear();
  }

   //Clasify other cols
  int colCurrClass = 0;
  while (numClassifiedCols < numCols){
    if (colCurrClass == depth)
      break;
    for (unsigned int i = 0; i < colCurrClassVec.size();i++){
      Eigen::SparseMatrix<double> currNode = inputGraph.block(colCurrClassVec[i],0,1,numRows+numCols);
      for (int k = 0; k < currNode.outerSize(); ++k)
	for (Eigen::SparseMatrix<double>::InnerIterator it(currNode,k); it; ++it){
	  if (colSet.count(it.col()) == 1 && classifiedCols[it.col()] == false){
	    colPos[colCurrClass + 1].push_back(it.col());
	    classifiedCols[it.col()] = true;
	    colNextClassVec.push_back(it.col());
	    numClassifiedCols ++;
	  }
	} 
    }
    colCurrClass ++;
    colCurrClassVec = colNextClassVec;
    colNextClassVec.clear();
  }
 
}

void createIdxFromBoundaryMap( std::map<int,std::vector<int> > & rowPos, std::map<int,std::vector<int> > & colPos, int depth,std::vector<int> & rowIdx,std::vector<int> & colIdx){
  assert(depth >= 0 );
  rowIdx = rowPos[0];
  colIdx = colPos[0];
  for (int i = 1; i <= depth; i++){
    rowIdx.insert(rowIdx.end(),rowPos[i].begin(),rowPos[i].end());
    colIdx.insert(colIdx.end(),colPos[i].begin(),colPos[i].end());
  }

}

void PS_Boundary_LowRankApprox(const Eigen::MatrixXd & matrixData,const Eigen::SparseMatrix<double> graphData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols,int & calculatedRank, const int depth){
   std::map<int,std::vector<int> > rowPos,colPos;
   std::set<int> rowSet,colSet;
   std::vector<int> rowIdx,colIdx;
   int max_i     = min_i + numRows - 1;
   int max_j     = min_j + numCols - 1;
   int minIdx    = std::min(min_i,min_j);
   int maxIdx    = std::max(max_i,max_j);
   int offset_i  = min_i - minIdx;
   int offset_j  = min_j - minIdx;
   int numPoints = maxIdx - minIdx + 1;
   
   for (int i = 0; i < numRows; i++)
     rowSet.insert(i + min_i - minIdx);
    
   for (int i = 0; i < numCols; i++)
     colSet.insert(i + min_j - minIdx);
    
   identifyBoundary(graphData.block(minIdx,minIdx,numPoints,numPoints),rowSet,colSet,rowPos,colPos,depth);
   createIdxFromBoundaryMap(rowPos,colPos,depth,rowIdx,colIdx);

   // Adjust for offsets
   for (unsigned int i = 0; i < rowIdx.size();i++)
     rowIdx[i] -= offset_i;
   for (unsigned int i = 0; i < colIdx.size();i++)
     colIdx[i] -= offset_j;
   //extractRowsCols(W,K,V,matrixData.block(min_i,min_j,numRows,numCols),rowIdx,colIdx);
   extractRowsCols(matrixData,min_i,min_j,numRows,numCols,W,K,V,rowIdx,colIdx);
   calculatedRank = K.cols();
}

