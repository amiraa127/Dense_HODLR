#include "lowRank.hpp"

const double pi = 3.14159265359;

// Some needed function prototypes
// ACA
//int chooseNNZRowIndex(const std::vector<bool> &chosenRows);
int chooseNextRowCol(const std::vector<bool> &chosenRowsCols, const Eigen::VectorXd &currColRow,const int minPivot);

// Boundary Identifier
int identifyBoundary(const Eigen::SparseMatrix<double> & inputGraph,const std::set<int> &rowSet,const std::set<int> &colSet,std::map<int,std::vector<int> > & rowPos,std::map<int,std::vector<int> > & colPos,int depth = -1);
void createIdxFromBoundaryMap(std::map<int,std::vector<int> > & rowPos, std::map<int,std::vector<int> > & colPos, int depth,std::vector<int> &rowIdx,std::vector<int> &colIdx);

// Psuedo Skeleton
//void extractRowsCols(Eigen::MatrixXd & W, Eigen::MatrixXd & K, Eigen::MatrixXd & V, const Eigen::MatrixXd &inputMatrix,const std::vector<int> & rowIndex,const std::vector<int> & colIndex);


template <typename T>
double fullPivACA_LowRankApprox(const T & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank,const int minPivot){

  int maxRank = std::min(numRows,numCols);
  int numColsW = 2;
  int numColsV = 2;
 
  W = Eigen::MatrixXd(numRows,numColsW);
  V = Eigen::MatrixXd(numCols,numColsV);
  
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
      W.conservativeResize(Eigen::NoChange,numColsW);
      V.conservativeResize(Eigen::NoChange,numColsV);
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

    W.col(k) = currPivot * residualMatrix.col(currColIdx);
    V.col(k) = residualMatrix.row(currRowIdx);

    
    // Update residual matrix
    residualMatrix -= W.col(k) * V.col(k).transpose();
 
    
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

  W.conservativeResize(Eigen::NoChange,calculatedRank);
  V.conservativeResize(Eigen::NoChange,calculatedRank);
  
  return epsilon;
}

template <typename T>
double partialPivACA_LowRankApprox(const T & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank,const int maxRank,const int minPivot){
  
  if (maxRank == 0){
    W = Eigen::MatrixXd::Zero(numRows,1);
    V = Eigen::MatrixXd::Zero(numCols,1);
    calculatedRank = 1;
    return 1;
  }
 
  int rankUpperBound = std::min(numRows,numCols);
  int numColsW = 2;
  int numColsV = 2;

  Eigen::MatrixXd tempW(numRows,numColsW);
  Eigen::MatrixXd tempV(numCols,numColsV);

  Eigen::VectorXd residualRow,residualCol;
  std::vector<bool> chosenRows(numRows,false),chosenCols(numCols,false);
  //for (int i = 0; i < numRows; i++)
  //  chosenRows[i] = false;
  //for (int i = 0; i < numCols; i++)
  //  chosenCols[i] = false;
  
  double frobNormSq = 0;
  double frobNorm   = 0;   
  double epsilon    = 1;
  int currRowIndex  = 0;
  int currColIndex  = 0;
  int nextRowIndex  = 0;
  int k = 0;
  
  while (((epsilon > tolerance) || (k < minRank)) && (k < rankUpperBound)){
    
    // increase the size of W and V if limit is reached
    if ( k == numColsW - 1){
      numColsW = 2 * numColsW;
      numColsV = 2 * numColsV;
      tempW.conservativeResize(Eigen::NoChange,numColsW);
      tempV.conservativeResize(Eigen::NoChange,numColsV);
    }

    chosenRows[currRowIndex] = true;
    int globalCurrRowIdx    = currRowIndex + min_i;
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
   
    /*
    if (maxValue <= minPivot){
      currRowIndex = chooseNNZRowIndex(chosenRows);
      if (currRowIndex == -1)
	break;
      continue;
      //absCurrRow = residualRow.cwiseAbs();
      //maxValue = absCurrRow.maxCoeff(&maxInd);
      //currColIndex = maxInd;
    }
    */

    if (maxValue <= minPivot){
      std::vector<bool>::iterator nextRowIter;
      nextRowIter = std::find(chosenRows.begin(),chosenRows.end(),false);
      if (nextRowIter == chosenRows.end())
	break;
      currRowIndex = nextRowIter - chosenRows.begin();
      continue;

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
    if (k == maxRank)
      break;
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
  if (k >= rankUpperBound - 1){
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

/*
int chooseNNZRowIndex(const std::vector<bool> &chosenRows){
  int n = chosenRows.size();
  for (int i = 0;i < n; i++){
    if (chosenRows[i] == false)
      return i;
  }
  return -1;
}
*/
int chooseNextRowCol(const std::vector<bool> &chosenRowsCols, const Eigen::VectorXd &currColRow, const int minPivot) {
  int n = currColRow.rows();
  for (int index = 0; index < n; index++){
    if  ((chosenRowsCols[index] == false) && (currColRow(index) > minPivot))
      return index;
  }
  return -1;
}



void PS_LowRankApprox_Sp(const Eigen::SparseMatrix<double> & matrixData_Sp,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j,const int numRows, const int numCols, const double tolerance, int &calculatedRank){
	

  std::set<int> rowIdxSet,colIdxSet;
 
  Eigen::SparseMatrix<double> lowRankMatrix_Sp = matrixData_Sp.block(min_i,min_j,numRows,numCols);
  //find numPoints
  if (lowRankMatrix_Sp.nonZeros() == 0){
    calculatedRank = 1;
    W = Eigen::MatrixXd::Zero(numRows,1);
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
  extractRowsCols(Eigen::MatrixXd(lowRankMatrix_Sp),0,0,numRows,numCols,W,V,rowIndex,colIndex,1e-10,calculatedRank);
  

}

void PS_LowRankApprox(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j,const int numRows, const int numCols, const double tolerance, int &calculatedRank, const std::string pointChoosingMethod,const int minRank){
	
  
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
     extractRowsCols(matrixData,min_i,min_j,numRows,numCols,W,V,rowIndex_Vec,colIndex_Vec,1e-10,calculatedRank);
    
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
      approxRowsTest.row(i) = W.row(rowIndexTest(i))*V.transpose();
      sampleColsTest.col(i) = lowRankMatrix.col(colIndexTest(i));
      approxRowsTest.row(i) = W.row(rowIndexTest(i))*V.transpose();
     
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

template <typename T>
void extractRowsCols(const T & matrixData, int min_i,int min_j,int numRows,int numCols,Eigen::MatrixXd & W,Eigen::MatrixXd & V,const std::vector<int> & rowIndex,const std::vector<int> & colIndex,const double tolerance,int & calculatedRank,const std::string mode){

    int numRowsSelect = rowIndex.size();
    int numColsSelect = colIndex.size();
    if (mode == "fullPivLU") {
        int numPoints = std::max(numRowsSelect,numColsSelect);
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
        lu.setThreshold(tolerance);
        int rank = lu.rank();
       
        if (rank > 0){
            V = ((lu.permutationP() * tempV.transpose()).transpose()).leftCols(rank);
            Eigen::MatrixXd L_Soln = lu.matrixLU().topLeftCorner(rank,rank).triangularView<Eigen::UnitLower>().solve(V.transpose());
            V = lu.matrixLU().topLeftCorner(rank,rank).triangularView<Eigen::Upper>().solve(L_Soln).transpose();
            W = (tempW * lu.permutationQ()).leftCols(rank);
            calculatedRank = rank;
        }else{
            W = Eigen::MatrixXd::Zero(numRows,1);
            V = Eigen::MatrixXd::Zero(numCols,1);
            calculatedRank = 1;
        }
    }else if(mode == "SVD"){
        Eigen::MatrixXd WTemp = Eigen::MatrixXd::Zero(numRows,numColsSelect);
        Eigen::MatrixXd VTemp = Eigen::MatrixXd::Zero(numCols,numRowsSelect);
        Eigen::MatrixXd KTemp = Eigen::MatrixXd::Zero(numRowsSelect,numColsSelect);
        
        //fill W
        for (int i = 0; i < numColsSelect; i++)
            WTemp.col(i) = matrixData.block(min_i,min_j + colIndex[i],numRows,1);
        
        //fill V
        for (int i = 0; i < numRowsSelect; i++)
            VTemp.col(i) = matrixData.block(min_i + rowIndex[i],min_j,1,numCols).transpose();
        
        
        //fill K
        for (int i = 0; i < numRowsSelect; i++)
            for (int j = 0; j < numColsSelect; j++)
                KTemp(i,j) = matrixData(min_i + rowIndex[i],min_j + colIndex[j]);
        
        
        Eigen::MatrixXd svdW,svdV,svdK;
        calculatedRank = SVD_LowRankApprox(KTemp,tolerance, &svdW, &svdV, &svdK);
        
        
        //calculate K
	Eigen::MatrixXd K = Eigen::MatrixXd::Zero(calculatedRank,calculatedRank);
        for (int i = 0; i < calculatedRank; i++)
            K(i,i) = 1./svdK(i,i);
	  
        //calculate W and V
        W = WTemp * svdV * K;
        V = VTemp * svdW;
      
    }else{
        std::cout<<"Error! Unknown operation mode."<<std::endl;
        exit(EXIT_FAILURE);
    }
  }


int SVD_LowRankApprox(const Eigen::MatrixXd & matrixData, const double accuracy, Eigen::MatrixXd* Wptr, Eigen::MatrixXd* Vptr, Eigen::MatrixXd* Kptr, int minRank){
  
  // Use svd to calculate the low-rank approximation
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrixData,Eigen::ComputeThinU|Eigen::ComputeThinV);
  
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
      *Wptr = Eigen::MatrixXd::Zero(matrixData.rows(),1);
    }
    if (Vptr != NULL){
      *Vptr = Eigen::MatrixXd::Zero(matrixData.cols(),1);
    }
    if (Kptr != NULL){
      *Kptr = Eigen::MatrixXd::Identity(1,1);
    }
    return 1;
  }
}


template <typename T>
void SVD_LowRankApprox(const T & matrixData, Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols,const double tolerance, int & calculatedRank, const int minRank ){
  Eigen::MatrixXd lowRankMatrix = matrixData.block(min_i,min_j,numRows,numCols);
  calculatedRank = SVD_LowRankApprox(lowRankMatrix, tolerance, &W, &V, &K, minRank);
}

int identifyBoundary(const Eigen::SparseMatrix<double> & inputGraph,const std::set<int> &rowSet,const std::set<int> &colSet,std::map<int,std::vector<int> > & rowPos,std::map<int,std::vector<int> > & colPos,int maxDepth){
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
  if (numClassifiedRows == 0){
    if ((numRows >= 2 ) && (numCols >= 2 )){
      int numSel = 2;
      std::vector<int> rowVec(rowSet.begin(),rowSet.end());
      std::vector<int> colVec(colSet.begin(),colSet.end());
      std::vector<int> rowSel = createUniqueRndIdx(0,numRows-1,numSel);
      std::vector<int> colSel = createUniqueRndIdx(0,numCols-1,numSel);
      for (int i = 0; i < numSel ; i++){
	int rowIdx = rowVec[rowSel[i]];
	int colIdx = colVec[colSel[i]];
	rowPos[0].push_back(rowIdx);
	colPos[0].push_back(colIdx);
	classifiedRows[rowIdx] = true;
	classifiedCols[colIdx] = true;
	rowCurrClassVec.push_back(rowIdx);
	colCurrClassVec.push_back(colIdx);
	numClassifiedRows ++;
	numClassifiedCols ++;
      }
    }else
      return 1;
  }
  //Clasify other rows
  int rowCurrClass = 0;
  while (numClassifiedRows < numRows){
    if (rowCurrClass == maxDepth)
      break;			 
    int currClassification = 0;
    for (unsigned int i = 0; i < rowCurrClassVec.size();i++){
      Eigen::SparseMatrix<double> currNode = inputGraph.block(rowCurrClassVec[i],0,1,numRows+numCols);
      for (int k = 0; k < currNode.outerSize(); ++k)
	for (Eigen::SparseMatrix<double>::InnerIterator it(currNode,k); it; ++it){
	  if (rowSet.count(it.col()) == 1 && classifiedRows[it.col()] == false){
	    rowPos[rowCurrClass + 1].push_back(it.col());
	    classifiedRows[it.col()] = true;
	    rowNextClassVec.push_back(it.col());
	    numClassifiedRows ++;
	    currClassification ++;
	  }
	}  
    }
    if (currClassification == 0){
      // plant a seed 
      int numSeeds = rowPos[rowCurrClass].size();
      for (std::set<int>::iterator iter = rowSet.begin(); iter != rowSet.end(); ++ iter)
	if (classifiedRows[*iter] == false){
	  rowPos[rowCurrClass + 1].push_back(*iter);
	  classifiedRows[*iter] = true;
	  rowNextClassVec.push_back(*iter);
	  numClassifiedRows ++;
	  currClassification ++;
	  if (currClassification == numSeeds)
	    break;
	    }
    }
    rowCurrClass ++;
    rowCurrClassVec = rowNextClassVec;
    rowNextClassVec.clear();
  }
  
  //Clasify other cols
  int colCurrClass = 0;
  while (numClassifiedCols < numCols){
    if (colCurrClass == maxDepth)
      break;
    int currClassification = 0;
    for (unsigned int i = 0; i < colCurrClassVec.size();i++){
      Eigen::SparseMatrix<double> currNode = inputGraph.block(colCurrClassVec[i],0,1,numRows+numCols);
      for (int k = 0; k < currNode.outerSize(); ++k)
	for (Eigen::SparseMatrix<double>::InnerIterator it(currNode,k); it; ++it){
	  if (colSet.count(it.col()) == 1 && classifiedCols[it.col()] == false){
	    colPos[colCurrClass + 1].push_back(it.col());
	    classifiedCols[it.col()] = true;
	    colNextClassVec.push_back(it.col());
	    numClassifiedCols ++;
	    currClassification ++;
	  }
	} 
    }
    if (currClassification == 0){
      // plant a seed
      int numSeeds = colPos[colCurrClass].size();
      for (std::set<int>::iterator iter = colSet.begin(); iter != colSet.end(); ++ iter)
	if (classifiedCols[*iter] == false){
	  colPos[colCurrClass + 1].push_back(*iter);
	  classifiedCols[*iter] = true;
	  colNextClassVec.push_back(*iter);
	  numClassifiedCols ++;
	  currClassification ++;
	  if (currClassification == numSeeds)
	    break;
	    }
    }
    colCurrClass ++;
    colCurrClassVec = colNextClassVec;
    colNextClassVec.clear();
  }
  return 0;
}

void createIdxFromBoundaryMap( std::map<int,std::vector<int> > & rowPos, std::map<int,std::vector<int> > & colPos, int depth,std::vector<int> & rowIdx,std::vector<int> & colIdx){
  assert(depth >= 0 );
  rowIdx = rowPos[0];
  colIdx = colPos[0];
  
  for (int i = 1; i <= (std::min(depth,(int)rowPos.size() - 1)); i++)
    rowIdx.insert(rowIdx.end(),rowPos[i].begin(),rowPos[i].end());

  for (int i = 1; i <= (std::min(depth,(int)colPos.size() - 1)); i++)
    colIdx.insert(colIdx.end(),colPos[i].begin(),colPos[i].end());


}

template <typename T>
void PS_Boundary_LowRankApprox(const T & matrixData,const Eigen::SparseMatrix<double> graphData,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j, const int numRows, const int numCols,double tolerance,int & calculatedRank, const int depth, const std::string savePath){

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
  
  int noInteraction = identifyBoundary(graphData.block(minIdx,minIdx,numPoints,numPoints),rowSet,colSet,rowPos,colPos,depth);
  if (noInteraction == 1){
    W = Eigen::MatrixXd::Zero(numRows,1);
    V = Eigen::MatrixXd::Zero(numCols,1);
    calculatedRank = 1;
    return;
  }
  
  createIdxFromBoundaryMap(rowPos,colPos,depth,rowIdx,colIdx);
  
  // Adjust for offsets
  for (unsigned int i = 0; i < rowIdx.size();i++)
    rowIdx[i] -= offset_i;
  for (unsigned int i = 0; i < colIdx.size();i++)
    colIdx[i] -= offset_j;

  extractRowsCols(matrixData,min_i,min_j,numRows,numCols,W,V,rowIdx,colIdx,tolerance,calculatedRank);
  
  /*
    for(std::map<int,std::vector<int> >::iterator iter = rowPos.begin(); iter != rowPos.end(); ++iter){
    std::cout<<iter->first<<":";
    for (unsigned int i = 0; i < iter->second.size();i++)
    std::cout<<iter->second[i]<<" ";
    std::cout<<std::endl;
    }
  */
 
  if (savePath != "none"){
    int numRowsSelect = rowIdx.size();
    int numColsSelect = colIdx.size();
    numPoints = std::max(numRowsSelect,numColsSelect);
    Eigen::MatrixXd tempK = Eigen::MatrixXd::Zero(numPoints,numPoints);
    
    //fill K
     for (int i = 0; i < numPoints; i++)
       for (int j = 0; j < numPoints; j++)
	 if (i < numRowsSelect && j < numColsSelect)
	   tempK(i,j) = matrixData(min_i + rowIdx[i],min_j + colIdx[j]);
   
     Eigen::FullPivLU<Eigen::MatrixXd> lu(tempK);
     lu.setThreshold(tolerance);
     
     Eigen::VectorXi colIdxVec = Eigen::VectorXi::Zero(numPoints);
     Eigen::VectorXi rowIdxVec = Eigen::VectorXi::Zero(numPoints);
     
     for (int i = 0; i < colIdx.size(); i++)
       colIdxVec[i] = colIdx[i] + offset_j;
     for (int i = 0; i < rowIdx.size(); i++)
       rowIdxVec[i] = rowIdx[i] + offset_i;
     Eigen::VectorXi colIdxVecPerm = lu.permutationQ() * colIdxVec ;
     Eigen::VectorXi rowIdxVecPerm = lu.permutationP() * rowIdxVec ;
     
     std::vector<double> rowDist = std::vector<double>(numPoints,0);
     std::vector<double> colDist = std::vector<double>(numPoints,0);
     
     for (int i = 0; i <= depth; i++)
       std::sort(rowPos[i].begin(),rowPos[i].end());
     for (int i = 0; i < numPoints; i++)
       for (int j = 0; j <= depth; j++){
	 if (std::binary_search(rowPos[j].begin(),rowPos[j].end(),rowIdxVecPerm[i])){
           //std::cout<<i<<" "<<j<<std::endl;
	   rowDist[i] = j;
	   break;
	 }
       }
     
         
     for (int i = 0; i <= depth; i++)
       std::sort(colPos[i].begin(),colPos[i].end());
     for (int i = 0; i < numPoints; i++)
       for (int j = 0; j <= depth; j++){
	 if (std::binary_search(colPos[j].begin(),colPos[j].end(),colIdxVecPerm[i])){
           //std::cout<<i<<" "<<j<<std::endl;
	   colDist[i] = j;
	   break;
	 }
       }
     std::string rowFileName = savePath + "_Row";
     std::string colFileName = savePath + "_Col";
     saveVectorAsText(rowFileName,rowDist);
     saveVectorAsText(colFileName,colDist);
   }
} 
	   
   
int getBoundaryRowColIdx(const Eigen::SparseMatrix<double>  & graphData,const int min_i, const int min_j,const int numRows,const int numCols,const int depth,std::vector<int> & rowIdx,std::vector<int> & colIdx){
  std::map<int,std::vector<int> > rowPos,colPos;
  std::set<int> rowSet,colSet;
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
  
  int noInteraction = identifyBoundary(graphData.block(minIdx,minIdx,numPoints,numPoints),rowSet,colSet,rowPos,colPos,depth);
  if (noInteraction == 1)
    return 1;
  
  createIdxFromBoundaryMap(rowPos,colPos,depth,rowIdx,colIdx);
    
  // Adjust for offsets
  for (unsigned int i = 0; i < rowIdx.size();i++)
    rowIdx[i] -= offset_i;
  for (unsigned int i = 0; i < colIdx.size();i++)
    colIdx[i] -= offset_j;
  
  std::sort(rowIdx.begin(),rowIdx.end());
  std::sort(colIdx.begin(),colIdx.end());
 
 return 0;
}


//int add_LR(Eigen::MatrixXd & result_U,Eigen::MatrixXd & result_K,Eigen::MatrixXd & result_V,const Eigen::MatrixXd & U1, const Eigen::MatrixXd & V1, const Eigen::MatrixXd & U2, const Eigen::MatrixXd & V2,double tol,std::string mode){

int add_LR(Eigen::MatrixXd & result_U,Eigen::MatrixXd & result_V,const Eigen::MatrixXd & U1, const Eigen::MatrixXd & V1, const Eigen::MatrixXd & U2, const Eigen::MatrixXd & V2,double tol,std::string mode){

  assert(U1.rows() == U2.rows());
  assert(V1.rows() == V2.rows());
  Eigen::MatrixXd Utot(U1.rows(),U1.cols() + U2.cols());
  Eigen::MatrixXd Vtot(V1.rows(),V1.cols() + V2.cols());
  Utot.leftCols(U1.cols())  = U1;
  Utot.rightCols(U2.cols()) = U2;
  Vtot.leftCols(V1.cols())  = V1;
  Vtot.rightCols(V2.cols()) = V2;
  
  if (mode == "Compress_QR"){
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_U(Utot);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr_V(Vtot);
    Eigen::MatrixXd thinQ_U,thinQ_V;
    thinQ_U.setIdentity(Utot.rows(),Utot.cols());
    thinQ_V.setIdentity(Vtot.rows(),Vtot.cols());
    qr_U.householderQ().applyThisOnTheLeft(thinQ_U);
    qr_V.householderQ().applyThisOnTheLeft(thinQ_V);
    int rank_U = qr_U.rank();
    int rank_V = qr_V.rank();
    rank_U = std::max(rank_U,1);
    rank_V = std::max(rank_V,1);
    
    Eigen::MatrixXd Q_U = thinQ_U.leftCols(rank_U);
    Eigen::MatrixXd Q_V = thinQ_V.leftCols(rank_V);
    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> permMatrix_U = qr_U.colsPermutation();
    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> permMatrix_V = qr_V.colsPermutation();
    Eigen::MatrixXd sigma = (Q_U.transpose() * Utot * Vtot.transpose() * Q_V);
    Eigen::MatrixXd sigma_W,sigma_V,sigma_K;
    assert(sigma.rows() * sigma.cols() > 0);
    ::SVD_LowRankApprox(sigma,tol,&sigma_W,&sigma_V,&sigma_K);
    //result_U = Q_U * sigma_W;
    //result_K = sigma_K;
    result_U = Q_U * sigma_W * sigma_K;
    result_V = Q_V * sigma_V;
    return sigma_K.rows();
  }else if (mode == "Compress_LU"){
    Eigen::MatrixXd U_U,V_U;
    int rank_U;
    ::fullPivACA_LowRankApprox(Utot,U_U,V_U,0,0,Utot.rows(),Utot.cols(),tol,rank_U);
    Eigen::MatrixXd U_V,V_V;
    int rank_V;
    ::fullPivACA_LowRankApprox(Vtot,U_V,V_V,0,0,Vtot.rows(),Vtot.cols(),tol,rank_V);
    Eigen::MatrixXd sigma = V_U.transpose() * V_V;
    Eigen::MatrixXd sigma_W,sigma_V,sigma_K;
    ::SVD_LowRankApprox(sigma,tol,&sigma_W,&sigma_V,&sigma_K);
    //result_U = U_U * sigma_W;
    //result_K = sigma_K;
    result_U = U_U * sigma_W * sigma_K;
    result_V = U_V * sigma_V;
    return sigma_K.rows();
  }else if (mode == "Exact"){
    int totRank = U1.cols() + U2.cols();
    result_U = Utot;
    result_V = Vtot;
    //result_K = Eigen::MatrixXd::Identity(totRank,totRank);
    return totRank;
  }else{
    std::cout<<"Error! Unknown operation mode"<<std::endl;
    exit(EXIT_FAILURE);
  }
}


int PS_PseudoInverse(Eigen::MatrixXd & colMatrix,Eigen::MatrixXd & rowMatrix, Eigen::MatrixXd & U, Eigen::MatrixXd & V,std::vector<int> rowIdxVec,const  double tol,const std::string mode){
  int numRowsSelect = rowMatrix.cols();
  int numColsSelect = colMatrix.cols();
  
  Eigen::MatrixXd tempK = Eigen::MatrixXd::Zero(numRowsSelect,numColsSelect);
   
  for (int i = 0; i < numRowsSelect; i++)
    for (int j = 0; j < numColsSelect; j++)
      if (i < (int)rowIdxVec.size()) 
	tempK(i,j) = colMatrix(rowIdxVec[i],j);
  int rank;
  if (mode == "fullPivLU"){
    Eigen::FullPivLU<Eigen::MatrixXd> lu(tempK);
    lu.setThreshold(tol);
    rank = lu.rank();
    double largestPivot = std::abs((lu.matrixLU())(0,0));
    
    if ((rank > 0) && (largestPivot >= 1e-6)){
      V = ((lu.permutationP() * rowMatrix.transpose()).transpose()).leftCols(rank);
      Eigen::MatrixXd L_Soln = lu.matrixLU().topLeftCorner(rank,rank).triangularView<Eigen::UnitLower>().solve(V.transpose());
      V = lu.matrixLU().topLeftCorner(rank,rank).triangularView<Eigen::Upper>().solve(L_Soln).transpose();
      U = (colMatrix * lu.permutationQ()).leftCols(rank);
    }else{
      U = Eigen::MatrixXd::Zero(colMatrix.rows(),1);
      V = Eigen::MatrixXd::Zero(rowMatrix.rows(),1);
      rank = 1;
    }
  }else{
    std::cout<<"Error! Unknown operation mode"<<std::endl;
    exit(EXIT_FAILURE);
  }
  return rank;
}

template void SVD_LowRankApprox<Eigen::MatrixXd>(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank, const int minRank = -1);
template void SVD_LowRankApprox<kernelMatrix>(const kernelMatrix & matrixData,Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank, const int minRank = -1);
template double partialPivACA_LowRankApprox<Eigen::MatrixXd>(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int maxRank = -1,const int minPivot = 0);
template double partialPivACA_LowRankApprox<kernelMatrix>(const kernelMatrix & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int maxRank = -1,const int minPivot = 0);
template double fullPivACA_LowRankApprox<Eigen::MatrixXd>(const Eigen::MatrixXd & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int minPivot = 0);
template double fullPivACA_LowRankApprox<kernelMatrix>(const kernelMatrix & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank = -1,const int minPivot = 0);
template void PS_Boundary_LowRankApprox<Eigen::MatrixXd>(const Eigen::MatrixXd & matrixData,const Eigen::SparseMatrix<double> graphData,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j, const int numRows, const int numCols,const double tolerance,int & calculatedRank, const int maxDepth = 2,const std::string savePath = "none");  
template void PS_Boundary_LowRankApprox<kernelMatrix>(const kernelMatrix & matrixData,const Eigen::SparseMatrix<double> graphData,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j, const int numRows, const int numCols,const double tolerance,int & calculatedRank, const int maxDepth = 2,const std::string savePath = "none");  
