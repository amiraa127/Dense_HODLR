

template <typename T>
double fullPivACA_LowRankApprox(const T & matrixData,Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int min_j, const int numRows, const int numCols, const double tolerance, int & calculatedRank,const int minRank,const int maxRank,const int minPivot){

  if (maxRank == 0){
    W = Eigen::MatrixXd::Zero(numRows,1);
    V = Eigen::MatrixXd::Zero(numCols,1);
    calculatedRank = 1;
    return 1;
  }
 
  int rankUpperBound = std::min(numRows,numCols);
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
  while (((epsilon > tolerance) || (k < minRank)) && (k < rankUpperBound)){
   
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
    //std::cout<<absMaxValue<<" "<<minPivot<<std::endl;
    if (absMaxValue <= minPivot){
      break;
    }
    double currPivot = 1./maxValue;
    // Write to W & V

    W.col(k) = currPivot * residualMatrix.col(currColIdx);
    V.col(k) = residualMatrix.row(currRowIdx);

    
    // Update residual matrix
    residualMatrix -= W.col(k) * V.col(k).transpose();
 
    
    // Calculate epsilon
    epsilon = residualMatrix.norm()/origMatrixNorm;
    
    // Set Values for next iteration
    k++;
    if (k == maxRank)
      break;
    
  }
  calculatedRank = k;
  // Return zero for zero matrix
  if (k == 0){
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


template <typename T>
void extractRowsCols(const T & matrixData, int min_i,int min_j,int numRows,int numCols,Eigen::MatrixXd & W,Eigen::MatrixXd & V,const std::vector<int> & rowIndex,const std::vector<int> & colIndex,const double tolerance,int & calculatedRank,const std::string mode,const int maxRank){
  
  
  int numRowsSelect = rowIndex.size();
  int numColsSelect = colIndex.size();

  int numPoints = std::max(numRowsSelect,numColsSelect);
  Eigen::MatrixXd tempW = Eigen::MatrixXd::Zero(numRows,numPoints);
  Eigen::MatrixXd tempV = Eigen::MatrixXd::Zero(numCols,numPoints);
  Eigen::MatrixXd tempK = Eigen::MatrixXd::Zero(numPoints,numPoints);
  
  //fill V and W
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

  if (mode == "fullPivACA" || maxRank > 0){
    Eigen::MatrixXd K_W,K_V;
    int kRank;
    fullPivACA_LowRankApprox(tempK,K_W,K_V,0,0,tempK.rows(),tempK.cols(),tolerance,kRank,-1,maxRank);
    W = ((K_V.transpose() * K_V).transpose().partialPivLu().solve(K_V.transpose()) * tempW.transpose()).transpose();
    V =((K_W.transpose() * K_W).partialPivLu().solve(K_W.transpose()) * tempV.transpose()).transpose();
    calculatedRank = kRank;
  }
  else if (mode == "fullPivLU") {
    /*
    int numPoints = std::max(numRowsSelect,numColsSelect);
    Eigen::MatrixXd tempW = Eigen::MatrixXd::Zero(numRows,numPoints);
    Eigen::MatrixXd tempV = Eigen::MatrixXd::Zero(numCols,numPoints);
    Eigen::MatrixXd tempK = Eigen::MatrixXd::Zero(numPoints,numPoints);
    
    //fill V and W
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
    */
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
    /*
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
    
    */
    Eigen::MatrixXd svdW,svdV,svdK;
    //calculatedRank = SVD_LowRankApprox_(KTemp,tolerance, &svdW, &svdV, &svdK);
    calculatedRank = SVD_LowRankApprox_(tempK,tolerance, &svdW, &svdV, &svdK);
   
    
    //calculate K
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(calculatedRank,calculatedRank);
    for (int i = 0; i < calculatedRank; i++)
      K(i,i) = 1./svdK(i,i);
    
    //calculate W and V
    //W = WTemp * svdV * K;
    //V = VTemp * svdW;
    W = tempW * svdV * K;
    V = tempV * svdW;
   
  }else{
    std::cout<<"Error! Unknown operation mode."<<std::endl;
    exit(EXIT_FAILURE);
  }
}


template <typename T>
void PS_Boundary_LowRankApprox(const T & matrixData,const Eigen::SparseMatrix<double> graphData,Eigen::MatrixXd & W, Eigen::MatrixXd & V,const int min_i, const int min_j, const int numRows, const int numCols,double tolerance,int & calculatedRank, const int depth, const std::string savePath,const int numSel,const int maxRank){

  if (maxRank == 0){
    W = Eigen::MatrixXd::Zero(numRows,1);
    V = Eigen::MatrixXd::Zero(numCols,1);
    calculatedRank = 1;
    return;
  }

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
  
  //int noInteraction = identifyBoundary(graphData.block(minIdx,minIdx,numPoints,numPoints),rowSet,colSet,rowPos,colPos,depth);
  //int noInteraction = identifyBoundary(graphData,minIdx,minIdx,numPoints,numPoints,rowSet,colSet,rowPos,colPos,depth,numSel);
  int noInteraction = identifyBoundary(graphData,min_i,min_j,numRows,numCols,rowPos,colPos,depth,numSel);

  
  if (noInteraction == 1){
    W = Eigen::MatrixXd::Zero(numRows,1);
    V = Eigen::MatrixXd::Zero(numCols,1);
    calculatedRank = 1;
    return;
  }
  
  createIdxFromBoundaryMap(rowPos,colPos,depth,rowIdx,colIdx);
  
  /*
  // Adjust for offsets
  for (unsigned int i = 0; i < rowIdx.size();i++)
    rowIdx[i] -= offset_i;
  for (unsigned int i = 0; i < colIdx.size();i++)
    colIdx[i] -= offset_j;
  */
  
  extractRowsCols(matrixData,min_i,min_j,numRows,numCols,W,V,rowIdx,colIdx,tolerance,calculatedRank,"fullPivLU",maxRank);
  //std::cout<<calculatedRank<<std::endl;
#if 0
    for(std::map<int,std::vector<int> >::iterator iter = rowPos.begin(); iter != rowPos.end(); ++iter){
    std::cout<<iter->first<<":";
    for (unsigned int i = 0; i < iter->second.size();i++)
    std::cout<<iter->second[i]<<" ";
    std::cout<<std::endl;
    }
#endif
 
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
     
     for (unsigned int i = 0; i < colIdx.size(); i++)
       colIdxVec[i] = colIdx[i] + offset_j;
     for (unsigned int i = 0; i < rowIdx.size(); i++)
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
	   
   