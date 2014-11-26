#include "HODLR_Matrix.hpp"
#include "user_IndexTree.hpp"
#include "helperFunctions.hpp"
#include "matrixIO.hpp"
#include "kernel.hpp"
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/extensions/HelperMacros.h>

/*This file contains various tests for the HODLR_Solver package*/
double quadraticKernel(int i, int j, void* pointsCoordPtr){
  if (i == j)
    return 0;
  double r = (*(Eigen::VectorXd*)pointsCoordPtr)(i) - (*(Eigen::VectorXd*)pointsCoordPtr)(j);
  return r * r + 1;
}

class HODLR_Matrix_Test: public CppUnit::TestCase
{
  /*----------------Creating a Test Suite----------------------*/
  CPPUNIT_TEST_SUITE(HODLR_Matrix_Test);
  
  CPPUNIT_TEST(recSM_Solver_Test_Random);

  CPPUNIT_TEST(recLU_Solver_Test);
  //CPPUNIT_TEST(extendedSp_Solver_Test);
  //CPPUNIT_TEST(extendedSp_Solver_Simple_Unbalanced_Test);
  //CPPUNIT_TEST(extendedSp_Solver_Schur_Unbalanced_Test);
  CPPUNIT_TEST(iterative_Solve_Test);
  CPPUNIT_TEST(assignment_Test_Simple);
  //CPPUNIT_TEST(assignment_Test_ExtendedSp);
  CPPUNIT_TEST(blockExtraction_Test);
  CPPUNIT_TEST(splitAtTop_Test);

  //CPPUNIT_TEST(boundaryFinder_Test);
  //CPPUNIT_TEST(boundaryFinder_lowRank_Test);
  
  CPPUNIT_TEST(kernelSolver_Test);
  CPPUNIT_TEST(determinant_Test);
  CPPUNIT_TEST_SUITE_END();

public:
  HODLR_Matrix_Test(): CppUnit::TestCase("HODLR Matrix Test"){}

  /* Function : recSM_Solver_Test_Radom
   * --------------------------------------
   * This function tests the recursive Sherman Morrison Solver ona 10kx10k random HODLR matrix and checks the accuracy.
   */
  void recSM_Solver_Test_Random(){
    std::cout<<"Testing recursive Sherman Morrison solver on a random matrix...."<<std::endl;
    int matrixSize = 10000;
    HODLR_Matrix sampleHODLR;
    Eigen::MatrixXd sampleMatrix = sampleHODLR.createExactHODLR(10,matrixSize,50);
    Eigen::VectorXd exactSoln  = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd sampleRHS  = sampleMatrix * exactSoln;
    Eigen::MatrixXd solverSoln = sampleHODLR.recSM_Solve(sampleRHS);
    double error = (solverSoln - exactSoln).norm()/exactSoln.norm();
    CPPUNIT_ASSERT(error < 1e-8);
  }

  /* Function : recLU_Solver_Test
   * ------------------------------
   * This function tests the recursive LU solver on a 10kx10k dense interaction matrix with an inverse multiquadratic kernel.
   * The functions checks if the solver solves the mentioned matrix with the expected accuracy for a given ACA tolerance.
   */
  void recLU_Solver_Test(){
    int intervalSize = 10000;
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,intervalSize,-2,2);
    std::cout<<"Testing recursive LU solver on a radial basis function matrix...."<<std::endl;
    std::cout<<"         Benchmarking inverse multi quadratic kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"data/recLU/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"recLU");
    std::ifstream testFile,refrenceFile;
    testFile.open("data/recLU/inverseMultiQuadratic_0");
    refrenceFile.open("data/recLU/inverseMultiQuadratic_ref_0");
    std::string testLine, refLine;
    double testTol, refTol, testErr, refErr;
    while ((!testFile.eof()) && (!refrenceFile.eof())){
      getline(testFile,testLine);
      getline(refrenceFile,refLine);
      sscanf(testLine.c_str(),"%lf %lf",&testTol,&testErr);
      sscanf(refLine.c_str(),"%lf %lf",&refTol,&refErr);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(refTol,testTol,1e-10);
      CPPUNIT_ASSERT(fabs((testErr - refErr)/refErr) < 10);
    }
    testFile.close();
    refrenceFile.close();
  }
  
  /* Function : extendedSp_Solver_Test
   * ------------------------------
   * This function tests the extended sparsification solver on a 10kx10k dense interaction matrix with an inverse multiquadratic kernel.
   * The functions checks if the solver solves the mentioned matrix with the expected accuracy for a given ACA tolerance.
   */
  void extendedSp_Solver_Test(){
    int intervalSize = 10000;
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,intervalSize,-2,2);
    std::cout<<"Testing extended sparsification solver on a radial basis function matrix...."<<std::endl;
    std::cout<<"         Benchmarking inverse multi quadratic kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"data/extendedSp/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"extendedSp");
    std::ifstream testFile,refrenceFile;
    testFile.open("data/extendedSp/inverseMultiQuadratic_0");
    refrenceFile.open("data/extendedSp/inverseMultiQuadratic_ref_0");
    std::string testLine, refLine;
    double testTol, refTol, testErr, refErr;
    while ((!testFile.eof()) && (!refrenceFile.eof())){
      getline(testFile,testLine);
      getline(refrenceFile,refLine);
      sscanf(testLine.c_str(),"%lf %lf",&testTol,&testErr);
      sscanf(refLine.c_str(),"%lf %lf",&refTol,&refErr);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(refTol,testTol,1e-10);
      CPPUNIT_ASSERT(fabs((testErr - refErr)/refErr) < 10);
    }
    testFile.close();
    refrenceFile.close();
  }

  /* Function : extendedSp_Solver_Simple_Unbalanced
   *------------------------------------------
   * This function tests the extended sparsification LU solver on a simple unbalanced tree. 
   * The HODLR tree is unbalanced on the left and right.
   * The function checks if the solver solves the mentioned matrix with the exoected accuracy.
   */
  void extendedSp_Solver_Simple_Unbalanced_Test(){
    std::cout<<"Testing extended sparsification solver on a simple unbalanced tree...."<<std::endl;
    int matrixSize = 10000;
    
    user_IndexTree usrTree;
    usrTree.rootNode = new user_IndexTree::node;
    usrTree.rootNode->splitIndex = 1500; 
    usrTree.rootNode->topOffDiag_minRank = -1;                                                               
    usrTree.rootNode->bottOffDiag_minRank = -1;
    usrTree.rootNode->LR_Method = "partialPiv_ACA";
    
    user_IndexTree::node* leftNode = new user_IndexTree::node;
    user_IndexTree::node* rightNode = new user_IndexTree::node;
    usrTree.rootNode->left = leftNode;
    usrTree.rootNode->right = rightNode;
    
    rightNode->splitIndex = 7500;
    rightNode->topOffDiag_minRank = -1;                                                               
    rightNode->bottOffDiag_minRank = -1;
    rightNode->LR_Method = "partialPiv_ACA";
    
    leftNode->splitIndex = 1000;
    leftNode->topOffDiag_minRank = -1;                                                               
    leftNode->bottOffDiag_minRank = -1;
    leftNode->LR_Method = "partialPiv_ACA";
    
    usrTree.setChildren_NULL(leftNode);
    usrTree.setChildren_NULL(rightNode);
    
    Eigen::MatrixXd sampleMatrix = makeMatrix1DUniformPts(-1,1,-1,1,matrixSize,matrixSize,0,inverseMultiQuadraticKernel);
    HODLR_Matrix sample_HODLR(sampleMatrix,30,usrTree);
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = sampleMatrix * exactSoln;
    sample_HODLR.set_LRTolerance(1e-7);  
    Eigen::VectorXd solverSoln = sample_HODLR.extendedSp_Solve(inputF);
    Eigen::VectorXd difference = solverSoln-exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    CPPUNIT_ASSERT(relError < 1e-6); 
  }

  
  /* Function : extendedSp_Solver_Schur_Unbalanced
   *------------------------------------------
   * This function tests the extended sparsification LU solver on an unbalanced tree. 
   * The HODLR tree is unbalanced on the left and right and has a structure similar to a Schur complement HODLR tree.
   * The function checks if the solver solves the mentioned matrix with the exoected accuracy.
   */   
  void extendedSp_Solver_Schur_Unbalanced_Test(){
    std::cout<<"Testing extended sparsification solver on a Schur complement unbalanced tree...."<<std::endl;
    // Create custom indexing tree
    user_IndexTree usrTree;
    usrTree.rootNode = new user_IndexTree::node;
    usrTree.rootNode->splitIndex = 2098;
    usrTree.rootNode->topOffDiag_minRank = -1;                                                               
    usrTree.rootNode->bottOffDiag_minRank = -1;
    usrTree.rootNode->LR_Method = "partialPiv_ACA";
    

    user_IndexTree::node* leftNode = new user_IndexTree::node;
    user_IndexTree::node* rightNode = new user_IndexTree::node;
    usrTree.rootNode->left = leftNode;
    usrTree.rootNode->right = rightNode;
    
    rightNode->splitIndex = 9097;
    rightNode->topOffDiag_minRank = -1;                                                               
    rightNode->bottOffDiag_minRank = -1;
    rightNode->LR_Method = "partialPiv_ACA";
    
    leftNode->splitIndex = 1048;
    leftNode->topOffDiag_minRank = -1;                                                               
    leftNode->bottOffDiag_minRank = -1;
    leftNode->LR_Method = "PS_Cheby";
     
    usrTree.setChildren_NULL(leftNode);
   
    user_IndexTree::node* leftRightNode = new user_IndexTree::node;
    user_IndexTree::node* rightRightNode = new user_IndexTree::node;
    rightNode->left = leftRightNode;
    rightNode->right = rightRightNode;
    
    leftRightNode->splitIndex = 5597;
    leftRightNode->topOffDiag_minRank = -1;                                                               
    leftRightNode->bottOffDiag_minRank = -1;
    leftRightNode->LR_Method = "PS_Cheby";
   

    rightRightNode->splitIndex = 11016;
    rightRightNode->topOffDiag_minRank = -1;                                                               
    rightRightNode->bottOffDiag_minRank = -1;
    rightRightNode->LR_Method = "PS_Cheby";
    
    usrTree.setChildren_NULL(leftRightNode);
    usrTree.setChildren_NULL(rightRightNode);

    // Create input Matrix 
    int matrixSize = 12936;
    Eigen::MatrixXd sampleMatrix = makeMatrix1DUniformPts(-1,1,-1,1,matrixSize,matrixSize,0,inverseMultiQuadraticKernel);

    // Initialize Solver
    HODLR_Matrix sampleHODLR(sampleMatrix,30,usrTree);
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,
-2,2);
    Eigen::VectorXd inputF = sampleMatrix * exactSoln;
    sampleHODLR.set_LRTolerance(1e-8);
    Eigen::VectorXd solverSoln = sampleHODLR.extendedSp_Solve(inputF);
    Eigen::VectorXd difference = solverSoln - exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    CPPUNIT_ASSERT(relError < 1e-6);    
  }


  void iterative_Solve_Test(){
    int matrixSize = 10000;
    std::cout<<"Testing iterative solver on a radial basis function matrix...."<<std::endl;
    Eigen::MatrixXd sampleMatrix = makeMatrix1DUniformPts(-1,1,-1,1,matrixSize,matrixSize,0,inverseMultiQuadraticKernel);
    HODLR_Matrix sample_HODLR(sampleMatrix,30);
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = sampleMatrix * exactSoln;
    //sample_HODLR.set_LRTolerance(1e-7);  
    sample_HODLR.printResultInfo = true;
    Eigen::VectorXd solverSoln = sample_HODLR.iterative_Solve(inputF,10000,1e-10,1e-2,"partialPiv_ACA","extendedSp");
    Eigen::VectorXd difference = solverSoln-exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    //std::cout<<relError<<std::endl;
    CPPUNIT_ASSERT(relError < 1e-6);
  }

  void assignment_Test_Simple(){
    std::cout<<"Testing assignment operator with recLU solver....."<<std::endl;
    int matrixSize = 12936;
    Eigen::MatrixXd sampleMatrix = makeMatrix1DUniformPts(-1,1,-1,1,matrixSize,matrixSize,0,inverseMultiQuadraticKernel);

    // Initialize Solver
    HODLR_Matrix sampleHODLR(sampleMatrix,30);
    sampleHODLR.set_LRTolerance(1e-8);
    HODLR_Matrix copy_sampleHODLR = sampleHODLR;
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = sampleMatrix * exactSoln;
    Eigen::VectorXd solverSoln = copy_sampleHODLR.recLU_Solve(inputF);
    Eigen::VectorXd difference = solverSoln - exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    //std::cout<<relError<<std::endl;
    CPPUNIT_ASSERT(relError < 1e-7);    
  }

  void assignment_Test_ExtendedSp(){
    std::cout<<"Testing assignment operator with extended sparsification solver....."<<std::endl;
    int matrixSize = 12936;
    Eigen::MatrixXd sampleMatrix = makeMatrix1DUniformPts(-1,1,-1,1,matrixSize,matrixSize,0,inverseMultiQuadraticKernel);
    
    // Initialize Solver
    HODLR_Matrix sampleHODLR(sampleMatrix,30);
    sampleHODLR.set_LRTolerance(1e-8);
    HODLR_Matrix copy_sampleHODLR = sampleHODLR;
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = sampleMatrix * exactSoln;
    Eigen::VectorXd solverSoln = copy_sampleHODLR.extendedSp_Solve(inputF);
    Eigen::VectorXd difference = solverSoln - exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    //std::cout<<relError<<std::endl;
    CPPUNIT_ASSERT(relError < 1e-6);    
  }

  void blockExtraction_Test(){
    std::cout<<"Testing HODLR block function....."<<std::endl;
    HODLR_Matrix sampleMatrix;
    int matrixSize = 10000;
    Eigen::MatrixXd exact_Matrix  = sampleMatrix.createExactHODLR(10,matrixSize,30);
    Eigen::MatrixXd extract_Full  = sampleMatrix.block(0,0,matrixSize,matrixSize);
    Eigen::MatrixXd extract_Part  = sampleMatrix.block(matrixSize/4,matrixSize/4,matrixSize/2,matrixSize/2);
    Eigen::MatrixXd extract_Row   = sampleMatrix.row(50);
    Eigen::MatrixXd extract_Col   = sampleMatrix.col(6532);
    CPPUNIT_ASSERT((extract_Full - exact_Matrix).norm() < 1e-16);
    CPPUNIT_ASSERT((extract_Part - exact_Matrix.block(matrixSize/4,matrixSize/4,matrixSize/2,matrixSize/2)).norm() < 1e-16);
    CPPUNIT_ASSERT((extract_Row - exact_Matrix.row(50)).norm() < 1e-16);
    CPPUNIT_ASSERT((extract_Col - exact_Matrix.col(6532)).norm() < 1e-16);
  }
  
  void splitAtTop_Test(){
    int matrixSize = 10000;
    HODLR_Matrix sampleMatrix;
    Eigen::MatrixXd parentMatrix = sampleMatrix.createExactHODLR(10,matrixSize,30);
    HODLR_Matrix topDiag,bottDiag;
    //sampleMatrix.splitAtTop(topDiag,bottDiag);
    splitAtTop(sampleMatrix,topDiag,bottDiag);
    int topDiagSize = topDiag.get_MatrixSize();
    std::cout<<(parentMatrix.topLeftCorner(topDiagSize,topDiagSize) - topDiag.block(0,0,topDiagSize,topDiagSize)).norm()<<std::endl;
  }


  void boundaryFinder_Test(){
    Eigen::MatrixXd inputMatrix = Eigen::MatrixXd::Zero(12,12);
    for (int i = 0; i < 12; i++)
      inputMatrix(i,i)   = 10;
    inputMatrix(0,1)   =  2; inputMatrix(0,4)  =  4;
    inputMatrix(1,2)   = -1; inputMatrix(1,5)  = -3;
    inputMatrix(2,3)   =  2 ; inputMatrix(2,6)  = -1;
    inputMatrix(3,7)   =  5 ;
    inputMatrix(4,5)   = -3; inputMatrix(4,8)  =  2;
    inputMatrix(5,6)   = -2; inputMatrix(5,9)  = -1;
    inputMatrix(6,7)   =  4; inputMatrix(6,10) = -2; inputMatrix(6,11) = 3 ; 
    inputMatrix(7,11)  = -1;
    inputMatrix(8,9)   = -3;
    inputMatrix(9,10)  =  5;
    inputMatrix(10,11) =  2;
    
    for (int i = 0; i < 12;i++)
      for (int j = i; j < 12; j++)
	inputMatrix(j,i) = inputMatrix(i,j);
  
    Eigen::SparseMatrix<double> inputSpMatrix = inputMatrix.sparseView();
    std::map<int,std::vector<int> > rowPos,colPos;
    std::set<int> rowSet,colSet;
    for (int i = 0; i <= 5; i++)
      rowSet.insert(i);
    
    for (int i = 6; i <= 11; i++)
      colSet.insert(i);
    /*
    identifyBoundary(inputSpMatrix,rowSet,colSet,rowPos,colPos);
      
    for(std::map<int,std::vector<int> >::iterator iter = rowPos.begin(); iter != rowPos.end(); ++iter){
      std::cout<<iter->first<<":";
      for (unsigned int i = 0; i < iter->second.size();i++)
	std::cout<<iter->second[i]<<" ";
      std::cout<<std::endl;
    }
    
    for(std::map<int,std::vector<int> >::iterator iter = colPos.begin(); iter != colPos.end(); ++iter){
      std::cout<<iter->first<<":";
      for (unsigned int i = 0; i < iter->second.size();i++)
	std::cout<<iter->second[i]<<" ";
      std::cout<<std::endl;
    }
    */
  }

  /*
  void boundaryFinder_lowRank_Test(){
    //Eigen::MatrixXd schurComplement = readBinaryIntoMatrixXd("data/SCHUR_Stiffness/unStructured/300k/300_front_num_2_level_2");
    //Eigen::SparseMatrix<double> graph = readMtxIntoSparseMatrix("data/SCHUR_Stiffness/unStructured/300k/300_front_num_2_level_2_Graph");
    Eigen::MatrixXd schurComplement = readBinaryIntoMatrixXd("data/SCHUR_FETI/Structured/400k/400_front_num_0_level_2");
    Eigen::SparseMatrix<double> graph = readMtxIntoSparseMatrix("data/SCHUR_FETI/Structured/400k/400_front_num_0_level_2_Graph");
    //Eigen::MatrixXd schurComplement1 = readBinaryIntoMatrixXd("../benchmarks/boundaryLR/stiffness/unstructured/2D/input/blade9000Schur.bin");
    //Eigen::SparseMatrix<double> graph1 = readMtxIntoSparseMatrix("../benchmarks/boundaryLR/stiffness/unstructured/2D/input/blade9000Graph");
    
    //Eigen::MatrixXd schurComplement = schurComplement1.block(0,0,8999,8999);
    //Eigen::SparseMatrix<double> graph = graph1.block(0,0,8999,8999);
    
    /*
    // identify boundary
    std::map<int,std::vector<int> > rowPos,colPos;
    std::set<int> rowSet,colSet;
    int split = schurComplement.rows()/2;
    for (int i = 0; i <= split; i++)
      rowSet.insert(i);
    
    for (int i = split + 1; i < schurComplement.cols() ; i++)
      colSet.insert(i);
    
    identifyBoundary(graph,rowSet,colSet,rowPos,colPos);
    
    
    for (int i = 0;i <= 10; i++){
      std::vector<int> rowIdx,colIdx;
      createIdxFromBoundaryMap(rowPos,colPos,i,rowIdx,colIdx);
      
      // adjust colIndex
      int minCol = split + 1;
      for (unsigned int i = 0; i < colIdx.size();i++)
	colIdx[i] = colIdx[i] - minCol;
      // Extract
      Eigen::MatrixXd W,K,V;
      Eigen::MatrixXd LR_Matrix = schurComplement.block(0,split+1,split+1,schurComplement.cols() - split - 1);
      extractRowsCols(W,K,V,LR_Matrix,rowIdx,colIdx);
      double absError = (LR_Matrix - W*K*V.transpose()).norm();
      double relError = absError/LR_Matrix.norm();
      std::cout<<"depth = "<<i<<" rank = "<<colIdx.size()<<" abs error = "<<absError<< " rel error = "<<relError<<std::endl;
    
    }
    
    
    Eigen::MatrixXd W,V;
    int rank;
    int split = schurComplement.rows()/2;
    
    //PS_Boundary_LowRankApprox(schurComplement,graph,W,V,K,split+1,0,schurComplement.rows() - split - 1,split + 1,1e-5,rank,1500);
    PS_Boundary_LowRankApprox(schurComplement,graph,W,V,split+1,0,schurComplement.rows() - split - 1,split + 1,1e-5,rank,1500);
    
    std::cout<<rank<<std::endl;
    Eigen::MatrixXd LR_Matrix = schurComplement.block(split+1,0,schurComplement.rows() - split - 1,split+1);
    
    //double absError = (LR_Matrix - W*K*V.transpose()).norm();
   double absError = (LR_Matrix - W * V.transpose()).norm();
    
    double relError = absError/LR_Matrix.norm();
    std::cout<<relError<<" "<<absError<<std::endl;
    
    HODLR_Matrix schur_HODLR(schurComplement,graph,100);
    //HODLR_Matrix schur_HODLR(schurComplement,30);
    
    schur_HODLR.printResultInfo = true;
    //schur_HODLR.printLevelAccuracy = true;
    
    int matrixSize = schurComplement.rows();
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = schurComplement * exactSoln;
    //schur_HODLR.printLevelAccuracy = true;
    schur_HODLR.set_BoundaryDepth(1);
    //schur_HODLR.set_LRTolerance(1e-1);
    //Eigen::VectorXd solverSoln = schur_HODLR.recLU_Solve(inputF);
   
    Eigen::VectorXd solverSoln = schur_HODLR.iterative_Solve(inputF,100,1e-10,1e-1,"PS_Boundary","recLU");
    Eigen::VectorXd difference = solverSoln - exactSoln;
    relError = difference.norm()/exactSoln.norm();
    std::cout<<relError<<std::endl;
  
    double startTime = clock();
    Eigen::PartialPivLU<Eigen::MatrixXd> LU (schurComplement);
    solverSoln = LU.solve(inputF);
    double endTime = clock();
    double time = (endTime - startTime)/CLOCKS_PER_SEC;
    std::cout<<time<<std::endl;
    std::cout<<schurComplement.rows()<<std::endl;
    
  }
  */

  void kernelSolver_Test(){
    int numPoints = 10000;
    Eigen::VectorXd X = Eigen::VectorXd::LinSpaced(Eigen::Sequential,numPoints,-2,2);
    HODLR_Matrix kernelHODLR(numPoints,numPoints,quadraticKernel,&X,30);
    kernelMatrix exactMatrixKernel(numPoints,numPoints,quadraticKernel,&X);
    Eigen::MatrixXd exactMatrix = exactMatrixKernel.block(0,0,numPoints,numPoints);
    Eigen::VectorXd inputF = exactMatrix * X;
    Eigen::VectorXd solverSoln = kernelHODLR.recLU_Solve(inputF);
    Eigen::VectorXd difference = solverSoln - X;
    double relError = difference.norm()/X.norm();
    std::cout<<relError<<std::endl;
  }

  void determinant_Test(){
    
    int matrixSize = 3000;
    HODLR_Matrix sample_HODLR;
    Eigen::MatrixXd sampleMatrix = sample_HODLR.createExactHODLR(100,matrixSize,50);
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(sampleMatrix);
    double logAbsDet = 0;
    Eigen::MatrixXd luMatrix = lu.matrixLU();
    for (int i = 0; i < luMatrix.rows(); i++)
      logAbsDet += log(fabs(luMatrix(i,i)));
    double error = fabs(logAbsDet - sample_HODLR.logAbsDeterminant())/fabs(logAbsDet);
    
    std::cout<<error<<std::endl;
    
  }

};   


int main(int argc, char* argv[]){
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(HODLR_Matrix_Test::suite());
  runner.run();
  return 0;
}
