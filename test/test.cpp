#include "HODLR_Matrix.hpp"
#include "user_IndexTree.hpp"
#include "helperFunctions.hpp"
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestCase.h>
#include <cppunit/TestCaller.h>
#include <cppunit/extensions/HelperMacros.h>

/*This file contains various tests for the HODLR_Solver package*/

class HODLR_Solver_Test: public CppUnit::TestCase
{
  /*----------------Creating a Test Suite----------------------*/
  CPPUNIT_TEST_SUITE(HODLR_Solver_Test);
  
  CPPUNIT_TEST(recLU_Solver_Test);
  CPPUNIT_TEST(extendedSp_Solver_Test);
  CPPUNIT_TEST(recLU_Solver_Schur_Test_9k);
  CPPUNIT_TEST(recLU_Solver_Schur_Test_12k);
  CPPUNIT_TEST(recLU_Solver_Schur_Test_16k);
  CPPUNIT_TEST(extendedSp_Solver_Simple_Unbalanced_Test);
  CPPUNIT_TEST(extendedSp_Solver_Schur_Unbalanced_Test);
  CPPUNIT_TEST(iterative_Solve_Test);
  CPPUNIT_TEST(assignment_Test_Simple);
  CPPUNIT_TEST(assignment_Test_ExtendedSp);
  CPPUNIT_TEST(blockExtraction_Test);
  CPPUNIT_TEST(extendAdd_LowRankToHODLR_Test);
  CPPUNIT_TEST(extendAdd_DenseToHODLR_Test);
  CPPUNIT_TEST(extend_HODLRtoHODLR_Test);
  CPPUNIT_TEST(extendAdd_HODLRtoHODLR_Test);
  

  CPPUNIT_TEST(boundaryFinder_Test);
  CPPUNIT_TEST(boundaryFinder_lowRank_Test);


  CPPUNIT_TEST_SUITE_END();

public:
  HODLR_Solver_Test(): CppUnit::TestCase("HODLR Solver Test"){}

  
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
      CPPUNIT_ASSERT(abs((testErr - refErr)/refErr) < 10);
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
      CPPUNIT_ASSERT(abs((testErr - refErr)/refErr) < 10);
    }
    testFile.close();
    refrenceFile.close();
  }



 /* Function : recLU_Solver_Schur_Test_9k
   * ------------------------------
   * This function tests the implicit recursive LU solver on a Schur complement matrix. The test is done on a 9kx9k Schur complement of a turbine blade geometry with two cooling paths. 
   * The functions checks if the solver solves the mentioned matrix with the expected accuracy.
   */ 
  void recLU_Solver_Schur_Test_9k(){
    std::cout<<"Testing recursive implicit LU solver on a 9k Schur complement matrix...."<<std::endl;
 
    // Create custom indexing tree
    user_IndexTree usrTree;
    usrTree.rootNode = new user_IndexTree::node;
    usrTree.rootNode->splitIndex = 1498;
    usrTree.rootNode->topOffDiag_minRank = -1;                                                               
    usrTree.rootNode->bottOffDiag_minRank = -1;
    usrTree.rootNode->LR_Method = "partialPiv_ACA";
    
    user_IndexTree::node* leftNode = new user_IndexTree::node;
    user_IndexTree::node* rightNode = new user_IndexTree::node;
    usrTree.rootNode->left = leftNode;
    usrTree.rootNode->right = rightNode;
    
    rightNode->splitIndex = 6497;
    rightNode->topOffDiag_minRank = -1;                                                               
    rightNode->bottOffDiag_minRank = -1;
    rightNode->LR_Method = "partialPiv_ACA";
    
    leftNode->splitIndex = 749;
    leftNode->topOffDiag_minRank = -1;                                                               
    leftNode->bottOffDiag_minRank = -1;
    leftNode->LR_Method = "PS_Cheby";
     
    usrTree.setChildren_NULL(leftNode);
   
    user_IndexTree::node* leftRightNode = new user_IndexTree::node;
    user_IndexTree::node* rightRightNode = new user_IndexTree::node;
    rightNode->left = leftRightNode;
    rightNode->right = rightRightNode;
    
    leftRightNode->splitIndex = 3998;
    leftRightNode->topOffDiag_minRank = -1;                                                               
    leftRightNode->bottOffDiag_minRank = -1;
    leftRightNode->LR_Method = "PS_Cheby";
    
    rightRightNode->splitIndex = 7868;
    rightRightNode->topOffDiag_minRank = -1;                                                               
    rightRightNode->bottOffDiag_minRank = -1;
    rightRightNode->LR_Method = "PS_Cheby";
    
    usrTree.setChildren_NULL(leftRightNode);
    usrTree.setChildren_NULL(rightRightNode);

    // Read input file 
    std::cout<<"         Reading input file....."<<std::endl;
    Eigen::MatrixXd schurComplement = readBinaryIntoMatrixXd("data/blade5000Schur.bin");
    int matrixSize = schurComplement.rows();
    std::cout<<"         Read successfull."<<std::endl;
    // Initialize Solver
    HODLR_Matrix schur_HODLR(schurComplement,30,usrTree);
    schur_HODLR.printResultInfo = true;
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = schurComplement * exactSoln;
    schur_HODLR.set_LRTolerance(1e-6);
    Eigen::VectorXd solverSoln = schur_HODLR.recLU_Solve(inputF);
    Eigen::VectorXd difference = solverSoln - exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    //std::cout<< relError<<std::endl;
    CPPUNIT_ASSERT(relError < 1e-5);    
  }

  /* Function : recLU_Solver_Schur_Test_12k
   * ------------------------------
   * This function tests the implicit recursive LU solver on a Schur complement matrix. The test is done on a 12kx12k Schur complement of a turbine blade geometry with two cooling paths. 
   * The functions checks if the solver solves the mentioned matrix with the expected accuracy.
   */ 
  void recLU_Solver_Schur_Test_12k(){
    std::cout<<"Testing recursive implicit LU solver on a 12k Schur complement matrix...."<<std::endl;
 
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

    // Read input file 
    std::cout<<"         Reading input file....."<<std::endl;
    Eigen::MatrixXd schurComplement = readBinaryIntoMatrixXd("data/blade7000Schur.bin");
    int matrixSize = schurComplement.rows();
    std::cout<<"         Read successfull."<<std::endl;
    
    // Initialize Solver
    HODLR_Matrix schur_HODLR(schurComplement,30,usrTree);
    schur_HODLR.printResultInfo = true;
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = schurComplement*exactSoln;
    schur_HODLR.set_LRTolerance(1e-7);
    Eigen::VectorXd solverSoln = schur_HODLR.recLU_Solve(inputF);
    Eigen::VectorXd difference = solverSoln - exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    //std::cout<< relError<<std::endl;
    CPPUNIT_ASSERT(relError < 1e-6);    
  }

 /* Function : recLU_Solver_Schur_Test_16k
   * ------------------------------
   * This function tests the implicit recursive LU solver on a Schur complement matrix. The test is done on a 16kx16k Schur complement of a turbine blade geometry with two cooling paths. 
   * The functions checks if the solver solves the mentioned matrix with the expected accuracy.
   */ 
  void recLU_Solver_Schur_Test_16k(){
    std::cout<<"Testing recursive implicit LU solver on a 16k Schur complement matrix...."<<std::endl;
 
    // Create custom indexing tree
    user_IndexTree usrTree;
    usrTree.rootNode = new user_IndexTree::node;
    usrTree.rootNode->splitIndex = 8998;
    usrTree.rootNode->topOffDiag_minRank = -1;                                                               
    usrTree.rootNode->bottOffDiag_minRank = -1;
    usrTree.rootNode->LR_Method = "partialPiv_ACA";
    
    user_IndexTree::node* leftNode = new user_IndexTree::node;
    user_IndexTree::node* rightNode = new user_IndexTree::node;
    usrTree.rootNode->left = leftNode;
    usrTree.rootNode->right = rightNode;
    
    rightNode->splitIndex = 13932;
    rightNode->topOffDiag_minRank = -1;                                                               
    rightNode->bottOffDiag_minRank = -1;
    rightNode->LR_Method = "partialPiv_ACA";
    
    leftNode->splitIndex = 4499;
    leftNode->topOffDiag_minRank = -1;                                                               
    leftNode->bottOffDiag_minRank = -1;
    leftNode->LR_Method = "PS_Cheby";
     
    usrTree.setChildren_NULL(leftNode);
   
    user_IndexTree::node* leftRightNode = new user_IndexTree::node;
    user_IndexTree::node* rightRightNode = new user_IndexTree::node;
    rightNode->left = leftRightNode;
    rightNode->right = rightRightNode;
    
    leftRightNode->splitIndex = 11465;
    leftRightNode->topOffDiag_minRank = -1;                                                               
    leftRightNode->bottOffDiag_minRank = -1;
    leftRightNode->LR_Method = "PS_Cheby";
    
    rightRightNode->splitIndex = 15282;
    rightRightNode->topOffDiag_minRank = -1;                                                               
    rightRightNode->bottOffDiag_minRank = -1;
    rightRightNode->LR_Method = "PS_Cheby";
    
    usrTree.setChildren_NULL(leftRightNode);
    usrTree.setChildren_NULL(rightRightNode);

    // Read input file 
    std::cout<<"         Reading input file....."<<std::endl;
    Eigen::MatrixXd schurComplement = readBinaryIntoMatrixXd("data/blade9000Schur.bin");
    int matrixSize = schurComplement.rows();
    std::cout<<"         Read successfull."<<std::endl;
    // Initialize Solver
    HODLR_Matrix schur_HODLR(schurComplement,30,usrTree);
    schur_HODLR.printResultInfo = true;
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = schurComplement * exactSoln;
    schur_HODLR.set_LRTolerance(0.3 * 1e-7);
    Eigen::VectorXd solverSoln = schur_HODLR.recLU_Solve(inputF);

    Eigen::VectorXd difference = solverSoln - exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    //std::cout<< relError<<std::endl;
    CPPUNIT_ASSERT(relError < 1e-7);    
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
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,
-2,2);
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
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,
-2,2);
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

  void extendAdd_LowRankToHODLR_Test(){
    std::cout<<"Testing low-rank to HODLR extend-add..."<<std::endl;
    int matrixSize = 10000;
    /*********************************Exact Compression**************************/
    std::cout<<"         Testing Exact Compression...."<<std::endl;
    HODLR_Matrix sampleMatrix;
    Eigen::MatrixXd exact_Matrix  = sampleMatrix.createExactHODLR(10,matrixSize,30);
    std::vector<int> idxVec;
    for (int i = 0; i < matrixSize; i++)
      idxVec.push_back(i);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(idxVec.begin(),idxVec.end(),std::default_random_engine(seed));
    int updateSize = matrixSize/2;
    std::vector<int> extendIdxVec = std::vector<int>(idxVec.begin(),idxVec.begin() + updateSize);
    std::sort(extendIdxVec.begin(),extendIdxVec.end());
    int rank = 10;
    Eigen::MatrixXd U = Eigen::MatrixXd::Random(updateSize,rank);
    Eigen::MatrixXd V = Eigen::MatrixXd::Random(updateSize,rank);
    sampleMatrix.extendAddUpdate(U,V,extendIdxVec,1e-16,"Exact");
    Eigen::MatrixXd HODLR_Result = sampleMatrix.block(0,0,matrixSize,matrixSize);
    Eigen::MatrixXd exact_Update = U * V.transpose();
    Eigen::MatrixXd exact_Result = exact_Matrix + extend(extendIdxVec,matrixSize,exact_Update,0,0,updateSize,updateSize,"RowsCols");
    double error =(HODLR_Result - exact_Result).norm();
    // std::cout<<error<<std::endl;
    CPPUNIT_ASSERT(error < 1e-11);
    /*********************************LU Compression**************************/
    std::cout<<"         Testing LU Compression...."<<std::endl;
    HODLR_Matrix sampleMatrix2;
    exact_Matrix  = sampleMatrix2.createExactHODLR(10,matrixSize,30);
    sampleMatrix2.extendAddUpdate(U,V,extendIdxVec,1e-6,"Compress_LU");
    HODLR_Result = sampleMatrix2.block(0,0,matrixSize,matrixSize);
    exact_Result = exact_Matrix + extend(extendIdxVec,matrixSize,exact_Update,0,0,updateSize,updateSize,"RowsCols");
    error =(HODLR_Result - exact_Result).norm();
    //std::cout<<error<<std::endl;
    CPPUNIT_ASSERT(error < 1e-5);
    /**********************************QR Compression**************************/
    std::cout<<"         Testing QR Compression...."<<std::endl;
    HODLR_Matrix sampleMatrix3;
    exact_Matrix  = sampleMatrix3.createExactHODLR(10,matrixSize,30);
    sampleMatrix3.extendAddUpdate(U,V,extendIdxVec,1e-6,"Compress_QR");
    HODLR_Result = sampleMatrix3.block(0,0,matrixSize,matrixSize);
    exact_Result = exact_Matrix + extend(extendIdxVec,matrixSize,exact_Update,0,0,updateSize,updateSize,"RowsCols");
    error =(HODLR_Result - exact_Result).norm();
    //std::cout<<error<<std::endl;
    CPPUNIT_ASSERT(error < 1e-5);
  }

  void extendAdd_DenseToHODLR_Test(){
    std::cout<<"Testing Dense to HODLR extend-add..."<<std::endl;
    int matrixSize = 1000;
    HODLR_Matrix sampleMatrix;
    Eigen::MatrixXd exact_Matrix  = sampleMatrix.createExactHODLR(10,matrixSize,30);
    std::vector<int> idxVec;
    for (int i = 0; i < matrixSize; i++)
      idxVec.push_back(i);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(idxVec.begin(),idxVec.end(),std::default_random_engine(seed));
    int updateSize = matrixSize/2;
    std::vector<int> extendIdxVec = std::vector<int>(idxVec.begin(),idxVec.begin() + updateSize);
    std::sort(extendIdxVec.begin(),extendIdxVec.end());
    Eigen::MatrixXd D = Eigen::MatrixXd::Random(updateSize,updateSize);
    sampleMatrix.extendAddUpdate(D,extendIdxVec,1e-6,"Compress_LU");
    Eigen::MatrixXd HODLR_Result = sampleMatrix.block(0,0,matrixSize,matrixSize);
    Eigen::MatrixXd exact_Result = exact_Matrix + extend(extendIdxVec,matrixSize,D,0,0,updateSize,updateSize,"RowsCols");
    double error =(HODLR_Result - exact_Result).norm();
    //std::cout<<error<<std::endl;
    CPPUNIT_ASSERT(error < 1e-5);
  }
  
  void extend_HODLRtoHODLR_Test(){
    std::cout<<"Testing HODLR to HODLR extend..."<<std::endl;
    int matrixSize = 10000;
    std::vector<int> idxVec;
    for (int i = 0; i < matrixSize; i++)
      idxVec.push_back(i);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(idxVec.begin(),idxVec.end(),std::default_random_engine(seed));
    int updateSize = matrixSize/2;
    std::vector<int> extendIdxVec = std::vector<int>(idxVec.begin(),idxVec.begin() + updateSize);
     std::sort(extendIdxVec.begin(),extendIdxVec.end());
    HODLR_Matrix D_HODLR;
    Eigen::MatrixXd D = D_HODLR.createExactHODLR(10,updateSize,30);
    D_HODLR.extend(extendIdxVec,matrixSize);
    Eigen::MatrixXd HODLR_Result = D_HODLR.block(0,0,matrixSize,matrixSize);
    Eigen::MatrixXd exact_Result = extend(extendIdxVec,matrixSize,D,0,0,updateSize,updateSize,"RowsCols");
    double error =(HODLR_Result - exact_Result).norm();
    //std::cout<<error<<std::endl;
    CPPUNIT_ASSERT(error < 1e-16);
  } 

  void extendAdd_HODLRtoHODLR_Test(){
    std::cout<<"Testing HODLR to HODLR extend-add..."<<std::endl;
    int matrixSize = 10000;
    HODLR_Matrix sampleMatrix;
    Eigen::MatrixXd exact_Matrix  = sampleMatrix.createExactHODLR(10,matrixSize,30);
    std::vector<int> idxVec;
    for (int i = 0; i < matrixSize; i++)
      idxVec.push_back(i);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(idxVec.begin(),idxVec.end(),std::default_random_engine(seed));
    int updateSize = matrixSize/2;
    std::vector<int> extendIdxVec = std::vector<int>(idxVec.begin(),idxVec.begin() + updateSize);
    std::sort(extendIdxVec.begin(),extendIdxVec.end());
    HODLR_Matrix D_HODLR;
    Eigen::MatrixXd D = D_HODLR.createExactHODLR(10,updateSize,30);
    sampleMatrix.extendAddUpdate(D_HODLR,extendIdxVec,1e-6,"Compress_LU");
    Eigen::MatrixXd HODLR_Result = sampleMatrix.block(0,0,matrixSize,matrixSize);
    Eigen::MatrixXd exact_Result = exact_Matrix + extend(extendIdxVec,matrixSize,D,0,0,updateSize,updateSize,"RowsCols");
    double error =(HODLR_Result - exact_Result).norm();
    //std::cout<<error<<std::endl;
    CPPUNIT_ASSERT(error < 1e-4);
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

  void boundaryFinder_lowRank_Test(){
    //Eigen::MatrixXd schurComplement = readBinaryIntoMatrixXd("data/SCHUR_Stiffness/unStructured/300k/300_front_num_2_level_2");
    //Eigen::SparseMatrix<double> graph = readMtxIntoSparseMatrix("data/SCHUR_Stiffness/unStructured/300k/300_front_num_2_level_2_Graph");
    Eigen::MatrixXd schurComplement = readBinaryIntoMatrixXd("data/SCHUR_FETI/Structured/400k/400_front_num_5_level_0");
    Eigen::SparseMatrix<double> graph = readMtxIntoSparseMatrix("data/SCHUR_FETI/Structured/400k/400_front_num_5_level_0_Graph");
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
    */
    
    Eigen::MatrixXd W,K,V;
    int rank;
    int split = schurComplement.rows()/2;
    PS_Boundary_LowRankApprox(schurComplement,graph,W,V,K,split+1,0,schurComplement.rows() - split - 1,split + 1,1e-5,rank,15);
   
    std::cout<<rank<<std::endl;
    Eigen::MatrixXd LR_Matrix = schurComplement.block(split+1,0,schurComplement.rows() - split - 1,split+1);
    
    double absError = (LR_Matrix - W*K*V.transpose()).norm();
    double relError = absError/LR_Matrix.norm();
    std::cout<<relError<<" "<<absError<<std::endl;
    
    HODLR_Matrix schur_HODLR(schurComplement,graph,30);
    //HODLR_Matrix schur_HODLR(schurComplement,30);
    
    schur_HODLR.printResultInfo = true;
    //schur_HODLR.printLevelAccuracy = true;
    
    int matrixSize = schurComplement.rows();
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    Eigen::VectorXd inputF = schurComplement * exactSoln;
    //Eigen::VectorXd solverSoln = schur_HODLR.recLU_Solve(inputF);
    Eigen::VectorXd solverSoln = schur_HODLR.iterative_Solve(inputF,10000,1e-10,1e-1,"PS_Boundary","recLU");
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
     

};   


int main(int argc, char* argv[]){
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(HODLR_Solver_Test::suite());
  runner.run();
  return 0;
}
