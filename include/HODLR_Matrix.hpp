#ifndef HODLR_MATRIX_HPP
#define HODLR_MATRIX_HPP

#include "HODLR_Tree.hpp"
#include "helperFunctions.hpp"
#include "user_IndexTree.hpp"
#include "recLU_FactorTree.hpp"
#include "lowRank.hpp"
#include "matrixIO.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>

/* Author : Amir Aminfar ---> aminfar@stanford.edu */
class HODLR_Matrix{

public:
  
  bool printLevelRankInfo;
  bool printLevelAccuracy;
  bool printLevelInfo;
  bool printResultInfo;

  /************************************ Constructors *********************************/
  HODLR_Matrix();
  
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix);
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix);
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,Eigen::SparseMatrix<double> &inputGraph);
  
  
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,int inputSizeThreshold);
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,int inputSizeThreshold);
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold);
  
  
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,int inputSizeThreshold,user_IndexTree &input_IndexTree);
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,int inputSizeThreshold,user_IndexTree &input_IndexTree);
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold,user_IndexTree &input_IndexTree);
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold,user_IndexTree &input_IndexTree);
  
  
  HODLR_Matrix(const HODLR_Matrix & rhs); //Copy Constructor

  ~HODLR_Matrix();
  
  /************************************* Create HODLR Structure ***************************************/
  void storeLRinTree();

  /************************************* Solve Methods **********************************/
  Eigen::MatrixXd recLU_Solve(const Eigen::MatrixXd & input_RHS);
  void recLU_Compute();
  Eigen::MatrixXd extendedSp_Solve(const Eigen::MatrixXd & input_RHS);
  Eigen::MatrixXd iterative_Solve(const Eigen::MatrixXd & input_RHS, const int maxIterations, const double stop_tolerance,const double init_LRTolerance,const std::string input_LR_Method, const std::string directSolve_Method);

  /************************************* Attribute Modification **********************************/
  void set_LRTolerance(double tolerance);  
  void set_minPivot(double minPivot);
  void set_def_LRMethod(std::string input_LRMethod);
  void set_FreeMatrixMemory(bool inputVal);
  void set_BoundaryDepth(int inputBoundaryDepth);
  void saveExtendedSp(std::string savePath);
 
  /************************************ Accessing Attributes ************************************/
  double get_recLU_FactorizationTime() const;
  double get_recLU_SolveTime() const;
  double get_recLU_TotalTime() const;
  double get_extendedSp_AssemblyTime() const; 
  double get_extendedSp_FactorizationTime() const;
  double get_extendedSp_SolveTime() const;
  double get_extendedSp_TotalTime() const;
  double get_LR_ComputationTime() const;
  double get_totalIter_SolveTime() const;
  double get_MatrixSize() const;

  int rows() const;
  int cols() const;
  double norm();
  HODLR_Tree::node* get_TreeRootNode();
  /************************************ Acessing HODLR Entries *******************************/
  Eigen::MatrixXd block(int min_i,int min_j,int numRows,int numCols);
  Eigen::MatrixXd row(int row);
  Eigen::MatrixXd col(int col);
 
  
  HODLR_Matrix  topDiag();
  HODLR_Matrix bottDiag();

  void keepTopDiag();
  void keepBottDiag();

  Eigen::MatrixXd& returnTopOffDiagU();
  Eigen::MatrixXd& returnTopOffDiagV();
  Eigen::MatrixXd& returnTopOffDiagK();
  Eigen::MatrixXd& returnBottOffDiagU();
  Eigen::MatrixXd& returnBottOffDiagV();
  Eigen::MatrixXd& returnBottOffDiagK();

  /******************************** Check ******************************************************/
  void check_Structure();
  Eigen::MatrixXd createExactHODLR(const int rank,int input_MatrixSize,const int inpt_SizeThreshold);

  void saveSolverInfo(const std::string outputFileName);
  
  void freeMatrixData();
  void recalculateSize();
private:

  int sizeThreshold;
  int extendedSp_Size;
  int matrixSize;
  int matrixNumRows;
  int matrixNumCols;

  double recLU_FactorizationTime;
  double recLU_SolveTime;
  double recLU_TotalTime;
  double LR_ComputationTime;
  double extendedSp_AssemblyTime;
  double extendedSp_FactorizationTime;
  double extendedSp_SolveTime;
  double extendedSp_TotalTime;
  double totalIter_SolveTime;
    
  std::vector<double> LR_ComputationLevelTimeVec;
  std::vector<double> recLU_FactorLevelTimeVec;
  std::vector<double> recLU_SolveLevelTimeVec;
  std::vector<double> iter_IterTimeVec;
  std::vector<double> iter_IterAccuracyVec;
  std::vector<double> levelRankAverageVec;
  std::vector<double> levelRankAverageVecCnt;
  
  bool LRStoredInTree;
  bool createdRecLUfactorTree;
  bool assembled_ExtendedSp;
  bool saveExtendedSp_Matrix;
  bool freeMatrixMemory;
  bool freeMatrixMemory_Sp;
  bool matrixDataAvail;
  bool matrixDataAvail_Sp;
  bool graphDataAvail;
  bool isSquareMatrix;

  double LR_Tolerance;
  double minPivot;
  int boundaryDepth;

  HODLR_Tree indexTree;
  recLU_FactorTree recLUfactorTree;
  Eigen::MatrixXd matrixData;
  Eigen::SparseMatrix<double> matrixData_Sp;
  Eigen::SparseMatrix<double> graphData;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > extendedSp_Solver;

  std::string extendedSp_SavePath;

  void setDefaultValues();
  void reset_attributes();
  void initializeInfoVecotrs(int numLevels);

  /****************************recLU Solver Functions*******************************/
 
  void storeLRinTree(HODLR_Tree::node* HODLR_Root);
  void recLU_Factorize();
  Eigen::MatrixXd recLU_Factorize(const Eigen::MatrixXd & input_RHS,const HODLR_Tree::node* HODLR_Root, recLU_FactorTree::node* factorRoot);
  Eigen::MatrixXd recLU_Solve(const Eigen::MatrixXd & input_RHS,const HODLR_Tree::node* HODLR_Root, const recLU_FactorTree::node* factorRoot);

  /**************************extendedSp Solver Functions***************************/
  void findNodesAtLevel(HODLR_Tree::node* HODLR_Root, const int level, std::vector<HODLR_Tree::node*> & outputVector);
  void findLeafNodes(HODLR_Tree::node* HODLR_Root, std::vector<HODLR_Tree::node*>& outputVector);
  void insertDenseBlockIntoSpMatrix(std::vector<Eigen::Triplet<double,int> > & Sp_TripletVec, const Eigen::MatrixXd & denseMatrix, const int start_i, const int start_j);
  void insertIdentityIntoSpMatrix(std::vector<Eigen::Triplet<double,int> > & Sp_TripletVec, const int startIndex_i,const  int startIndex_j, const int matrixSize, const int constant);
  int sumRanks(HODLR_Tree::node* HODLR_Root);
  Eigen::SparseMatrix<double> assembleExtendedSPMatrix();
  
  /***************************Iterative Solver Functions****************************/
  Eigen::MatrixXd oneStep_Iterate(const Eigen::MatrixXd &  prevStep_result,const Eigen::MatrixXd & RHS,const Eigen::MatrixXd & initSolveGuess,Eigen::MatrixXd & prevStep_Product,const std::string directSolve_Method);


  /**********************************Accessing Matrix Entries***************************/
  void fill_Block(Eigen::MatrixXd & blkMatrix,HODLR_Tree::node* root,int min_i,int min_j,int max_i,int max_j);
  void fill_BlockWithLRProduct(Eigen::MatrixXd & blkMatrix,int LR_Min_i,int LR_Min_j, int LR_numRows, int LR_numCols,Eigen::MatrixXd & LR_U,Eigen::MatrixXd & LR_V,int blk_Min_i,int blk_Min_j);

  /******************************Memory Management Functions***********************************/ 
  void freeDenseMatMem();
  void freeSparseMatMem();
  /******************************** Check ******************************************************/
  void check_Structure(HODLR_Tree::node* HODLR_Root);
  void createExactHODLR(HODLR_Tree::node* HODLR_Root,const int rank,Eigen::MatrixXd & result);

};


#endif
