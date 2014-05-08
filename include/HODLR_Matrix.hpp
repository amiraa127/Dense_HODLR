#ifndef HODLR_MATRIX_HPP
#define HODLR_MATRIX_HPP

#include "HODLR_Tree.hpp"
#include "helperFunctions.hpp"
#include "user_IndexTree.hpp"
#include "recLU_FactorTree.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <cmath>
#include <vector>
#include <set>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>
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
  
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,int inputSizeThreshold);
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,int inputSizeThreshold);

  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,int inputSizeThreshold,user_IndexTree &input_IndexTree);
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,int inputSizeThreshold,user_IndexTree &input_IndexTree);
  
  HODLR_Matrix(const HODLR_Matrix & rhs); //Copy Constructor

  ~HODLR_Matrix();
  


  /************************************* Solve Methods **********************************/
  Eigen::MatrixXd recLU_Solve(const Eigen::MatrixXd & input_RHS);
  void recLU_Compute();
  Eigen::MatrixXd extendedSp_Solve(const Eigen::MatrixXd & input_RHS);
  Eigen::MatrixXd iterative_Solve(const Eigen::MatrixXd & input_RHS, const int maxIterations, const double stop_tolerance,const double init_LRTolerance,const std::string input_LR_Method, const std::string directSolve_Method);

  /**************************** Low-Rank Approximation Methods **************************/
  double partialPivACA_LowRankApprox(Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int max_i, const int min_j, const int max_j, const double tolerance, int & calculatedRank,const int minRank = -1);

  double fullPivACA_LowRankApprox(Eigen::MatrixXd & W,Eigen::MatrixXd & V, const int min_i, const int max_i, const int min_j, const int max_j, const double tolerance, int & calculatedRank,const int minRank = -1);
  
  void PS_LowRankApprox(Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int max_i, const int min_j, const int max_j, const double tolerance, int & calculatedRank, const std::string pointChoosingMethod = "Chebyshev",const int minRank = -1) const;
  
  void SVD_LowRankApprox(Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int max_i, const int min_j, const int max_j, const double tolerance, int & calculatedRank, const int minRank = -1) const;
  
  void PS_LowRankApprox_Sp(Eigen::MatrixXd & W, Eigen::MatrixXd & V, Eigen::MatrixXd & K,const int min_i, const int max_i,const int min_j, const int max_j, const double tolerance, int &calculatedRank)const;

  /************************************* Attribute Modification **********************************/
  void set_LRTolerance(double tolerance);  
  void set_MinValueACA(double minValue);
  void set_def_LRMethod(std::string input_LRMethod);
  void set_FreeMatrixMemory(bool inputVal);
  void saveExtendedSp(std::string savePath);
 
  /************************************ Accessing Attributes ************************************/
  double get_recLU_FactorizationTime() const;
  double get_recLU_SolveTime() const;
  double get_extendedSp_AssemblyTime() const; 
  double get_extendedSp_FactorizationTime() const;
  double get_extendedSp_SolveTime() const;
  double get_LR_ComputationTime() const;
  double get_iter_SolveTime() const;
  double get_MatrixSize() const;

  /************************************ Acessing HODLR Entries *******************************/
  Eigen::MatrixXd get_Block(int min_i,int min_j,int numRows,int numCols);
  Eigen::MatrixXd get_Row(int row);
  Eigen::MatrixXd get_Col(int col);
 
  
  HODLR_Matrix  topDiag();
  HODLR_Matrix bottDiag();

  void keepTopDiag();
  void keepBottDiag();
  /**********************************Extend Add Functions **************************************/
  void extend(std::vector<int> & extendIdxVec,int parentSize);
  void extendAddUpdate(Eigen::MatrixXd & D,std::vector<int> & updateIdxVec,std::string mode);
  void extendAddUpdate(Eigen::MatrixXd & updateExtendU,Eigen::MatrixXd & updateExtendV);

  Eigen::MatrixXd& returnTopOffDiagU();
  Eigen::MatrixXd& returnTopOffDiagV();
  Eigen::MatrixXd& returnTopOffDiagK();
  Eigen::MatrixXd& returnBottOffDiagU();
  Eigen::MatrixXd& returnBottOffDiagV();
  Eigen::MatrixXd& returnBottOffDiagK();

private:

  int sizeThreshold;
  int extendedSp_Size;
  int matrixSize;
  int matrixNumRows;
  int matrixNumCols;

  double recLU_FactorizationTime;
  double recLU_SolveTime;
  double LR_ComputationTime;
  double extendedSp_AssemblyTime;
  double extendedSp_FactorizationTime;
  double extendedSp_SolveTime;
  double iter_SolveTime;

  bool LRStoredInTree;
  bool createdRecLUfactorTree;
  bool assembled_ExtendedSp;
  bool saveExtendedSp_Matrix;
  bool freeMatrixMemory;
  bool freeMatrixMemory_Sp;
  bool matrixDataAvail;
  bool matrixDataAvail_Sp;
  bool isSquareMatrix;

  double LR_Tolerance;
  double minValueACA;

  HODLR_Tree indexTree;
  recLU_FactorTree recLUfactorTree;
  Eigen::MatrixXd matrixData;
  Eigen::SparseMatrix<double> matrixData_Sp;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > extendedSp_Solver;

  std::string extendedSp_SavePath;

  void setDefaultValues();
  void reset_attributes();

  /****************************recLU Solver Functions*******************************/
  void storeLRinTree();
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

  /***************************LR Approximation Related Functions********************/  
  int chooseNNZRowIndex(const std::vector<bool> &chosenRows) const;
  int chooseNextRowCol(const std::vector<bool> &chosenRowsCols, const Eigen::VectorXd &currColRow) const;
  
  void extractRowsCols(Eigen::MatrixXd & W, Eigen::MatrixXd & K, Eigen::MatrixXd & V, const Eigen::MatrixXd & inputMatrix,const std::vector<int> & rowIndex,const std::vector<int> & colIndex)const;
  
  void LUDecompose(const Eigen::MatrixXd &inputMatrix,Eigen::MatrixXd &LU,Eigen::MatrixXd &P) const;

  int SVD_LowRankApprox(const Eigen::MatrixXd & inputMatrix, const double accuracy, Eigen::MatrixXd* Wptr = NULL, Eigen::MatrixXd* Vptr = NULL, Eigen::MatrixXd* Kptr = NULL, int minRank = -1) const;
  
  /**********************************Accessing Matrix Entries***************************/
  void fill_Block(Eigen::MatrixXd & blkMatrix,HODLR_Tree::node* root,int min_i,int min_j,int max_i,int max_j);

  void fill_BlockWithLRProduct(Eigen::MatrixXd & blkMatrix,int LR_Min_i,int LR_Min_j, int LR_numRows, int LR_numCols,Eigen::MatrixXd & LR_U,Eigen::MatrixXd & LR_K,Eigen::MatrixXd & LR_V,int blk_Min_i,int blk_Min_j);

  /******************************Memory Management Functions***********************************/ 
  void freeDenseMatMem();
  void freeSparseMatMem();
  


  /***********************************Extend Add Functions *******************************/
  void extend(HODLR_Tree::node* HODLR_Root,std::vector<int> & extendIdxVec,int parentSize);
  void extendAddLRinTree(HODLR_Tree::node* HODLR_Root,const Eigen::MatrixXd & updateExtendU,const Eigen::MatrixXd & updateExtendV);
  void extendAddLRinTree(HODLR_Tree::node* HODLR_Root,HODLR_Matrix & extendD_HODLR,std::vector<int> & updateIdxVec);
  void extendAddLRinTree(HODLR_Tree::node* HODLR_Root,Eigen::MatrixXd & extendD,std::vector<int> & updateIdxVec);

  int add_LR(Eigen::MatrixXd & result_U,Eigen::MatrixXd & result_K,Eigen::MatrixXd & result_V,const Eigen::MatrixXd & U1, const Eigen::MatrixXd & V1, const Eigen::MatrixXd & U2, const Eigen::MatrixXd & V2);


};

Eigen::MatrixXd extend(std::vector<int> & extendIdxVec,int parentSize,Eigen::MatrixXd & child,int min_i,int min_j,int numRows,int numCols,std::string mode);


#endif
