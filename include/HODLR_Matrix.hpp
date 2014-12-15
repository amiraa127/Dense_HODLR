#ifndef HODLR_MATRIX_HPP
#define HODLR_MATRIX_HPP

//Standard C++
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

//External Dependencies

#include <Eigen/Dense>
#include <Eigen/Sparse>

//Custom Dependencies
#include "helperFunctions.hpp"
#include "HODLR_Tree.hpp"
#include "kernel.hpp"
#include "lowRank.hpp"
#include "matrixIO.hpp"
#include "perturbI.hpp"
#include "user_IndexTree.hpp"
#include "recLU_FactorTree.hpp"

/**
 * \author Amirhossein Aminfar
 */

class HODLR_Matrix{

public:
 
  bool printLevelRankInfo;
  bool printLevelAccuracy;
  bool printLevelInfo;
  bool printResultInfo;
  
  
  /**
   *  \defgroup Constructors
   *  The HODLR matrix class provides a variety of options for initialization.
   *  I'm planning to merge some of these constructors into a single constructor. However, some constructors may become depreciated in the feature.
   *  @{
   */

  /**
   * \defgroup defConstruct Default Constructor
   * @{
   */
  HODLR_Matrix();
  /** @} */ 

  /**
   * \defgroup denseConstruct Dense Matrix Constructors
   * \brief Thse constructors build an HODLR matrix from a preallocated dense matrix. 
   * \note Currently, the dense matrix is being passed by a non const reference variable. So if you want your original matrix, copy it elsewhere befor passing it to the constructor.
   * \note Only use the partial pivoting ACA ("partialPiv_ACA"), full pivoting ACA ("fullPiv_ACA"), singualr value Decomposition ("SVD"), pseudo skeleton with Chebyshev point selection ("PS_Cheby") or BDLR ("PS_Boundary") as the low-rank approximation schemes for dense matrices.
   * @{
   */
  
  /**
   * \param[in] inputMatrix : Preallocated dense matrix as Eigen::MatrixXd matrix class.
   * \param[in] inputSizeThreshold : Leaf size threshold. If no value is provided, it will set the leaf size to the default value of 30.
   * \param[in] LR_Method : Low-rank approximation scheme to be used in calculating the off-diagonal low-rank approximations. If no value is provided, it will set the LR_Method parameter to "partialPiv_ACA".
   * \brief This constructor initializes the class with a dense matrix and an optional leaf size threshold. 
   */
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,int inputSizeThreshold = 30, std::string LR_Method = "partialPiv_ACA");

  /**
   * \param[in] inputMatrix : Preallocated dense matrix as Eigen::MatrixXd matrix class.
   * \param[in] inputSizeThreshold : Leaf size threshold.
   * \param[in] input_IndexTree : User defined splitting scheme stored as a user_IndexTree class. 
   * \brief This constructor initializes the class with a dense matrix and a user specified indexing scheme which will be used to create the HODLR index tree.
   * This constructor initializes the class with a dense matrix. 
   */
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,int inputSizeThreshold,user_IndexTree &input_IndexTree);
  
  /**
   * \param[in] inputMatrix : Preallocated dense matrix as Eigen::MatrixXd matrix class.
   * \param[in] inputGraph  : Preallocated sparse matrix. This sparse matrix will be used as the interaction graph in the BDLR low-rank approximation scheme.
   * \param[in] inputSizeThreshold : Leaf size threshold. If no value is provided, it will set the leaf size to the default value of 30.
   * \brief This constructor initializes the class with a dense matrix and an interaction graph (sparse matrix). The grpah is going to be used primarily as the interaction graph in the BDLR scheme.
   */
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold = 30);
  
  /**
   * \param[in] inputMatrix : Preallocated dense matrix as Eigen::MatrixXd matrix class.
   * \param[in] inputGraph  : Preallocated sparse matrix. This sparse matrix will be used as the interaction graph in the BDLR low-rank approximation scheme. 
   * \param[in] inputSizeThreshold : Leaf size threshold.
   * \param[in] input_IndexTree : User defined splitting scheme stored as a user_IndexTree class. 
   * \brief This constructor initializes the class with a dense matrix and an interaction graph (sparse matrix) and a user defined indexing schemes. The grpah is going to be used primarily as the interaction graph in the BDLR scheme. The user defined indexing scheme will be used to create the HODLR index tree.
   */
  HODLR_Matrix(Eigen::MatrixXd &inputMatrix,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold,user_IndexTree &input_IndexTree);
  
  
  /** @} */ 
  
  
  /**
   * \defgroup sparseConstruct Sparse Matrix Constructors
   * \brief Thse constructors build an HODLR matrix from a preallocated sparse matrix. 
   * \note Currently, the sparse matrix is being passed by a non const reference variable. So if you want your original matrix, copy it elsewhere befor passing it to the constructor.
   * \note Only use "PS_Sparse" or "identifyBoundary" as low-rank approximation methods when initializing the class with a sparse matrix.
   * @{
   */
  
  /**
   * \param[in] inputMatrix : Preallocated Sparse Matrix.
   * \param[in] inputSizeThreshold : Leaf size threshold. If no value is provided, it will set the leaf size to the default value of 30.
   * \param[in] LR_Method   : Low-rank approximation scheme to be used in calculating the off-diagonal low-rank approximations. If no value is provided, it will set the LR_Method parameter to "PS_Sparse".
   * \brief This constructor initializes the class with a sparse matrix and optional leaf size and low-rank approximation method parameters.
   */
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,int inputSizeThreshold = 30,std::string LR_Method = "PS_Sparse");

  /**
   * \param[in] inputMatrix : Preallocated Sparse Matrix.
   * \param[in] inputSizeThreshold : Leaf size threshold. If no value is provided, it will set the leaf size to the default value of 30.
   * \param[in] input_IndexTree : User defined splitting scheme stored as a user_IndexTree class. 
   * \param[in] LR_Method   : Low-rank approximation scheme to be used in calculating the off-diagonal low-rank approximations. If no value is provided, it will set the LR_Method parameter to "PS_Sparse".
   * \brief This constructor initializes the class with a sparse matrix and a user defined indexing scheme used to create the HODLR indexing tree.
   */
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,int inputSizeThreshold,user_IndexTree &input_IndexTree,std::string LR_Method = "PS_Sparse");
 

  /**
   * \param[in] inputMatrix : Preallocated Sparse Matrix.
   * \param[in] inputGraph  : Preallocated sparse matrix. This sparse matrix will be used as the interaction graph in the BDLR low-rank approximation scheme.
   * \param[in] inputSizeThreshold : Leaf size threshold. If no value is provided, it will set the leaf size to the default value of 30.
   * \param[in] input_IndexTree : User defined splitting scheme stored as a user_IndexTree class. 
   * \param[in] LR_Method   : Low-rank approximation scheme to be used in calculating the off-diagonal low-rank approximations. If no value is provided, it will set the LR_Method parameter to "PS_Sparse".
   * \brief This constructor initializes the class with a sparse matrix and an interaction graph (sparse matrix) and a user defined indexing schemes. The grpah is going to be used primarily as the interaction graph in the BDLR scheme. The user defined indexing scheme will be used to create the HODLR index tree.  
   */
  HODLR_Matrix(Eigen::SparseMatrix<double> &inputMatrix,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold,user_IndexTree &input_IndexTree,std::string LR_Method = "PS_Sparse");
  
  /** @} */ 


  /**
   * \defgroup kernelConstruct Kernel Matrix Constructors
   * \brief These constructors build an HODLR matrix from a kernel function. 
   * @{
   */
  

  /**
   * \param[in] numRows : Number of rows of the dense matrix.
   * \param[in] numCols : Number of columns of the dense matrix.
   * \param[in] inputKernel: Kernel function defining the matrix.
   * \param[in] inputKernelData: Additional information needed by the kernel function in excess of the row and column index of the entry.
   * \param[in] inputSizeThreshold: Leaf size threshold. If no value is provided, it will set the leaf size to the default value of 30.
   * \brief This constructor initializes the class with a kernel function and a size threshold. 
   */
  HODLR_Matrix(int numRows, int numCols,double (*inputKernel)(int i,int j,void* inputKernelData),void* inputKernelData,int inputSizeThreshold = 30);
  
    /**
   * \param[in] numRows : Number of rows of the dense matrix.
   * \param[in] numCols : Number of columns of the dense matrix.
   * \param[in] inputKernel: Kernel function defining the matrix.
   * \param[in] inputKernelData: Additional information needed by the kernel function in excess of the row and column index of the entry.
   * \param[in] inputGraph  : Preallocated sparse matrix. This sparse matrix will be used as the interaction graph in the BDLR low-rank approximation scheme.
   * \param[in] inputSizeThreshold: Leaf size threshold. If no value is provided, it will set the leaf size to the default value of 30.
   * \brief This constructor initializes the class with a kernel function, a size threshold and an interaction graph.  The grpah is going to be used primarily as the interaction graph in the BDLR scheme. 
   */
  HODLR_Matrix(int numRows, int numCols,double (*inputKernel)(int i,int j,void* inputKernelData),void* inputKernelData,Eigen::SparseMatrix<double> &inputGraph,int inputSizeThreshold = 30);
  
  /**
   * \param[in] numRows : Number of rows of the dense matrix.
   * \param[in] numCols : Number of columns of the dense matrix.
   * \param[in] inputKernel: Kernel function defining the matrix.
   * \param[in] inputKernelData: Additional information needed by the kernel function in excess of the row and column index of the entry.
   * \param[in] inputSizeThreshold: Leaf size threshold.
   * \param[in] input_IndexTree : User defined splitting scheme stored as a user_IndexTree class. 
   * \brief This constructor initializes the class with a kernel function a user defined indexing scheme used to create the HODLR indexing tree. 
   */
  HODLR_Matrix(int numRows, int numCols,double (*inputKernel)(int i,int j,void* inputKernelData),void* inputKernelData,int inputSizeThreshold,user_IndexTree &input_IndexTree);

  /**
   * \param[in] numRows : Number of rows of the dense matrix.
   * \param[in] numCols : Number of columns of the dense matrix.
   * \param[in] inputKernel: Kernel function defining the matrix.
   * \param[in] inputKernelData: Additional information needed by the kernel function in excess of the row and column index of the entry.
   * \param[in] inputGraph  : Preallocated sparse matrix. This sparse matrix will be used as the interaction graph in the BDLR low-rank approximation scheme.
   * \param[in] inputSizeThreshold: Leaf size threshold.
   * \param[in] input_IndexTree : User defined splitting scheme stored as a user_IndexTree class. 
   * \brief This constructor initializes the class with a kernel function and an interaction graph (sparse matrix) and a user defined indexing schemes. The grpah is going to be used primarily as the interaction graph in the BDLR scheme. The user defined indexing scheme will be used to create the HODLR index tree.  
   */  
  HODLR_Matrix(int numRows, int numCols,double (*inputKernel)(int i,int j,void* inputKernelData),void* inputKernelData, Eigen::SparseMatrix<double> &inputGraph, int inputSizeThreshold, user_IndexTree &input_IndexTree);
  /** @} */

 
  /**
   * \param[in] inputMatrix
   * This constructor initializes the class with a dense matrix. 
   */
  HODLR_Matrix(const HODLR_Matrix & rhs); //Copy Constructor
  /** @} */

  
  ~HODLR_Matrix();
  
  /************************************* Create HODLR Structure ***************************************/
    
  void storeLRinTree();
  
  /************************************* Solve Methods **********************************/
  Eigen::MatrixXd recLU_Solve(const Eigen::MatrixXd & input_RHS);
  Eigen::MatrixXd recSM_Solve(const Eigen::MatrixXd & input_RHS);

  void recLU_Compute();
  Eigen::MatrixXd extendedSp_Solve(const Eigen::MatrixXd & input_RHS);
  Eigen::MatrixXd iterative_Solve(const Eigen::MatrixXd & input_RHS, const int maxIterations, const double stop_tolerance,const double init_LRTolerance,const std::string input_LR_Method, const std::string directSolve_Method);

  /************************************* Attribute Modification **********************************/
  void set_LRTolerance(double tolerance);  
  void set_minPivot(double minPivot);
  void set_pastix_MinPivot(double minPivot);
  void set_LRMethod(std::string input_LRMethod);
  void set_FreeMatrixMemory(bool inputVal);
  void set_BoundaryDepth(int inputBoundaryDepth);  
  void set_recLUFactorizedFlag(bool factorized);
  void set_numSel(int numSel_);
  void set_LeafConst();
  
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
  double determinant();
  double logAbsDeterminant();

  /************************************ Acessing HODLR Entries *******************************/
  Eigen::MatrixXd block(int min_i,int min_j,int numRows,int numCols);
  Eigen::MatrixXd row(int row);
  Eigen::MatrixXd col(int col);
 
  HODLR_Matrix  topDiag();
  HODLR_Matrix bottDiag();

  void keepTopDiag();
  void keepBottDiag();
  friend void splitAtTop(HODLR_Matrix& self,HODLR_Matrix& topHODLR, HODLR_Matrix& bottHODLR);


  Eigen::MatrixXd& returnTopOffDiagU();
  Eigen::MatrixXd& returnTopOffDiagV();
  Eigen::MatrixXd& returnTopOffDiagK();
  Eigen::MatrixXd& returnBottOffDiagU();
  Eigen::MatrixXd& returnBottOffDiagV();
  Eigen::MatrixXd& returnBottOffDiagK();

  /******************************** Check ******************************************************/
  void check_Structure();
  double calcAbsDiff();
  Eigen::MatrixXd createExactHODLR(const int rank,int input_MatrixSize,const int inpt_SizeThreshold);

  void saveSolverInfo(const std::string outputFileName);
  
  void freeMatrixData();
  void destroyAllData();
  void recalculateSize();
  void correctIndices();
  void initInfoVectors();

private:

  int sizeThreshold;
  int extendedSp_Size;
  int matrixSize;
  int matrixNumRows;
  int matrixNumCols;
  int constLeafSize;

  double recLU_FactorizationTime;
  double recLU_SolveTime;
  double recLU_TotalTime;
  double recSM_FactorizationTime;
  double recSM_SolveTime;
  double recSM_TotalTime;

  double LR_ComputationTime;
  double extendedSp_AssemblyTime;
  double extendedSp_FactorizationTime;
  double extendedSp_SolveTime;
  double extendedSp_TotalTime;
  double totalIter_SolveTime;
  double determinant_;
  double logAbsDeterminant_;

  std::vector<double> LR_ComputationLevelTimeVec;
  std::vector<double> recLU_FactorLevelTimeVec;
  std::vector<double> recLU_SolveLevelTimeVec;
  std::vector<double> iter_IterTimeVec;
  std::vector<double> iter_IterAccuracyVec;
  std::vector<double> levelRankAverageVec;
  std::vector<double> levelRankAverageVecCnt;
  
  bool LRStoredInTree;
  bool recLU_Factorized;
  bool recSM_Factorized;
  bool assembled_ExtendedSp;
  bool saveExtendedSp_Matrix;
  bool freeMatrixMemory;
  bool freeMatrixMemory_Sp;
  bool freeGraphMemmory;
  bool matrixDataAvail;
  bool matrixDataAvail_Sp;
  bool graphDataAvail;
  bool kernelDataAvail;
  bool isSquareMatrix;
  bool isLeafConst;
  bool constLeafSet;
  bool constLeafFactorized;
  bool calculatedDet;

  double LR_Tolerance;
  double minPivot;
  double pastix_MinPivot;
  int boundaryDepth;
  int numSel;
  
  HODLR_Tree indexTree;
  recLU_FactorTree recLUfactorTree;
  Eigen::MatrixXd matrixData;
  Eigen::SparseMatrix<double> matrixData_Sp;
  Eigen::SparseMatrix<double> graphData;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > extendedSp_Solver;
  Eigen::MatrixXd constLeaf;
  Eigen::PartialPivLU<Eigen::MatrixXd> constLeafLU;
  kernelMatrix kernelMatrixData;

  std::string extendedSp_SavePath;

  void setDefaultValues();
  void initialize(Eigen::MatrixXd& inputMatrix);
  void initialize(Eigen::SparseMatrix<double>& inputMatrix);
  void initialize(int numRows, int numCols,double (*inputKernel)(int i,int j,void* kernelData),void* inputKernelData);
  void reset_attributes();
  void initializeInfoVecotrs(int numLevels);

  /****************************recLU Solver Functions*******************************/
 
  void storeLRinTree(HODLR_Tree::node* HODLR_Root);
  void recLU_Factorize();
  void recSM_Factorize();
  Eigen::MatrixXd recLU_Factorize(const Eigen::MatrixXd & input_RHS,const HODLR_Tree::node* HODLR_Root, recLU_FactorTree::node* factorRoot);
  void recSM_Factorize(HODLR_Tree::node* HODLR_Root,std::vector<HODLR_Tree::node*> &leftChildren, std::vector<HODLR_Tree::node*> &rightChildren,int desLevel);

  Eigen::MatrixXd recLU_Solve(const Eigen::MatrixXd & input_RHS,const HODLR_Tree::node* HODLR_Root, const recLU_FactorTree::node* factorRoot);

  void recSM_Solve(HODLR_Tree::node* HODLR_Root,Eigen::MatrixXd &RHS);
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
  
  /******************************** Determinant ************************************************/    void calcDeterminant();
  void calcDeterminant(HODLR_Tree::node* HODLR_Root);

  /******************************** Friend Functions *******************************************/
  friend void extendAddUpdate(HODLR_Matrix & parentHODLR, std::vector<Eigen::MatrixXd*> D_Array,std::vector<HODLR_Matrix*> D_HODLR_Array,std::vector<Eigen::MatrixXd*> U_Array,std::vector<Eigen::MatrixXd*> V_Array,std::vector<std::vector<int> > & updateIdxVec_Array_D,std::vector<std::vector<int> > & updateIdxVec_Array_D_HODLR,double tol,std::string mode);
  friend void extendAddUpdate(HODLR_Matrix & parentHODLR,HODLR_Matrix & D_HODLR,std::vector<int> & updateIdxVec,double tol,std::string mode);
  friend void extendAddUpdate(HODLR_Matrix & parentHODLR,Eigen::MatrixXd & U,Eigen::MatrixXd & V,std::vector<int>& updateIdxVec,double tol,std::string mode);
  friend HODLR_Matrix extend(std::vector<int> & extendIdxVec, int parentSize, HODLR_Matrix & childHODLR);
  friend void extendAddUpdate(HODLR_Matrix & parentHODLR,Eigen::MatrixXd & D,std::vector<int> & updateIdxVec,double tol,std::string mode);
};

/** \class HODLR_Matrix HODLR_Matrix.hpp 
 *  \brief This is the main HODLR class that includes all the fast HODLR solvers.
 *
 *  This class can be used as a sibling of the Eigen::MatrixXd class in cases where access to matrix entries is required.
 */


#endif
