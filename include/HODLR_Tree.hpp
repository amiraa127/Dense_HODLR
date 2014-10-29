#ifndef HODLR_TREE_HPP
#define HODLR_TREE_HPP

#include <string>
#include <iostream>
#include <vector>
#include "assert.h"
#include "user_IndexTree.hpp"
#include "Eigen/Dense"

class HODLR_Tree{

public:
  struct node{
    node* left;
    node* right;
    int currLevel;
    int min_i;
    int min_j;
    int max_i;
    int max_j;
    int splitIndex_i;
    int splitIndex_j;
    int topOffDiag_minRank;
    int bottOffDiag_minRank;
    int topOffDiag_maxRank;
    int bottOffDiag_maxRank;
    bool isLeaf;
    std::string LR_Method;
    Eigen::MatrixXd topOffDiagU;
    Eigen::MatrixXd topOffDiagV;
    Eigen::MatrixXd bottOffDiagU;
    Eigen::MatrixXd bottOffDiagV;
    Eigen::MatrixXd leafMatrix;
    int topOffDiagRank;
    int bottOffDiagRank;
    std::vector<int> topOffDiagRowIdx;
    std::vector<int> topOffDiagColIdx;
    std::vector<int> bottOffDiagRowIdx;
    std::vector<int> bottOffDiagColIdx;
  };
  
  HODLR_Tree();
  HODLR_Tree(const int matrixSize, const std::string LRMethod = "ACA");
  ~HODLR_Tree();
  
  HODLR_Tree& operator = (const HODLR_Tree & rhs);

  
  void createDefaultTree(const int matrixSize);
  void createFromUsrTree(const int matrixSize,const user_IndexTree &usrIndexTree);
  void printTree() const; 

  void set_sizeThreshold(const int input_sizeThreshold);
  void set_LRMethod(const std::string input_LRMethod);

  int get_sizeThreshold()const;
  int get_numLevels()const;
  std::string get_def_LRMethod()const;

  node* rootNode;
  node* copyTree(std::vector<std::vector<node*> > &rhs_nodeLevelVec,std::vector<node*> &rhs_leafNodesVec) const;
  void freeTree(node* root);
  void correctIndices();
  void correctIndices(node *root,int offset_i,int offset_j);
  
  std::vector<std::vector<node*> > nodeLevelVec;
  std::vector<node*> leafNodesVec;

private: 
  
  int sizeThreshold;
  int numLevels;
  std::string def_LRMethod;
  
  void copyNode(const node* oldNode, node* newNode,std::vector<std::vector<node*> > &rhs_nodeLevelVec,std::vector<node*> &rhs_leafNodesVec)const;
  void createDefaultTree(node* root);
  void printTree(const node* root) const;
  void userTree_To_HODLRTree(const int currLevel,const int min_i,const int max_i,const int min_j,const int max_j,const user_IndexTree::node* user_IndexRoot, node* HODLR_IndexRoot);
  
  void nodeLevelVec_Assign(unsigned int level,node* root);
  void set_LRMethod(node* root,const std::string input_LRMethod);
};


#endif
