#ifndef RECLU_FACTORTREE_HPP
#define RECLU_FACTORTREE_HPP

#include "Eigen/Dense"



class recLU_FactorTree{

public:
  struct node {
    node* left;
    node* right;
    bool isLeaf;
    Eigen::MatrixXd P_S, LU_S;
    Eigen::MatrixXd P_leaf, LU_leaf;
    Eigen::MatrixXd topDiagSoln_LR;
    Eigen::MatrixXd bottDiagSoln_LR;
  };

  recLU_FactorTree();
  ~recLU_FactorTree();

  node* rootNode;

private:
  void freeTree(node* root);

};


#endif
