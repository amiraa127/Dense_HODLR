#include "recLU_FactorTree.hpp"

recLU_FactorTree::recLU_FactorTree(){
  rootNode = NULL;
}

recLU_FactorTree::~recLU_FactorTree(){
  if (rootNode !=NULL)
    freeTree(rootNode);

}

void recLU_FactorTree::freeTree(node* root){
  if (root->isLeaf == true){
    delete root;
    return;
  }
  freeTree(root->left);
  freeTree(root->right);
  delete root;

}
