#include "user_IndexTree.hpp"

user_IndexTree::user_IndexTree(){
  rootNode = NULL;

}

user_IndexTree::~user_IndexTree(){
  if (rootNode != NULL)
    freeTree(rootNode);
  
}


void user_IndexTree::setChildren_NULL(node* parent){
  parent->left = NULL;
  parent->right = NULL;
}

void user_IndexTree::printTree() const{
  if (rootNode != NULL)
    printTree(rootNode);
  else
    std::cout<<"Error! user_IndexTree root node is NULL"<<std::endl;
}

void user_IndexTree::printTree(const node* root)const{
  if (root == NULL){
    return;
  }
  std::cout<<"Split Index = "<<root->splitIndex<<std::endl;
  std::cout<<"*********************"<<std::endl;
  printTree(root->left);
  printTree(root->right);
}

void user_IndexTree::freeTree(node* root){
  if (root == NULL)
    return;
  freeTree(root->left);
  freeTree(root->right);
  delete root;

}
