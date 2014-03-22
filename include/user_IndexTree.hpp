#ifndef USER_INDEXTREE_HPP
#define USER_INDEXTREE_HPP

#include <string>
#include <iostream>

class user_IndexTree{                                                                     

public:                                                                                  
  struct node{                                                                           
    node* left;
    node* right;
    int splitIndex; 
    int topOffDiag_minRank;
    int bottOffDiag_minRank;
    std::string LR_Method;        
  };
  
  node* rootNode;
  
  user_IndexTree();
  ~user_IndexTree();
  
  void setChildren_NULL(node* parent);
  void printTree() const;

private:
  void printTree(const node* root)const;
  void freeTree(node* root);
};   

#endif
