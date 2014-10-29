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
    int topOffDiag_maxRank;
    int bottOffDiag_maxRank;
    std::string LR_Method;        
    node(){
      topOffDiag_minRank  = -1;
      bottOffDiag_minRank = -1; 
      topOffDiag_maxRank  = -1;
      bottOffDiag_maxRank = -1;
    }
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
