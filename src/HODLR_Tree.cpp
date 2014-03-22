#include "HODLR_Tree.hpp"

HODLR_Tree::HODLR_Tree(){
  rootNode = NULL;
  def_LRMethod = "partialPiv_ACA";
  numLevels = 0;
  sizeThreshold = 30;
}

HODLR_Tree::HODLR_Tree(const int matrixSize, const std::string LRMethod){
  rootNode = NULL;
  numLevels = 0;
  sizeThreshold = 30;
  def_LRMethod = LRMethod;
  createDefaultTree(matrixSize);
}
  
HODLR_Tree::~HODLR_Tree(){
  if (rootNode != NULL)
    freeTree(rootNode);
}

void HODLR_Tree::createDefaultTree(const int matrixSize){
  numLevels = 0;
  rootNode = new node;
  rootNode->min_i = 0;
  rootNode->max_i = matrixSize - 1;
  rootNode->min_j = 0;
  rootNode->max_j = matrixSize - 1;
  rootNode->splitIndex_i = matrixSize/2 - 1;
  rootNode->splitIndex_j = matrixSize/2 - 1;
  rootNode->LR_Method = def_LRMethod;
  rootNode->currLevel = 0;
  numLevels = 1;
  rootNode->topOffDiag_minRank = -1;
  rootNode->bottOffDiag_minRank = -1;

  if ( matrixSize < sizeThreshold ){
    rootNode->isLeaf = true;
    return;
  }
  rootNode->isLeaf = false;
  createDefaultTree(rootNode);  
}

void HODLR_Tree::createDefaultTree(node* root){
  
  int topDiagSize = (root->splitIndex_i - root->min_i + 1);
  int bottDiagSize = (root->max_i - root->splitIndex_i);

  if (root->currLevel > (numLevels - 1))
     numLevels = root->currLevel + 1;
 
  node* leftNode  = new node;
  node* rightNode = new node;

  // Allocate left node 
  leftNode->min_i = root->min_i;
  leftNode->max_i = root->min_i + topDiagSize - 1;
  leftNode->min_j = root->min_j;
  leftNode->max_j = root->min_j + topDiagSize - 1;
  leftNode->LR_Method = root->LR_Method;
  leftNode->currLevel = root->currLevel + 1;
  leftNode->topOffDiag_minRank = -1;
  leftNode->bottOffDiag_minRank = -1;

  root->left = leftNode;
 
  if (topDiagSize > sizeThreshold){
    leftNode->isLeaf = false;
    leftNode->splitIndex_i = root->min_i + topDiagSize/2 - 1;
    leftNode->splitIndex_j = root->min_j + topDiagSize/2 - 1;
    createDefaultTree(leftNode);
  }else{
    leftNode->isLeaf = true;
    leftNode->splitIndex_i = -1;
    leftNode->splitIndex_j = -1;
    leftNode->left = NULL;
    leftNode->right = NULL;
  }

  // Allocate right node
  rightNode->min_i = root->min_i + topDiagSize;
  rightNode->max_i = root->max_i;
  rightNode->min_j = root->min_j + topDiagSize;
  rightNode->max_j = root->max_j;
  rightNode->LR_Method = root->LR_Method;
  rightNode->currLevel = root->currLevel + 1;
  rightNode->topOffDiag_minRank = -1;
  rightNode->bottOffDiag_minRank = -1;

  root->right = rightNode;
  
  if (bottDiagSize > sizeThreshold){
    rightNode->isLeaf = false;
    rightNode->splitIndex_i = root->max_i - bottDiagSize/2;
    rightNode->splitIndex_j = root->max_j - bottDiagSize/2;
    createDefaultTree(rightNode);
  }else{
    rightNode->isLeaf = true;
    rightNode->splitIndex_i = -1;
    rightNode->splitIndex_j = -1;
    leftNode->left = NULL;
    leftNode->right = NULL;
  }
  return;
}

void HODLR_Tree::printTree() const{
  if (rootNode != NULL)
    printTree(rootNode);
  else
    std::cout<<"Error! HODLR tree root node is NULL."<<std::endl;
}

void HODLR_Tree::printTree(const node* root) const{
  if (root->isLeaf == true){
    std::cout<<"level = "<<root->currLevel<<" (leaf)"<<std::endl;
    std::cout<<"min_i = "<<root->min_i<<" -- "<<"splitIndex_i = "<<root->splitIndex_i<<" -- "<<"max_i = "<<root->max_i<<std::endl;
    std::cout<<"min_j = "<<root->min_j<<" -- "<<"splitIndex_j = "<<root->splitIndex_j<<" -- "<<"max_j = "<<root->max_j<<std::endl;
    std::cout<<"**************************************************"<<std::endl;
    return;
  }
  printTree(root->left);
  printTree(root->right);
  std::cout<<"level = "<<root->currLevel<<std::endl;
  std::cout<<"min_i = "<<root->min_i<<" -- "<<"splitIndex_i = "<<root->splitIndex_i<<" -- "<<"max_i = "<<root->max_i<<std::endl;
  std::cout<<"min_j = "<<root->min_j<<" -- "<<"splitIndex_j = "<<root->splitIndex_j<<" -- "<<"max_j = "<<root->max_j<<std::endl;
  std::cout<<"**************************************************"<<std::endl;
}

void HODLR_Tree::freeTree(node* root){
  if (root->isLeaf == true){
    delete root;
    return;
  }
  freeTree(root->left);
  freeTree(root->right);
  delete root;
}

void HODLR_Tree::userTree_To_HODLRTree(const int currLevel,const int min_i,const int max_i,const int min_j,const int max_j,const user_IndexTree::node* user_IndexRoot, node* HODLR_IndexRoot){   

  int matrixSize = max_i - min_i + 1;
 
  if (currLevel > (numLevels - 1))
    numLevels = currLevel + 1;
 
  if (matrixSize <= sizeThreshold){
    HODLR_IndexRoot->isLeaf = true;
    HODLR_IndexRoot->min_i = min_i;
    HODLR_IndexRoot->min_j = min_j;
    HODLR_IndexRoot->max_i = max_i;
    HODLR_IndexRoot->max_j = max_j;
    HODLR_IndexRoot->left = NULL;
    HODLR_IndexRoot->right = NULL;
    HODLR_IndexRoot->currLevel = currLevel;
    HODLR_IndexRoot->topOffDiag_minRank = user_IndexRoot->topOffDiag_minRank;
    HODLR_IndexRoot->bottOffDiag_minRank = user_IndexRoot->bottOffDiag_minRank;
    HODLR_IndexRoot->LR_Method = user_IndexRoot->LR_Method;
    HODLR_IndexRoot->splitIndex_i = -1;
    HODLR_IndexRoot->splitIndex_j = -1;
    return;
  }
    HODLR_IndexRoot->isLeaf = false;
    HODLR_IndexRoot->min_i = min_i;
    HODLR_IndexRoot->min_j = min_j;
    HODLR_IndexRoot->max_i = max_i;
    HODLR_IndexRoot->max_j = max_j;
    HODLR_IndexRoot->currLevel = currLevel;
    HODLR_IndexRoot->topOffDiag_minRank = user_IndexRoot->topOffDiag_minRank;
    HODLR_IndexRoot->bottOffDiag_minRank = user_IndexRoot->bottOffDiag_minRank;
    HODLR_IndexRoot->LR_Method = user_IndexRoot->LR_Method;
    HODLR_IndexRoot->splitIndex_i = user_IndexRoot->splitIndex;
    HODLR_IndexRoot->splitIndex_j = user_IndexRoot->splitIndex;
    
    int leftNodeSize = user_IndexRoot->splitIndex - min_i + 1;
    int rightNodeSize = max_i - user_IndexRoot->splitIndex;
    assert(leftNodeSize > 0);
    assert(rightNodeSize > 0);

    //handle children
    node* leftNode = new node;  
    node* rightNode = new node;  
    HODLR_IndexRoot->left = leftNode;
    HODLR_IndexRoot->right = rightNode;
    if ((user_IndexRoot->left == NULL) && (user_IndexRoot->right ==NULL)){
      //create default left node
      leftNode->currLevel = currLevel + 1;
      leftNode->min_i = min_i;
      leftNode->min_j = min_j;
      leftNode->max_i = user_IndexRoot->splitIndex;
      leftNode->max_j = user_IndexRoot->splitIndex;
      leftNode->topOffDiag_minRank = -1;
      leftNode->bottOffDiag_minRank = -1;
      leftNode->LR_Method = def_LRMethod;
      if (leftNodeSize <= sizeThreshold){
	leftNode->isLeaf = true;
	leftNode->splitIndex_i = -1;
	leftNode->splitIndex_j = -1;
	leftNode->left = NULL;
	leftNode->right = NULL;
      }else{
	leftNode->isLeaf = false;
	leftNode->splitIndex_i = min_i + leftNodeSize/2 - 1;
	leftNode->splitIndex_j = min_j + leftNodeSize/2 - 1;
	createDefaultTree(leftNode);
      }
      //create default right node
      rightNode->currLevel = currLevel + 1;
      rightNode->min_i = user_IndexRoot->splitIndex + 1;
      rightNode->min_j = user_IndexRoot->splitIndex + 1;
      rightNode->max_i = max_i;
      rightNode->max_j = max_j;
      rightNode->topOffDiag_minRank = -1;
      rightNode->bottOffDiag_minRank = -1;
      rightNode->LR_Method = def_LRMethod;
      if (rightNodeSize <= sizeThreshold){
	rightNode->isLeaf = true;
	rightNode->splitIndex_i = -1;
	rightNode->splitIndex_j = -1;
	rightNode->left = NULL;
	rightNode->right = NULL;
      }else{
	rightNode->isLeaf = false;
	rightNode->splitIndex_i = user_IndexRoot->splitIndex + rightNodeSize/2;
	rightNode->splitIndex_j = user_IndexRoot->splitIndex + rightNodeSize/2;
	createDefaultTree(rightNode);
      }
    }else{
      userTree_To_HODLRTree(currLevel + 1,min_i,user_IndexRoot->splitIndex,min_j,user_IndexRoot->splitIndex,user_IndexRoot->left,HODLR_IndexRoot->left); 
      userTree_To_HODLRTree(currLevel + 1,user_IndexRoot->splitIndex+1,max_i,user_IndexRoot->splitIndex + 1,max_j,user_IndexRoot->right,HODLR_IndexRoot->right); 
 }   
}
  

void HODLR_Tree::createFromUsrTree(const int matrixSize,const user_IndexTree &usrIndexTree){
  if (usrIndexTree.rootNode == NULL)
    createDefaultTree(matrixSize);
  else{
    rootNode = new node;
    numLevels = 1;
    userTree_To_HODLRTree(0,0,matrixSize - 1,0,matrixSize - 1,usrIndexTree.rootNode,rootNode);
  }
}


void HODLR_Tree::set_sizeThreshold(const int input_sizeThreshold){
  if (rootNode == NULL)
    sizeThreshold = input_sizeThreshold;
  else{
    std::cout<<"Error! Cannot perform this operation at this time."<<std::endl;
    exit(EXIT_FAILURE);
  }
}

void HODLR_Tree::set_def_LRMethod(const std::string input_LRMethod){
  def_LRMethod = input_LRMethod;
  return;
}

int HODLR_Tree::get_sizeThreshold()const{
  return sizeThreshold;
}

int HODLR_Tree::get_numLevels()const{
  return numLevels;
}

std::string HODLR_Tree::get_def_LRMethod()const{
  return def_LRMethod;
}
  
