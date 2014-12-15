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

HODLR_Tree& HODLR_Tree::operator = (const HODLR_Tree & rhs){

  sizeThreshold = rhs.sizeThreshold;
  numLevels     = rhs.numLevels;
  def_LRMethod  = rhs.def_LRMethod;
  nodeLevelVec.resize(rhs.nodeLevelVec.size());
  rootNode = rhs.copyTree(nodeLevelVec,leafNodesVec);
  return *this;
}

void HODLR_Tree::createDefaultTree(const int matrixSize){
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
  rootNode->topOffDiag_minRank  = -1;
  rootNode->bottOffDiag_minRank = -1;
  rootNode->topOffDiag_maxRank  = -1;
  rootNode->bottOffDiag_maxRank = -1;

  if ( matrixSize < sizeThreshold ){
    rootNode->isLeaf = true;
    leafNodesVec.push_back(rootNode);
  }else{
    rootNode->isLeaf = false;
    nodeLevelVec_Assign(rootNode->currLevel,rootNode);
    createDefaultTree(rootNode);  
  }
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
  leftNode->topOffDiag_minRank  = -1;
  leftNode->bottOffDiag_minRank = -1;
  leftNode->topOffDiag_maxRank  = -1;
  leftNode->bottOffDiag_maxRank = -1;

  root->left = leftNode;
  
  if (topDiagSize > sizeThreshold){
    leftNode->isLeaf = false;
    leftNode->splitIndex_i = root->min_i + topDiagSize/2 - 1;
    leftNode->splitIndex_j = root->min_j + topDiagSize/2 - 1;
    nodeLevelVec_Assign(leftNode->currLevel,leftNode);
    createDefaultTree(leftNode);
  }else{
    leftNode->isLeaf = true;
    leftNode->splitIndex_i = -1;
    leftNode->splitIndex_j = -1;
    leftNode->left = NULL;
    leftNode->right = NULL;
    leafNodesVec.push_back(leftNode);
    if (leftNode->currLevel  > (numLevels - 1))
      numLevels = leftNode->currLevel + 1;
  }

  // Allocate right node
  rightNode->min_i = root->min_i + topDiagSize;
  rightNode->max_i = root->max_i;
  rightNode->min_j = root->min_j + topDiagSize;
  rightNode->max_j = root->max_j;
  rightNode->LR_Method = root->LR_Method;
  rightNode->currLevel = root->currLevel + 1;
  rightNode->topOffDiag_minRank  = -1;
  rightNode->bottOffDiag_minRank = -1;
  rightNode->topOffDiag_maxRank  = -1;
  rightNode->bottOffDiag_maxRank = -1;

  root->right = rightNode;
  
  if (bottDiagSize > sizeThreshold){
    rightNode->isLeaf = false;
    rightNode->splitIndex_i = root->max_i - bottDiagSize/2;
    rightNode->splitIndex_j = root->max_j - bottDiagSize/2;
    nodeLevelVec_Assign(rightNode->currLevel,rightNode);
    createDefaultTree(rightNode);
  }else{
    rightNode->isLeaf = true;
    rightNode->splitIndex_i = -1;
    rightNode->splitIndex_j = -1;
    rightNode->left = NULL;
    rightNode->right = NULL;
    leafNodesVec.push_back(rightNode);
    if (rightNode->currLevel  > (numLevels - 1))
      numLevels = rightNode->currLevel + 1;
  }
  return;
}

void HODLR_Tree::nodeLevelVec_Assign(unsigned int level,node* root){
  if (nodeLevelVec.size() > level)
    nodeLevelVec[level].push_back(root);
  else{
    nodeLevelVec.resize(level + 1);
    nodeLevelVec[level].push_back(root);
  }
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
    root = NULL;
    return;
  }
  freeTree(root->left);
  freeTree(root->right);
  delete root;
  root = NULL;
}


HODLR_Tree::node* HODLR_Tree::copyTree(std::vector<std::vector<node*> > &rhs_nodeLevelVec,std::vector<node*> &rhs_leafNodesVec ) const{
  if (rootNode == NULL)
    return NULL;
  node *cpyRoot = new node;
  copyNode(rootNode,cpyRoot,rhs_nodeLevelVec,rhs_leafNodesVec);
  return cpyRoot;
}

void HODLR_Tree::copyNode(const node* oldNode, node* newNode, std::vector<std::vector<node*> > &rhs_nodeLevelVec,std::vector<node*> &rhs_leafNodesVec)const{
  *(newNode) = *(oldNode);
  if (oldNode->isLeaf == true){
    newNode->left = NULL;
    newNode->right = NULL; 
    rhs_leafNodesVec.push_back(newNode);
    return;
  }
  else{
    node* leftCpy = new node;
    node* rightCpy = new node;
    newNode->left = leftCpy;
    newNode->right = rightCpy;
    rhs_nodeLevelVec[newNode->currLevel].push_back(newNode);
    copyNode(oldNode->left,leftCpy,rhs_nodeLevelVec,rhs_leafNodesVec);
    copyNode(oldNode->right,rightCpy,rhs_nodeLevelVec,rhs_leafNodesVec);
  }
}

void HODLR_Tree::userTree_To_HODLRTree(const int currLevel,const int min_i,const int max_i,const int min_j,const int max_j,const user_IndexTree::node* user_IndexRoot, node* HODLR_IndexRoot){   
  
  int matrixSize = max_i - min_i + 1;
  
  if (currLevel > (numLevels - 1))
    numLevels = currLevel + 1;
  
  if ((matrixSize <= sizeThreshold) && (user_IndexRoot->splitIndex == -1)){
 
    HODLR_IndexRoot->isLeaf = true;
    HODLR_IndexRoot->min_i = min_i;
    HODLR_IndexRoot->min_j = min_j;
    HODLR_IndexRoot->max_i = max_i;
    HODLR_IndexRoot->max_j = max_j;
    HODLR_IndexRoot->left = NULL;
    HODLR_IndexRoot->right = NULL;
    HODLR_IndexRoot->currLevel = currLevel;
    HODLR_IndexRoot->topOffDiag_minRank  = user_IndexRoot->topOffDiag_minRank;
    HODLR_IndexRoot->bottOffDiag_minRank = user_IndexRoot->bottOffDiag_minRank;
    HODLR_IndexRoot->topOffDiag_maxRank  = user_IndexRoot->topOffDiag_maxRank;
    HODLR_IndexRoot->bottOffDiag_maxRank = user_IndexRoot->bottOffDiag_maxRank;
    HODLR_IndexRoot->LR_Method = user_IndexRoot->LR_Method;
    HODLR_IndexRoot->splitIndex_i = -1;
    HODLR_IndexRoot->splitIndex_j = -1;
    leafNodesVec.push_back(HODLR_IndexRoot);
    return;
  }
  
  HODLR_IndexRoot->isLeaf = false;
  HODLR_IndexRoot->min_i  = min_i;
  HODLR_IndexRoot->min_j  = min_j;
  HODLR_IndexRoot->max_i  = max_i;
  HODLR_IndexRoot->max_j  = max_j;
  HODLR_IndexRoot->currLevel = currLevel;
  HODLR_IndexRoot->topOffDiag_minRank  = user_IndexRoot->topOffDiag_minRank;
  HODLR_IndexRoot->bottOffDiag_minRank = user_IndexRoot->bottOffDiag_minRank;
  HODLR_IndexRoot->topOffDiag_maxRank  = user_IndexRoot->topOffDiag_maxRank;
  HODLR_IndexRoot->bottOffDiag_maxRank = user_IndexRoot->bottOffDiag_maxRank;
  
  HODLR_IndexRoot->LR_Method = user_IndexRoot->LR_Method;
  HODLR_IndexRoot->splitIndex_i = user_IndexRoot->splitIndex;
  HODLR_IndexRoot->splitIndex_j = user_IndexRoot->splitIndex;
  nodeLevelVec_Assign(HODLR_IndexRoot->currLevel,HODLR_IndexRoot);
  
  int leftNodeSize = user_IndexRoot->splitIndex - min_i + 1;
  int rightNodeSize = max_i - user_IndexRoot->splitIndex;
  
  assert(leftNodeSize > 0);
  assert(rightNodeSize > 0);
  
  //handle children
  node* leftNode = new node;  
  node* rightNode = new node;  
  HODLR_IndexRoot->left = leftNode;
  HODLR_IndexRoot->right = rightNode;
  
  if ((user_IndexRoot->left == NULL) && (user_IndexRoot->right == NULL)){
    //create default left node
    leftNode->currLevel = currLevel + 1;
    leftNode->min_i = min_i;
    leftNode->min_j = min_j;
    leftNode->max_i = user_IndexRoot->splitIndex;
    leftNode->max_j = user_IndexRoot->splitIndex;
    leftNode->topOffDiag_minRank  = -1;
    leftNode->bottOffDiag_minRank = -1;
    leftNode->topOffDiag_maxRank  = -1;
    leftNode->bottOffDiag_maxRank = -1;
    
    leftNode->LR_Method = def_LRMethod;
    if (leftNodeSize <= sizeThreshold){
      leftNode->isLeaf = true;
      leftNode->splitIndex_i = -1;
      leftNode->splitIndex_j = -1;
      leftNode->left = NULL;
      leftNode->right = NULL;
      leafNodesVec.push_back(leftNode);
      if (leftNode->currLevel  > (numLevels - 1))
	numLevels = leftNode->currLevel + 1;
    }else{
      leftNode->isLeaf = false;
      leftNode->splitIndex_i = min_i + leftNodeSize/2 - 1;
      leftNode->splitIndex_j = min_j + leftNodeSize/2 - 1;
      nodeLevelVec_Assign(leftNode->currLevel,leftNode);
      createDefaultTree(leftNode);
    }
    //create default right node
    rightNode->currLevel = currLevel + 1;
    rightNode->min_i = user_IndexRoot->splitIndex + 1;
    rightNode->min_j = user_IndexRoot->splitIndex + 1;
    rightNode->max_i = max_i;
    rightNode->max_j = max_j;
    rightNode->topOffDiag_minRank  = -1;
    rightNode->bottOffDiag_minRank = -1;
    rightNode->topOffDiag_maxRank  = -1;
    rightNode->bottOffDiag_maxRank = -1;
    
    
    rightNode->LR_Method = def_LRMethod;
    if (rightNodeSize <= sizeThreshold){
      rightNode->isLeaf = true;
      rightNode->splitIndex_i = -1;
      rightNode->splitIndex_j = -1;
      rightNode->left = NULL;
      rightNode->right = NULL;
      leafNodesVec.push_back(rightNode);
      if (rightNode->currLevel  > (numLevels - 1))
	numLevels = rightNode->currLevel + 1;
    }else{
      rightNode->isLeaf = false;
      rightNode->splitIndex_i = user_IndexRoot->splitIndex + rightNodeSize/2;
      rightNode->splitIndex_j = user_IndexRoot->splitIndex + rightNodeSize/2;
      nodeLevelVec_Assign(rightNode->currLevel,rightNode);
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


void HODLR_Tree::correctIndices(){
  assert(rootNode != NULL);
  int offset_i,offset_j;
  offset_i = rootNode->min_i;
  offset_j = rootNode->min_j;
  rootNode->currLevel = 0;
  numLevels = 1;
  correctIndices(rootNode,offset_i,offset_j);
}

void HODLR_Tree::correctIndices(node* root,int offset_i,int offset_j){
  if (root->currLevel > (numLevels - 1))
    numLevels = root->currLevel + 1; 
  root->min_i -=  offset_i;
  root->min_j -=  offset_j;
  root->max_i -=  offset_i;
  root->max_j -=  offset_j;
  if (root->isLeaf == true)
    return;
  root->splitIndex_i -=  offset_i;
  root->splitIndex_j -=  offset_j;
  root->left->currLevel = root->currLevel + 1;
  root->right->currLevel = root->currLevel + 1;
  correctIndices(root->left,offset_i,offset_j);
  correctIndices(root->right,offset_i,offset_j);
}
  

void HODLR_Tree::set_sizeThreshold(const int input_sizeThreshold){
  if (rootNode == NULL)
    sizeThreshold = input_sizeThreshold;
  else{
    std::cout<<"Error! Cannot perform this operation at this time."<<std::endl;
    exit(EXIT_FAILURE);
  }
}

void HODLR_Tree::set_LRMethod(const std::string input_LRMethod){
  def_LRMethod = input_LRMethod;
  if (rootNode != NULL)
    set_LRMethod(rootNode,input_LRMethod);
  return;
}

void HODLR_Tree::set_LRMethod(node* root,const std::string input_LRMethod){
  root->LR_Method = input_LRMethod;
  if (root->isLeaf == true)
    return;
  set_LRMethod(root->left,input_LRMethod);
  set_LRMethod(root->right,input_LRMethod);
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
  
