#include <math.h>
#include <Eigen/Dense>
#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

#define isRPY 1
#define NDIM 3 //Dimension 
#define MU 1.0 //Fluid viscosity
#define HYD_RAD 1.47//Hydrodynamic Raduis for standard IB6 delta function

#ifndef KRON
#define KRON(i, j) ((i==j)?1:0) //Kronecker symbol
#endif

static double
get_sqnorm(
    const double *a_vec)
{
#if (NDIM==3)
  return a_vec[0]*a_vec[0]+a_vec[1]*a_vec[1]+a_vec[2]*a_vec[2];
#elif(NDIM==2)
  return a_vec[0]*a_vec[0]+a_vec[1]*a_vec[1];
#endif
}

static double 
F_R_INF(const double r)
/*!
  This functions returns the value of normilized f(r) function for empirical fit
*/
{
    double temp_f_r= (1.0/(3.0*HYD_RAD/4.0+r));
    const double numer  = 0.336317*r+0.0023172*r*r*r;
    double denom;
    if (r < 1.0)  denom = 0.132882*r*r*r*r+0.132882*r*r*r+0.0911596*r*r+1.0;
    else denom = 1.0+0.0911596*r*r+0.132882*r*r*r+0.0419235*r*r*r*r;
    temp_f_r += (numer/denom);

    return temp_f_r;
}// _F_R_INF

static double 
G_R_INF(const double r)
/*!
  This functions returns the value of normilized g(r)function for empirical fit
*/
{
    const double numer = r*r;
    double denom;
    if (r < 1.0 ) denom = r*r*r + 0.253716*r*r + 17.4211;
    else denom =  17.4211 + 0.253716*r*r + r*r*r;
    return (numer/denom);
}// _G_R_INF


void
get_RPY_SubMatrix(const double *X,
		  const double DX,
		  const int i,
		  const int j,
		  double* subM)
/*!
  This function function generates and stores the the RPY submatrix for i,j blobs into array subM
  X must be initialized before calling this function,
  subM must have a size non less than sizeof(double)*NDIM*NDIM.
*/
{
    const double R_H = HYD_RAD*DX;
    const double mu_tt = 1./(6.0*M_PI*MU*R_H);
    double r_vec[NDIM];
    if (i==j)
    {
	for (int idir = 0; idir < NDIM; idir++) 
	    for (int jdir = 0; jdir < NDIM; jdir++) 
		subM[idir*NDIM+jdir] = mu_tt*KRON(idir,jdir);
    }
    else
    {
	for (int cdir = 0; cdir < NDIM; cdir++) 
	{
	    r_vec[cdir] = X[i*NDIM+cdir] - X[j*NDIM+cdir]; //r(i) - r(j)
	}
	  
	const double rsq = get_sqnorm(r_vec); 
	const double r   = sqrt(rsq);
	for (int idir = 0; idir < NDIM; idir++) 
	    for (int jdir = idir; jdir < NDIM; jdir++) 
	    {
		if (r<=2.0*R_H)
		{
		    subM[idir*NDIM+jdir] = mu_tt*(1-9.0/32.0*r/R_H)*KRON(idir,jdir)
			+mu_tt*r_vec[idir]*r_vec[jdir]/rsq*3.0*r/32./R_H;
		    if (idir != jdir)
			subM[jdir*NDIM+idir] = mu_tt*(1-9.0/32.0*r/R_H)*KRON(idir,jdir)
			    +mu_tt*r_vec[idir]*r_vec[jdir]/rsq*3.0*r/32./R_H;

		}
		else
		{
		    double cube=R_H*R_H*R_H/r/r/r;
		    subM[idir*NDIM+jdir] = mu_tt*(3.0/4.0*R_H/r+1.0/2.0*cube)*KRON(idir,jdir)
			+ mu_tt*r_vec[idir]*r_vec[jdir]/rsq*(3.0/4.0*R_H/r-3.0/2.0*cube);
		    if (idir != jdir)
			subM[jdir*NDIM+idir] = mu_tt*(3.0/4.0*R_H/r+1.0/2.0*cube)*KRON(idir,jdir)
			    + mu_tt*r_vec[idir]*r_vec[jdir]/rsq*(3.0/4.0*R_H/r-3.0/2.0*cube);
		}
	    }
    }
}

void
get_FIT_SubMatrix(const double *X,
		  const double DX,
		  const int i,
		  const int j,
		  double* subM)
/*!
  This function function generates and stores the the FIT submatrix for i,j blobs into array subM
  X must be initialized before calling this function,
  subM must have a size non less than sizeof(double)*NDIM*NDIM.
*/
{
    double mu_tt = 1./(8.0*M_PI*MU*DX);
    double r_vec[NDIM];
    for (int cdir = 0; cdir < NDIM; cdir++) 
    {
	r_vec[cdir] = X[i*NDIM+cdir] - X[j*NDIM+cdir]; //r(i) - r(j)
    }
    
    const double rsq = get_sqnorm(r_vec); 
    const double r   = sqrt(rsq);
	for (int idir = 0; idir < NDIM; idir++) 
	    for (int jdir = 0; jdir < NDIM; jdir++) 
	    {
		subM[idir*NDIM+jdir] = mu_tt*F_R_INF(r/DX)*KRON(idir,jdir);

		if (i!=j) 
		{
		    subM[idir*NDIM+jdir] += mu_tt*r_vec[idir]*r_vec[jdir]/rsq*G_R_INF(r/DX);
		}

	    }
}

void
get_MatrixVectorProduct(const double *X,
			const double DX,
			const int N,
			const double *F,
			double* U)
/*!
  This function function computes a product of RPY or FIT mobility matrix to vector F 
  X must be initialized before calling this function,
  U must have a size non less than sizeof(double)*N*NDIM.
*/
{
    //initilize U
    for (int ipart=0; ipart<NDIM*N;++ipart) U[ipart]=0.;

    double *SubMM=new double[NDIM*NDIM];

    for (int ipart=0; ipart<N;++ipart)
    	for (int jpart=ipart; jpart<N;++jpart)
	{
	    if (isRPY) get_RPY_SubMatrix(X,DX,ipart,jpart, SubMM);
	    else get_FIT_SubMatrix(X,DX,ipart,jpart, SubMM);

	    //compute part of product for i,j submatrix
	    for (int idir=0; idir<NDIM;++idir)
		for (int jdir=0; jdir<NDIM;++jdir)
		{
		    U[ipart*NDIM+idir] +=SubMM[idir*NDIM+jdir]*F[jpart*NDIM+jdir];
		    if (ipart!=jpart)
			U[jpart*NDIM+idir] +=SubMM[idir*NDIM+jdir]*F[ipart*NDIM+jdir];
		}
	}
}


double
get_RPY_Matrix_Entry(const double *X,
		     const double DX,
		     const int i,
		     const int j)
/*!
  This function function generates and returns the (i,j)th element of the RPY matrix based on coordinates of points (array X).
  X must be initialized before calling this function.
*/
{
    const double R_H = HYD_RAD*DX;
    const double mu_tt = 1./(6.0*M_PI*MU*R_H);
    double RPY_entry;
    const int row=i/NDIM;
    const int col=j/NDIM;
    
    double r_vec[NDIM];
    if (row==col)
    {
	const int idir=i%NDIM;
	const int jdir=j%NDIM;
	RPY_entry = mu_tt*KRON(idir,jdir);
	  
    }
    else
    {
	for (int cdir = 0; cdir < NDIM; cdir++) 
	{
	    r_vec[cdir] = X[row*NDIM+cdir] - X[col*NDIM+cdir]; //r(i) - r(j)
	}
	  
	const double rsq = get_sqnorm(r_vec); 
	const double r   = sqrt(rsq);
	const int idir=i%NDIM;
	const int jdir=j%NDIM;
	if (r<=2.0*R_H)
	{
	    RPY_entry = mu_tt*(1-9.0/32.0*r/R_H)*KRON(idir,jdir)
		+mu_tt*r_vec[idir]*r_vec[jdir]/rsq*3.0*r/32./R_H;
	}
	else
	{
	    double cube=R_H*R_H*R_H/r/r/r;
	    RPY_entry = mu_tt*(3.0/4.0*R_H/r+1.0/2.0*cube)*KRON(idir,jdir)
		+ mu_tt*r_vec[idir]*r_vec[jdir]/rsq*(3.0/4.0*R_H/r-3.0/2.0*cube);
	}
    }
    return RPY_entry;
}

double
get_FIT_Matrix_Entry(const double *X,
		     const double DX,
		     const int i,
		     const int j)
/*!
  This function function generates and returns the (i,j)th element of the empirical FIT matrix based on coordinates of points (array X).
  X must be initialized before calling this function.
*/
{
    double mu_tt = 1./(8.0*M_PI*MU*DX);
   
    double FIT_entry=0.;
    const int row=i/NDIM;
    const int col=j/NDIM;
    double r_vec[NDIM];

    for (int cdir = 0; cdir < NDIM; cdir++) 
    {
	r_vec[cdir] = X[row*NDIM+cdir] - X[col*NDIM+cdir]; //r(i) - r(j)
    }
    
    const double rsq = get_sqnorm(r_vec); 
    const double r   = sqrt(rsq);
    const int idir=i%NDIM;
    const int jdir=j%NDIM;

		  
    FIT_entry = mu_tt*F_R_INF(r/DX)*KRON(idir,jdir);
    if (row != col) FIT_entry += mu_tt*r_vec[idir]*r_vec[jdir]/rsq*G_R_INF(r/DX);

    return FIT_entry;
} 
void 
mergeSortedLists(
    Eigen::MatrixXd& list1, 
    Eigen::MatrixXd& list2, 
    unsigned index, 
    Eigen::MatrixXd& finalList) 
{
    unsigned N1	=	list1.rows();
    unsigned N2	=	list2.rows();
    unsigned j1	=	0;
    unsigned j2	=	0;
    unsigned j	=	0;
    while (j1<N1 && j2 <N2) {
	if (list1(j1,index) < list2(j2,index)) {
	    finalList.row(j)	=	list1.row(j1);
	    ++j1;
	}
	else {
	    finalList.row(j)	=	list2.row(j2);
	    ++j2;
	}
	++j;
    }
    while (j1<N1) {
	finalList.row(j)	=	list1.row(j1);
	++j1;
	++j;
    }
    while (j2<N2) {
	finalList.row(j)	=	list2.row(j2);
	++j2;
	++j;
    }
}

void 
mergeSort(
    Eigen::MatrixXd& locations, 
    unsigned index) 
{
    unsigned N		=	locations.rows();
    if (N==1) {
	return;
    }
    else {
	///	Number of points in the left cluster.
	unsigned Nleft		=	N/2;
	
	///	Number of points in the right cluster.
	unsigned Nright		=	N-Nleft;
	
	///	Dimension of the space.
	unsigned nDimensions	=	locations.cols();
	
	///	Left locations.
	Eigen::MatrixXd leftLocations	=	locations.block(0,0,Nleft,nDimensions);
	
	///	Right locations.
	Eigen::MatrixXd rightLocations	=	locations.block(Nleft,0,Nright,nDimensions);
	
	///	Mergesort for the left.
	mergeSort(leftLocations, index);
	
	///	Mergesort for the right.
	mergeSort(rightLocations, index);
	
	///	Merge the sorted left and right lists.
	mergeSortedLists(leftLocations, rightLocations, index, locations);
    }
}

void 
get_KDTree_Sorted(
    Eigen::MatrixXd& locations, 
    unsigned index) 
{
    ///	Get the total number of points.
    unsigned N		=	locations.rows();
    if (N==1) {
	return;
    }
    else {
	///	Number of points in the left cluster.
	unsigned Nleft		=	N/2;
	
	///	Number of points in the right cluster.
	unsigned Nright		=	N-Nleft;
	
	///	Dimension of the space.
	unsigned nDimensions	=	locations.cols();
	
	///	Merge sort on the input locations based on the coordinate index%2.
	mergeSort(locations, index%nDimensions);
	
	///	Obtain the left and right locations.
	Eigen::MatrixXd leftLocations	=	locations.block(0,0,Nleft,nDimensions);
	Eigen::MatrixXd rightLocations	=	locations.block(Nleft,0,Nright,nDimensions);
	
	///	Sort the left and right locations based on a KDTree.
	get_KDTree_Sorted(leftLocations, index+1);
	get_KDTree_Sorted(rightLocations, index+1);
	
	///	Output the locations.
	locations.block(0,0,Nleft,nDimensions)		=	leftLocations;
	locations.block(Nleft,0,Nright,nDimensions)	=	rightLocations;
    }
}

void get_KDTree_Sorted(Eigen::MatrixXd& locations,unsigned index,int globalStartIdx,user_IndexTree::node* usrTreeNode,int pointsPerSphere,int rankUpperBound,std::string LR_Method){
 
  ///	Get the total number of points.  
  usrTreeNode->topOffDiag_minRank  = -1;
  usrTreeNode->bottOffDiag_minRank = -1; 
  usrTreeNode->topOffDiag_maxRank  = rankUpperBound;
  usrTreeNode->bottOffDiag_maxRank = rankUpperBound;

  usrTreeNode->LR_Method = LR_Method;
  
  unsigned N = locations.rows();
  if (N == 1) {
    usrTreeNode->left = NULL;
    usrTreeNode->right = NULL;
    usrTreeNode->splitIndex = -1;
    return;
  }
  
  ///	Number of points in the left cluster.
  unsigned Nleft  = N/2;

  ///	Number of points in the right cluster.
  unsigned Nright = N - Nleft;

  ///   Detremine splitIndex
  int splitIdx = NDIM * pointsPerSphere *(Nleft + globalStartIdx) - 1;
  usrTreeNode->splitIndex = splitIdx;
  
  ///	Dimension of the space.
  unsigned nDimensions = locations.cols();
  
  ///	Merge sort on the input locations based on the coordinate index%2.
  mergeSort(locations, index%nDimensions);
  
  user_IndexTree::node* left  = new user_IndexTree::node;
  user_IndexTree::node* right = new user_IndexTree::node;
  usrTreeNode->left  = left;
  usrTreeNode->right = right;
  
  ///	Obtain the left and right locations.
  Eigen::MatrixXd leftLocations  = locations.block(0,0,Nleft,nDimensions);
  Eigen::MatrixXd rightLocations = locations.block(Nleft,0,Nright,nDimensions);
    
  ///	Sort the left and right locations based on a KDTree.
  get_KDTree_Sorted(leftLocations, index+1,globalStartIdx,usrTreeNode->left,pointsPerSphere,rankUpperBound,LR_Method);
  get_KDTree_Sorted(rightLocations, index+1,globalStartIdx + Nleft,usrTreeNode->right,pointsPerSphere,rankUpperBound,LR_Method);
  
  ///	Output the locations.
  locations.block(0,0,Nleft,nDimensions) = leftLocations;
  locations.block(Nleft,0,Nright,nDimensions)  = rightLocations;

}

user_IndexTree get_KDTree_Sorted(Eigen::MatrixXd& locations,int pointsPerSphere,int rankUpperBound,std::string LR_Method = "partialPiv_ACA"){
  user_IndexTree result;
  result.rootNode = new user_IndexTree::node;
  int numPoints = locations.rows();
  get_KDTree_Sorted(locations,0,0,result.rootNode,pointsPerSphere,rankUpperBound,LR_Method);
  return result;
}


class hodlr_mobility : public HODLR_Matrix {
public:
    hodlr_mobility( 
	double *inX,
	const double DX)
	:d_DX(DX)
    {
	X = inX;
    }
  
  double get_Matrix_Entry(const unsigned i, const unsigned j) 
    {
	//check if RPY or FIT Mobility matrix to generate
	if (isRPY) return get_RPY_Matrix_Entry(X, d_DX, i, j);
	else return get_FIT_Matrix_Entry(X, d_DX, i, j);	
    }
  
private:
    double  d_DX; 
    double  *X;
};

