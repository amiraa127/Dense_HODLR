#include "MobFunctions.h"
#include "HODLR_Matrix.hpp"
#include "kernel.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <ctime>

struct kernelDataInput{
  double DX;
  double* X;
};

double RPY_Kernel(int i,int j,void* inputData){
  return get_RPY_Matrix_Entry(((kernelDataInput*)inputData)->X,((kernelDataInput*)inputData)->DX, i, j);
};

Eigen::SparseMatrix<double> createGraph(double* X, int numPoints,double distanceThreshold){
  std::vector<Eigen::Triplet<double,int> > tripletVec;
  for (int i = 0; i < numPoints; i++)
    for (int j = 0; j < numPoints; j++){
      double dist = 0;
      for (int dim = 0; dim < NDIM; dim++)
	dist += pow(X[i*NDIM + dim] - X[j*NDIM + dim],2);
      dist = sqrt(dist);
      if (dist <= distanceThreshold)
	for (int dim = 0; dim < NDIM; dim++){
	  Eigen::Triplet<double,int> currEntry(i*NDIM + dim,j*NDIM + dim,1);
	  tripletVec.push_back(currEntry);
	}
    }
  Eigen::SparseMatrix<double> result(numPoints*NDIM,numPoints*NDIM);
  result.setFromTriplets(tripletVec.begin(),tripletVec.end());
  return result;
}

int main (int argc, char* argv[]){
  
    //First read the input file to save all coordinates to array X
    char* filename= argv[1];
    std::ifstream myfile (filename);
    if (!myfile.is_open()) 
    {
      std::cout<<"Error! Cannot open file."<<std::endl;
      exit(EXIT_FAILURE);
    }

    int numPoints; //number of blobs
    double DX; //grid resolution step
    std::string line, line_u;
    
    getline(myfile,line);
    std::istringstream n_line(line);
    n_line >> numPoints;
    n_line >> DX;
    const int matrixSize = NDIM * numPoints; 
    
    //allocate space for positions of all blobs
    double *X = new double[matrixSize];

    for (int row = 0; row < numPoints; row++)
    {
	if (!std::getline(myfile,line)) 
	{
	    std::cout <<"Error in input file particle.vertex"<<std::endl;
	    exit(EXIT_FAILURE);
	}
	std::istringstream iss(line);
	std::istringstream iss_u(line_u);
	for (int rdir = 0;rdir < NDIM; rdir++) 
	{
	    iss>>X[row * NDIM + rdir];
	}
    }
    myfile.close();

    //Sorting according to HODLR.
    Eigen::MatrixXd XX(numPoints, NDIM+1);
    
    for (int i = 0; i < matrixSize; ++i){
	XX(i/NDIM,i%NDIM) = X[i];
	XX(i/NDIM,NDIM)   = i/NDIM;
    }
    get_KDTree_Sorted(XX,0);
    
    int *permut=new int[matrixSize];
    for (int i = 0; i < matrixSize; ++i) 
      permut[((int)XX(i/NDIM,NDIM))*NDIM+i%NDIM] = i;
    
    //Reorder X according to HODLR
    for (int i = 0; i < matrixSize; ++i){
	X[i]=XX(i/NDIM,i%NDIM);
    }

    kernelDataInput RPY_Matrix_Input;
    RPY_Matrix_Input.X =  X;
    RPY_Matrix_Input.DX = DX;


    Eigen::SparseMatrix<double> interactionGraph = createGraph(X,numPoints,5*DX);

    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
    //HODLR_Matrix kernelHODLR(numPoints,numPoints,RPY_Kernel,&RPY_Matrix_Input,50);
    HODLR_Matrix kernelHODLR(matrixSize,matrixSize,RPY_Kernel,&RPY_Matrix_Input,interactionGraph,100); 
    kernelMatrix exactMatrixKernel(matrixSize,matrixSize,RPY_Kernel,&RPY_Matrix_Input);
    Eigen::VectorXd inputF = exactMatrixKernel * exactSoln;
    kernelHODLR.printResultInfo = true;
    std::cout<<"Solving"<<std::endl;
    //Eigen::VectorXd solverSoln = kernelHODLR.recLU_Solve(inputF);
    kernelHODLR.set_BoundaryDepth(4);
    Eigen::VectorXd solverSoln = kernelHODLR.iterative_Solve(inputF,10,1e-10,1e-5,"PS_Boundary","recLU");
    Eigen::VectorXd difference = solverSoln - exactSoln;
    double relError = difference.norm()/exactSoln.norm();
    std::cout<<relError<<std::endl;



    
    delete[] X;
    delete[] permut;
}
