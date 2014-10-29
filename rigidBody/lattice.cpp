//C version of matlab file for sphere triangulation

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <fstream>
#include <iostream>
#include "MobFunctions.h"


struct kernelDataInput{
  double DX;
  double* X;
};

double RPY_Kernel(int i,int j,void* inputData){
  return get_RPY_Matrix_Entry(((kernelDataInput*)inputData)->X,((kernelDataInput*)inputData)->DX, i, j);
};

void print_help()
{
  printf("Usage for 3D only: vertex_filename \t length_between_objects \t N_x (size of the lattice in x direction) \t N_y \t N_z \n");
  printf("Default centering is the origin (0,0,0);\n");
  exit(1);
}

//create coordinates of point for 3D cube
int main(int arc, char**argv)
{
  
  if (arc<6) print_help();

  double center[3];
  //change this in case you need different centering
  for (int idir=0;idir<3;++idir)  center[idir]=0.0;

  char* filename= argv[1];
  std::ifstream myfile (filename);
  if (!myfile.is_open()) 
  {
    std::cout<<"Wrong input file name"<<std::endl;
    exit(1);
  }
  int pointsPerSphere; //number of blobs
  std::string line;
  double DX;//grid resoultion step
  
  getline(myfile,line);
  std::istringstream n_line(line);
  n_line >> pointsPerSphere;
  n_line >> DX;

  //allocate space for positions of all blobs
  double *X = new double[3 * pointsPerSphere];

  for (int row=0;row < pointsPerSphere;row++){
    if (!std::getline(myfile,line)) {
      std::cout <<"Error in input  file particle.vertex"<<std::endl;
      exit(1);
    } 
    std::istringstream iss(line);
    for (int rdir=0;rdir<3;rdir++){
	iss>>X[row*3+rdir];
      }
  }
  myfile.close();
  
  double length = strtod(argv[2],NULL);
  int size[3];
  size[0] = strtol(argv[3],NULL,10); // Number of markers x
  size[1] = strtol(argv[4],NULL,10); // Number of markers y
  size[2] = strtol(argv[5],NULL,10); // Number of markers z


  std::cout<<"Total number of vertices is "<<pointsPerSphere * size[0] * size[1] * size[2]<<std::endl;
  int matrixSize  = pointsPerSphere * size[0] * size[1] * size[2] * 3;
  double* globalX = new double[matrixSize];
  int numSpheres  = size[0] * size[1] * size[2];
  int diagBlkSize = pointsPerSphere * 3;
  int sphereIdx   = 0;
  Eigen::MatrixXd sphereCenters(numSpheres,NDIM + 1);
  for(int ipart=0;ipart<size[0];ipart++) 
    for(int jpart=0;jpart<size[1];jpart++) 
      for(int kpart=0;kpart<size[2];kpart++) 
	{
	  double offset[3];
	  
	  offset[0]=center[0]+((double)ipart-size[0]/2.0)*length;
	  offset[1]=center[1]+((double)jpart-size[1]/2.0)*length;
	  offset[2]=center[2]+((double)kpart-size[2]/2.0)*length;
	  for (int rdir = 0; rdir < NDIM; rdir++)
	    sphereCenters(sphereIdx,rdir) = offset[rdir];
	  sphereCenters(sphereIdx,NDIM) = sphereIdx;
	  sphereIdx ++;
	}
 
 
  kernelDataInput RPY_Matrix_Input;
  RPY_Matrix_Input.X =  globalX;
  RPY_Matrix_Input.DX = DX;
  
  int numPoints = matrixSize/3;

  //Sorting according to HODLR.

  user_IndexTree usrTree = get_KDTree_Sorted(sphereCenters,pointsPerSphere,1000);
  
  int globalX_Idx = 0;
  
  for (int i = 0; i < numSpheres; i++){
    for (int row = 0;row < pointsPerSphere;row++)
      for (int rdir = 0;rdir < 3;rdir++){
	globalX[globalX_Idx] = X[row * 3 + rdir] + sphereCenters(i,rdir);
	globalX_Idx++;
      } 
  }

  
  Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,matrixSize,-2,2);
  HODLR_Matrix kernelHODLR(matrixSize,matrixSize,RPY_Kernel,&RPY_Matrix_Input,diagBlkSize + 1,usrTree); 
  kernelMatrix exactMatrixKernel(matrixSize,matrixSize,RPY_Kernel,&RPY_Matrix_Input);
  Eigen::VectorXd inputF = exactMatrixKernel * exactSoln; // Computing the RHS 
  kernelHODLR.printResultInfo = true; // Comment this line out if you don't want the l2 error
  std::cout<<"Solving"<<std::endl;
  kernelHODLR.set_LRTolerance(1e-4); // Low-rank approximation tolerance
  kernelHODLR.set_LeafConst();
  Eigen::VectorXd solverSoln = kernelHODLR.recLU_Solve(inputF);
  Eigen::VectorXd difference = solverSoln - exactSoln;
  double relError = difference.norm()/exactSoln.norm();
  std::cout<<relError<<std::endl;
  

  delete[] X;
  delete[] globalX;
} 
