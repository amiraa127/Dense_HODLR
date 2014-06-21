#include "HODLR_Matrix.hpp"
#include "helperFunctions.hpp"

int main(int argc,char* argv[]){
  bool ACA_TestConvergence = false;
  bool ACA_TestSpeed =false;
  bool boundaryLR_Test = false;
  switch (argc){
  case 1:
    {
      ACA_TestConvergence = true;
      ACA_TestSpeed = true;
      boundaryLR_Test = true;
      break;
    }
  case 2:
   {
     if (strcmp(argv[1],"-cACA") == 0)
       ACA_TestConvergence = true;
     else if (strcmp(argv[1],"-sACA") == 0) 
       ACA_TestSpeed = true;
     else if (strcmp(argv[1],"-boundary") == 0)
       boundaryLR_Test = true;
     else{
       std::cout<<"Error! Use -cACA flag for ACA convergence or -sACA for ACA speed analysis."<<std::endl;
       std::cout<<"Use -boundary for boundary low-rank approximation method analysis."<<std::endl;
       exit(EXIT_FAILURE);
     }
     break;
   }
  default:
    {
      std::cout<<"Error! invalid number of inputs. Only use single flags."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  /**************************************************Benchmarking ACA convegence on 1D interval***********************************************/  
  if (ACA_TestConvergence){
    int intervalSize = 10000;
    Eigen::VectorXd exactSoln = Eigen::VectorXd::LinSpaced(Eigen::Sequential,intervalSize,-2,2);
    
    std::cout<<"Benchmarking recLU solver convergence on 1D interval with ACA LR...."<<std::endl;
    
    std::cout<<"     Benchmarking quadratic kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/quadratiOBc_0",quadraticKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/quadratic_10",quadraticKernel,"recLU");
    
    std::cout<<"     Benchmarking multi quadratic kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/multiQuadratic_0",multiQuadraticKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/multiQuadratic_10",multiQuadraticKernel,"recLU");
    
    std::cout<<"     Benchmarking inverse quadratic kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/inverseQuadratic_0",inverseQuadraticKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/inverseQuadratic_10",inverseQuadraticKernel,"recLU");
    
    std::cout<<"     Benchmarking inverse multi quadratic kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/inverseMultiQuadratic_10",inverseMultiQuadraticKernel,"recLU");
    
    std::cout<<"     Benchmarking gaussian kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/gaussian_0",gaussianKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/gaussian_10",gaussianKernel,"recLU");
    
    std::cout<<"     Benchmarking exponential kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/exponential_0",exponentialKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/exponential_10",exponentialKernel,"recLU");
    
    std::cout<<"     Benchmarking logarithmic kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/logarithmic_0",logarithmicKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/logarithmic_10",logarithmicKernel,"recLU");
    
    std::cout<<"     Benchmarking one over r kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/oneOverR_0",oneOverRKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/oneOverR_10",oneOverRKernel,"recLU");
  
    std::cout<<"     Benchmarking one over r squared kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/oneOverSqR_0",oneOverSqRKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/oneOverSqR_10",oneOverSqRKernel,"recLU");
    
    std::cout<<"     Benchmarking log r kernel..."<<std::endl;
    testACASolverConv1DUnifromPts(-1,1,intervalSize,0,exactSoln,"ACA/recLU/convergence/logR_0",logRKernel,"recLU");
    testACASolverConv1DUnifromPts(-1,1,intervalSize,10,exactSoln,"ACA/recLU/convergence/logR_10",logRKernel,"recLU");
  }
  /**************************************************Benchmarking ACA solver Speed on 1D interval**************************************************/
  if (ACA_TestSpeed){
    
    std::cout<<"***********************************************************************"<<std::endl;
    std::cout<<"Benchmarking recLU solver speed on 1D interval with ACA LR...."<<std::endl;
      
    std::cout<<"     Benchmarking quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/quadratic_0",quadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/quadratic_10",quadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/quadratic_0",quadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/quadratic_10",quadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/quadratic_0",quadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/quadratic_10",quadraticKernel,"recLU");

    std::cout<<"     Benchmarking multi quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/multiQuadratic_0",multiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/multiQuadratic_10",multiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/multiQuadratic_0",multiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/multiQuadratic_10",multiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/multiQuadratic_0",multiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/multiQuadratic_10",multiQuadraticKernel,"recLU");

    std::cout<<"     Benchmarking inverse quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/inverseQuadratic_0",inverseQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/inverseQuadratic_10",inverseQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/inverseQuadratic_0",inverseQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/inverseQuadratic_10",inverseQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/inverseQuadratic_0",inverseQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/inverseQuadratic_10",inverseQuadraticKernel,"recLU");
    
    std::cout<<"     Benchmarking inverse multi quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/inverseMultiQuadratic_10",inverseMultiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/inverseMultiQuadratic_10",inverseMultiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/inverseMultiQuadratic_10",inverseMultiQuadraticKernel,"recLU");

    std::cout<<"     Benchmarking gaussian kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/gaussian_0",gaussianKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/gaussian_10",gaussianKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/gaussian_0",gaussianKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/gaussian_10",gaussianKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/gaussian_0",gaussianKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/gaussian_10",gaussianKernel,"recLU");

    std::cout<<"     Benchmarking exponential kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/exponential_0",exponentialKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/exponential_10",exponentialKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/exponential_0",exponentialKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/exponential_10",exponentialKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/exponential_0",exponentialKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/exponential_10",exponentialKernel,"recLU");
    
    std::cout<<"     Benchmarking logarithmic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/logarithmic_0",logarithmicKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/logarithmic_10",logarithmicKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/logarithmic_0",logarithmicKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/logarithmic_10",logarithmicKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/logarithmic_0",logarithmicKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/logarithmic_10",logarithmicKernel,"recLU");

    std::cout<<"     Benchmarking one over r kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/oneOverR_0",oneOverRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/oneOverR_10",oneOverRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/oneOverR_0",oneOverRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/oneOverR_10",oneOverRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/oneOverR_0",oneOverRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/oneOverR_10",oneOverRKernel,"recLU");
    
    std::cout<<"     Benchmarking one over r squared kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/oneOverSqR_0",oneOverSqRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/oneOverSqR_10",oneOverSqRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/oneOverSqR_0",oneOverSqRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/oneOverSqR_10",oneOverSqRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/oneOverSqR_0",oneOverSqRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/oneOverSqR_10",oneOverSqRKernel,"recLU");
    
    std::cout<<"     Benchmarking log r kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/recLU/speed/1e-7/logR_0",logRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/recLU/speed/1e-7/logR_10",logRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/recLU/speed/1e-5/logR_0",logRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/recLU/speed/1e-5/logR_10",logRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/recLU/speed/1e-10/logR_0",logRKernel,"recLU");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/recLU/speed/1e-10/logR_10",logRKernel,"recLU");
    
    std::cout<<"***********************************************************************"<<std::endl;
    std::cout<<"Benchmarking recLU solver speed on 1D interval with ACA LR for a fixed size matrix...."<<std::endl;

    std::cout<<"     Benchmarking quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/quadratic_0",quadraticKernel,10000,"recLU");

    std::cout<<"     Benchmarking multi quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/multiQuadratic_0",multiQuadraticKernel,10000,"recLU");

    std::cout<<"     Benchmarking inverse quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/inverseQuadratic_0",inverseQuadraticKernel,10000,"recLU");

    std::cout<<"     Benchmarking inverse multi quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,10000,"recLU");

    std::cout<<"     Benchmarking gaussian kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/gaussian_0",gaussianKernel,10000,"recLU");

    std::cout<<"     Benchmarking exponential kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/exponential_0",exponentialKernel,10000,"recLU");
    
    std::cout<<"     Benchmarking logarithmic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/logarithmic_0",logarithmicKernel,10000,"recLU");

    std::cout<<"     Benchmarking one over r kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/oneOverR_0",oneOverRKernel,10000,"recLU");
    
    std::cout<<"     Benchmarking one over r squared kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/oneOverSqR_0",oneOverSqRKernel,10000,"recLU");
    
    std::cout<<"     Benchmarking log r kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts_FixedSize(-1,1,0,1e-10,"ACA/recLU/speed/fixedSize/10k/1e-10/logR_0",logRKernel,10000,"recLU");
    

    std::cout<<"***********************************************************************"<<std::endl;
    std::cout<<"Benchmarking extendedSp solver speed on 1D interval with ACA LR...."<<std::endl;

    std::cout<<"     Benchmarking quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/quadratic_0",quadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/quadratic_10",quadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/quadratic_0",quadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/quadratic_10",quadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/quadratic_0",quadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/quadratic_10",quadraticKernel,"extendedSp");
   
    std::cout<<"     Benchmarking multi quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/multiQuadratic_0",multiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/multiQuadratic_10",multiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/multiQuadratic_0",multiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/multiQuadratic_10",multiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/multiQuadratic_0",multiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/multiQuadratic_10",multiQuadraticKernel,"extendedSp");

    std::cout<<"     Benchmarking inverse quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/inverseQuadratic_0",inverseQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/inverseQuadratic_10",inverseQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/inverseQuadratic_0",inverseQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/inverseQuadratic_10",inverseQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/inverseQuadratic_0",inverseQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/inverseQuadratic_10",inverseQuadraticKernel,"extendedSp");

    std::cout<<"     Benchmarking inverse multi quadratic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/inverseMultiQuadratic_10",inverseMultiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/inverseMultiQuadratic_10",inverseMultiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/inverseMultiQuadratic_0",inverseMultiQuadraticKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/inverseMultiQuadratic_10",inverseMultiQuadraticKernel,"extendedSp");

    std::cout<<"     Benchmarking gaussian kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/gaussian_0",gaussianKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/gaussian_10",gaussianKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/gaussian_0",gaussianKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/gaussian_10",gaussianKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/gaussian_0",gaussianKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/gaussian_10",gaussianKernel,"extendedSp");

    std::cout<<"     Benchmarking exponential kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/exponential_0",exponentialKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/exponential_10",exponentialKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/exponential_0",exponentialKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/exponential_10",exponentialKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/exponential_0",exponentialKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/exponential_10",exponentialKernel,"extendedSp");
    
    std::cout<<"     Benchmarking logarithmic kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/logarithmic_0",logarithmicKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/logarithmic_10",logarithmicKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/logarithmic_0",logarithmicKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/logarithmic_10",logarithmicKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/logarithmic_0",logarithmicKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/logarithmic_10",logarithmicKernel,"extendedSp");

    std::cout<<"     Benchmarking one over r kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/oneOverR_0",oneOverRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/oneOverR_10",oneOverRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/oneOverR_0",oneOverRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/oneOverR_10",oneOverRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/oneOverR_0",oneOverRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/oneOverR_10",oneOverRKernel,"extendedSp");
    
    std::cout<<"     Benchmarking one over r squared kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/oneOverSqR_0",oneOverSqRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/oneOverSqR_10",oneOverSqRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/oneOverSqR_0",oneOverSqRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/oneOverSqR_10",oneOverSqRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/oneOverSqR_0",oneOverSqRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/oneOverSqR_10",oneOverSqRKernel,"extendedSp");
    
    std::cout<<"     Benchmarking log r kernel..."<<std::endl;
    testACASolverSpeed1DUniformPts(-1,1,0,1e-7,"ACA/extendedSp/speed/1e-7/logR_0",logRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-7,"ACA/extendedSp/speed/1e-7/logR_10",logRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-5,"ACA/extendedSp/speed/1e-5/logR_0",logRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-5,"ACA/extendedSp/speed/1e-5/logR_10",logRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,0,1e-10,"ACA/extendedSp/speed/1e-10/logR_0",logRKernel,"extendedSp");
    testACASolverSpeed1DUniformPts(-1,1,10,1e-10,"ACA/extendedSp/speed/1e-10/logR_10",logRKernel,"extendedSp");    
  }
  if (boundaryLR_Test){
    std::cout<<"Testing solver for structured meshes.."<<std::endl;
    testBoundaryLR("boundaryLR/FETI/structured/400k/input/400_front_num_5_level_0","boundaryLR/FETI/structured/400k/input/400_front_num_5_level_0_Graph","boundaryLR/FETI/structured/400k/results/400_front_num_5_level_0_",1e-1,2);   
    testBoundaryLR("boundaryLR/FETI/structured/400k/input/400_front_num_4_level_1","boundaryLR/FETI/structured/400k/input/400_front_num_4_level_1_Graph","boundaryLR/FETI/structured/400k/results/400_front_num_4_level_1_",1e-1,2);   
    testBoundaryLR("boundaryLR/FETI/structured/400k/input/400_front_num_2_level_2","boundaryLR/FETI/structured/400k/input/400_front_num_2_level_2_Graph","boundaryLR/FETI/structured/400k/results/400_front_num_2_level_2_",1e-1,2);   
    std::cout<<"Testing solver for unstructured meshes.."<<std::endl;
    testBoundaryLR("boundaryLR/stiffness/unstructured/300k/input/300_front_num_6_level_0","boundaryLR/stiffness/unstructured/300k/input/300_front_num_6_level_0_Graph","boundaryLR/stiffness/unstructured/300k/results/300_front_num_6_level_0_",1e-1,2);
    testBoundaryLR("boundaryLR/stiffness/unstructured/300k/input/300_front_num_4_level_1","boundaryLR/stiffness/unstructured/300k/input/300_front_num_4_level_1_Graph","boundaryLR/stiffness/unstructured/300k/results/300_front_num_4_level_1_",1e-1,2);
    testBoundaryLR("boundaryLR/stiffness/unstructured/300k/input/300_front_num_2_level_2","boundaryLR/stiffness/unstructured/300k/input/300_front_num_2_level_2_Graph","boundaryLR/stiffness/unstructured/300k/results/300_front_num_2_level_2_",1e-1,2);
    
  }
  return 0;
}
