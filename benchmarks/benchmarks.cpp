#include "HODLR_Matrix.hpp"
#include "helperFunctions.hpp"

int main(int argc,char* argv[]){
  bool ACA_TestConvergence = false;
  bool ACA_TestSpeed       = false;
  bool boundaryLR_Test     = false;
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
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"Testing solver for 3D FETI structured meshes.."<<std::endl;
    std::cout<<"Timing tests........."<<std::endl;
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_5_level_0"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/structured/400k/input/400_front_num_5_level_0","boundaryLR/FETI/structured/400k/input/400_front_num_5_level_0_Graph","boundaryLR/FETI/structured/400k/results_timing/400_front_num_5_level_0_",1e-1,30,1);   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_4_level_1"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/structured/400k/input/400_front_num_4_level_1","boundaryLR/FETI/structured/400k/input/400_front_num_4_level_1_Graph","boundaryLR/FETI/structured/400k/results_timing/400_front_num_4_level_1_",1e-1,30,1); 
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_3_level_1"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/structured/400k/input/400_front_num_3_level_1","boundaryLR/FETI/structured/400k/input/400_front_num_3_level_1_Graph","boundaryLR/FETI/structured/400k/results_timing/400_front_num_3_level_1_",1e-1,30,1);   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_2_level_2"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/structured/400k/input/400_front_num_2_level_2","boundaryLR/FETI/structured/400k/input/400_front_num_2_level_2_Graph","boundaryLR/FETI/structured/400k/results_timing/400_front_num_2_level_2_",1e-1,30,1);
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_1_level_2"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/structured/400k/input/400_front_num_1_level_2","boundaryLR/FETI/structured/400k/input/400_front_num_1_level_2_Graph","boundaryLR/FETI/structured/400k/results_timing/400_front_num_1_level_2_",1e-1,100,1);   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_0_level_2"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/structured/400k/input/400_front_num_0_level_2","boundaryLR/FETI/structured/400k/input/400_front_num_0_level_2_Graph","boundaryLR/FETI/structured/400k/results_timing/400_front_num_0_level_2_",1e-1,100,1);
    
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"Testing solver for 3D FETI unstructured meshes.."<<std::endl;
    std::cout<<"Timing tests........."<<std::endl;
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_0_level_13"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/unstructured/400k/input/400_front_num_0_level_13","boundaryLR/FETI/unstructured/400k/input/400_front_num_0_level_13_Graph","boundaryLR/FETI/unstructured/400k/results_timing/400_front_num_0_level_13_",1e-1,30,1);   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_1_level_11"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/unstructured/400k/input/400_front_num_1_level_11","boundaryLR/FETI/unstructured/400k/input/400_front_num_1_level_11_Graph","boundaryLR/FETI/unstructured/400k/results_timing/400_front_num_1_level_11_",1e-1,30,1); 
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_2_level_7"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/unstructured/400k/input/400_front_num_2_level_7","boundaryLR/FETI/unstructured/400k/input/400_front_num_2_level_7_Graph","boundaryLR/FETI/unstructured/400k/results_timing/400_front_num_2_level_7_",1e-1,30,1);   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_3_level_6"<<std::endl;
    testBoundaryLRSolver("boundaryLR/FETI/unstructured/400k/input/400_front_num_3_level_6","boundaryLR/FETI/unstructured/400k/input/400_front_num_3_level_6_Graph","boundaryLR/FETI/unstructured/400k/results_timing/400_front_num_3_level_6_",1e-1,30,1); 
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"Testing solver for 3D unstructured meshes.."<<std::endl;
    std::cout<<"Timing tests........."<<std::endl;
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_6_level_0"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_6_level_0","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_6_level_0_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_timing/300_front_num_6_level_0_",1e-1,30,1);
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_5_level_1"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_5_level_1","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_5_level_1_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_timing/300_front_num_5_level_1_",1e-1,30,1); 
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_4_level_1"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_4_level_1","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_4_level_1_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_timing/300_front_num_4_level_1_",1e-1,30,1);
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_3_level_2"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_3_level_2","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_3_level_2_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_timing/300_front_num_3_level_2_",1e-1,30,1);
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_2_level_2"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_2_level_2","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_2_level_2_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_timing/300_front_num_2_level_2_",1e-1,30,1);    
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_1_level_2"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_1_level_2","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_1_level_2_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_timing/300_front_num_1_level_2_",1e-1,30,1);
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_0_level_2"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_0_level_2","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_0_level_2_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_timing/300_front_num_0_level_2_",1e-1,30,1);    
    

    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"*******************************************************"<<std::endl;
    // Create usrTrees fro turbine blade results
    std::cout<<"Testing solver for 2D unstructured meshes(turbine blade)"<<std::endl;
    
    // Create 9k user tree
    // Create custom indexing tree
    
    user_IndexTree usrTree;
    usrTree.rootNode = new user_IndexTree::node;
    usrTree.rootNode->splitIndex = 1498;
    usrTree.rootNode->topOffDiag_minRank = -1;
    usrTree.rootNode->bottOffDiag_minRank = -1;
    usrTree.rootNode->LR_Method = "PS_Boundary";

    user_IndexTree::node* leftNode = new user_IndexTree::node;
    user_IndexTree::node* rightNode = new user_IndexTree::node;
    usrTree.rootNode->left = leftNode;
    usrTree.rootNode->right = rightNode;

    rightNode->splitIndex = 6497;
    rightNode->topOffDiag_minRank = -1;
    rightNode->bottOffDiag_minRank = -1;
    rightNode->LR_Method = "PS_Boundary";

    leftNode->splitIndex = 749;
    leftNode->topOffDiag_minRank = -1;
    leftNode->bottOffDiag_minRank = -1;
    leftNode->LR_Method = "PS_Boundary";

    usrTree.setChildren_NULL(leftNode);
    usrTree.setChildren_NULL(rightNode);
   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"blade_9k_level_0"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/2D/input/blade_9k_level_0","boundaryLR/stiffness/unstructured/2D/input/blade_9k_level_0_Graph","boundaryLR/stiffness/unstructured/2D/results_timing/blade_9k_level_0_",1e-1,30,2,usrTree);
  
    // Create 12k userTree
    usrTree.rootNode = new user_IndexTree::node;
    usrTree.rootNode->splitIndex = 2098;
    usrTree.rootNode->topOffDiag_minRank = -1;
    usrTree.rootNode->bottOffDiag_minRank = -1;
    usrTree.rootNode->LR_Method = "PS_Boundary";

    leftNode  = new user_IndexTree::node;
    rightNode = new user_IndexTree::node;
    usrTree.rootNode->left = leftNode;
    usrTree.rootNode->right = rightNode;

    rightNode->splitIndex = 9097;
    rightNode->topOffDiag_minRank = -1;
    rightNode->bottOffDiag_minRank = -1;
    rightNode->LR_Method = "PS_Boundary";

    leftNode->splitIndex = 1048;
    leftNode->topOffDiag_minRank = -1;
    leftNode->bottOffDiag_minRank = -1;
    leftNode->LR_Method = "PS_Boundary";

    usrTree.setChildren_NULL(leftNode);
    usrTree.setChildren_NULL(rightNode);
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"blade_12k_level_0"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/2D/input/blade_12k_level_0","boundaryLR/stiffness/unstructured/2D/input/blade_12k_level_0_Graph","boundaryLR/stiffness/unstructured/2D/results_timing/blade_12k_level_0_",1e-1,30,2,usrTree);
     
    // Create 16K userTree
    usrTree.rootNode = new user_IndexTree::node;
    usrTree.rootNode->splitIndex = 8998;
    usrTree.rootNode->topOffDiag_minRank = -1;
    usrTree.rootNode->bottOffDiag_minRank = -1;
    usrTree.rootNode->LR_Method = "PS_Boundary";

    leftNode = new user_IndexTree::node;
    rightNode = new user_IndexTree::node;
    usrTree.rootNode->left = leftNode;
    usrTree.rootNode->right = rightNode;
   
    rightNode->splitIndex = 13932;
    rightNode->topOffDiag_minRank = -1;
    rightNode->bottOffDiag_minRank = -1;
    rightNode->LR_Method = "PS_Boundary";

    leftNode->splitIndex = 4499;
    leftNode->topOffDiag_minRank = -1;
    leftNode->bottOffDiag_minRank = -1;
    leftNode->LR_Method = "PS_Boundary";

    usrTree.setChildren_NULL(leftNode);
    usrTree.setChildren_NULL(rightNode);
   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"blade_16k_level_0"<<std::endl;
    testBoundaryLRSolver("boundaryLR/stiffness/unstructured/2D/input/blade_16k_level_0","boundaryLR/stiffness/unstructured/2D/input/blade_16k_level_0_Graph","boundaryLR/stiffness/unstructured/2D/results_timing/blade_16k_level_0_",1e-1,30,2,usrTree);

    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"Low-rank accuracy tests for a 3D FETI structured mesh ........"<<std::endl;
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_5_level_0"<<std::endl;
    analyzeRank("boundaryLR/FETI/structured/400k/input/400_front_num_5_level_0","boundaryLR/FETI/structured/400k/input/400_front_num_5_level_0_Graph","boundaryLR/FETI/structured/400k/results_LR/400_front_num_5_level_0_",0,0,0,0,"topOffDiag");   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_4_level_1"<<std::endl;
    analyzeRank("boundaryLR/FETI/structured/400k/input/400_front_num_4_level_1","boundaryLR/FETI/structured/400k/input/400_front_num_4_level_1_Graph","boundaryLR/FETI/structured/400k/results_LR/400_front_num_4_level_1_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_3_level_1"<<std::endl;
    analyzeRank("boundaryLR/FETI/structured/400k/input/400_front_num_3_level_1","boundaryLR/FETI/structured/400k/input/400_front_num_3_level_1_Graph","boundaryLR/FETI/structured/400k/results_LR/400_front_num_3_level_1_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_2_level_2"<<std::endl;
    analyzeRank("boundaryLR/FETI/structured/400k/input/400_front_num_2_level_2","boundaryLR/FETI/structured/400k/input/400_front_num_2_level_2_Graph","boundaryLR/FETI/structured/400k/results_LR/400_front_num_2_level_2_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_1_level_2"<<std::endl;
    analyzeRank("boundaryLR/FETI/structured/400k/input/400_front_num_1_level_2","boundaryLR/FETI/structured/400k/input/400_front_num_1_level_2_Graph","boundaryLR/FETI/structured/400k/results_LR/400_front_num_1_level_2_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_0_level_2"<<std::endl;
    analyzeRank("boundaryLR/FETI/structured/400k/input/400_front_num_0_level_2","boundaryLR/FETI/structured/400k/input/400_front_num_0_level_2_Graph","boundaryLR/FETI/structured/400k/results_LR/400_front_num_0_level_2_",0,0,0,0,"topOffDiag");
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"Low-rank accuracy tests for a 3D FETI unstructured mesh ........"<<std::endl;
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_0_level_13"<<std::endl;
    analyzeRank("boundaryLR/FETI/unstructured/400k/input/400_front_num_0_level_13","boundaryLR/FETI/unstructured/400k/input/400_front_num_0_level_13_Graph","boundaryLR/FETI/unstructured/400k/results_LR/400_front_num_0_level_13_",0,0,0,0,"topOffDiag");   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_1_level_11"<<std::endl;
    analyzeRank("boundaryLR/FETI/unstructured/400k/input/400_front_num_1_level_11","boundaryLR/FETI/unstructured/400k/input/400_front_num_1_level_11_Graph","boundaryLR/FETI/unstructured/400k/results_LR/400_front_num_1_level_11_",0,0,0,0,"topOffDiag"); 
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_2_level_7"<<std::endl;
    analyzeRank("boundaryLR/FETI/unstructured/400k/input/400_front_num_2_level_7","boundaryLR/FETI/unstructured/400k/input/400_front_num_2_level_7_Graph","boundaryLR/FETI/unstructured/400k/results_LR/400_front_num_2_level_7_",0,0,0,0,"topOffDiag");   
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"400_front_num_3_level_6"<<std::endl;
    analyzeRank("boundaryLR/FETI/unstructured/400k/input/400_front_num_3_level_6","boundaryLR/FETI/unstructured/400k/input/400_front_num_3_level_6_Graph","boundaryLR/FETI/unstructured/400k/results_LR/400_front_num_3_level_6_",0,0,0,0,"topOffDiag"); 
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"Low-rank accuracy tests for a 3D unstructured mesh........."<<std::endl;
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_6_level_0"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_6_level_0","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_6_level_0_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_LR/300_front_num_6_level_0_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_5_level_1"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_5_level_1","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_5_level_1_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_LR/300_front_num_5_level_1_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_4_level_1"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_4_level_1","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_4_level_1_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_LR/300_front_num_4_level_1_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_3_level_2"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_3_level_2","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_3_level_2_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_LR/300_front_num_3_level_2_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_2_level_2"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_2_level_2","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_2_level_2_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_LR/300_front_num_2_level_2_",0,0,0,0,"topOffDiag");    
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_1_level_2"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_1_level_2","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_1_level_2_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_LR/300_front_num_1_level_2_",0,0,0,0,"topOffDiag");
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"300_front_num_0_level_2"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_0_level_2","boundaryLR/stiffness/unstructured/3D/300k/input/300_front_num_0_level_2_Graph","boundaryLR/stiffness/unstructured/3D/300k/results_LR/300_front_num_0_level_2_",0,0,0,0,"topOffDiag");    
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"*******************************************************"<<std::endl;
    std::cout<<"Low-rank accuracy tests for a 2D unstructured mesh........."<<std::endl;
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"blade_9k_level_0"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/2D/input/blade_9k_level_0","boundaryLR/stiffness/unstructured/2D/input/blade_9k_level_0_Graph","boundaryLR/stiffness/unstructured/2D/results_LR/blade_9k_level_0_",1499,3997,3996-1499+1,6497-3997+1);
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"blade_12k_level_0"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/2D/input/blade_12k_level_0","boundaryLR/stiffness/unstructured/2D/input/blade_12k_level_0_Graph","boundaryLR/stiffness/unstructured/2D/results_LR/blade_12k_level_0_",2099,5598,5597-2099+1,9072-5598+1);
    std::cout<<"-----------------------"<<std::endl;
    std::cout<<"blade_16k_level_0"<<std::endl;
    analyzeRank("boundaryLR/stiffness/unstructured/2D/input/blade_16k_level_0","boundaryLR/stiffness/unstructured/2D/input/blade_16k_level_0_Graph","boundaryLR/stiffness/unstructured/2D/results_LR/blade_16k_level_0_",0,4500,4500,8998-4500+1);

  }
  return 0;
}
