#ifndef EIGEN_PASTIX_SUPPORT_MODULE_HPP
#define EIGEN_PASTIX_SUPPORT_MODULE_HPP

#include "Eigen/SparseCore"
#include "Eigen/src/Core/util/DisableStupidWarnings.h"

#include <complex.h>
extern "C" {
#include <pastix_nompi.h>
#include <pastix.h>
}

#ifdef complex
#undef complex
#endif

#include "Eigen/src/misc/Solve.h"
#include "Eigen/src/misc/SparseSolve.h"

#define  COMPLEX  std::complex<float>
#define  DCOMPLEX std::complex<double>

#include "Eigen/src/PaStiXSupport/PaStiXSupport.h"
#include "Eigen/src/Core/util/ReenableStupidWarnings.h"

#endif
