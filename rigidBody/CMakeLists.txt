# Directories to look for include files
include_directories(${DenseHODLR_SOURCE_DIR}/include)
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${PASTIX_INCLUDE_DIRS})

# Directories to look for libraries
link_directories(${DenseHODLR_BINARY_DIR}/src)

# Compile tests from singleBody.cpp
add_executable(singleBodyTest singleBody.cpp)
add_executable(latticeTest lattice.cpp)

# test should link to solver libraries
target_link_libraries(singleBodyTest HODLRSolver ${PASTIX_LIBRARIES})
target_link_libraries(latticeTest HODLRSolver ${PASTIX_LIBRARIES})

