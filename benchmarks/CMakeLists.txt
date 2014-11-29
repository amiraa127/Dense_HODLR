# Directories to look for include files
include_directories(${DenseHODLR_SOURCE_DIR}/include)
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${PASTIX_INCLUDE_DIRS})


# Directories to look for libraries
link_directories(${DenseHODLR_BINARY_DIR}/src)

# Compile benchmarks from benchmarks.cpp
add_executable(HODLRBenchmarks benchmarks.cpp)

# test should link to solver libraries
target_link_libraries(HODLRBenchmarks HODLRSolver ${PASTIX_LIBRARIES})


