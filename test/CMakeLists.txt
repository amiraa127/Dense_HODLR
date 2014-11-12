# Directories to look for include files
include_directories(${DenseHODLR_SOURCE_DIR}/include)
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${CPPUNIT_INCLUDE_DIRS})

# Directories to look for libraries
link_directories(${DenseHODLR_BINARY_DIR}/src)

# Compile tests from test.cpp
add_executable(HODLRTest test.cpp)

# test should link to solver libraries
target_link_libraries(HODLRTest HODLRSolver CppUnit)

# Tests to run
add_test(BasicTests HODLRTest)

