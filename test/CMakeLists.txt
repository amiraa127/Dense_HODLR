# Directories to look for include files
include_directories(${DenseHODLR_SOURCE_DIR}/include)
include_directories(${DenseHODLR_SOURCE_DIR}/IML/include/)
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${CPPUNIT_INCLUDE_DIRS})
include_directories(${PASTIX_INCLUDE_DIRS})

# Directories to look for libraries
link_directories(${DenseHODLR_BINARY_DIR}/src)

# Compile tests from test.cpp
add_executable(HODLRTest test.cpp)

# test should link to solver libraries
target_link_libraries(HODLRTest HODLRSolver ${CPPUNIT_LIBRARIES} ${PASTIX_LIBRARIES})

# Tests to run
add_test(BasicTests HODLRTest)

