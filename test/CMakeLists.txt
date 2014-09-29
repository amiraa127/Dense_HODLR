# Directories to look for include files
include_directories(${DenseHODLR_SOURCE_DIR}/include)
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${CPPUNIT_INCLUDE_DIRS})

# Add CPPUNIT
add_library(libCPPUNIT STATIC IMPORTED)
set_target_properties(libCPPUNIT PROPERTIES IMPORTED_LOCATION /opt/local/lib/libcppunit.a)

# Directories to look for libraries
link_directories(${DenseHODLR_BINARY_DIR}/src)

# Compile tests from test.cpp
add_executable(HODLRTest test.cpp)

# test should link to solver libraries
target_link_libraries(HODLRTest HODLRSolver libCPPUNIT)

# Tests to run
add_test(BasicTests HODLRTest)

