cmake_minimum_required(VERSION 2.8)
project(DenseHODLR)
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Enable test running via 'make test' command
enable_testing()

set(GCC_COVERAGE_COMPILE_FLAGS " -std=c++11 -O3 -Wall")

add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

if (NOT PASTIX)
   add_definitions(-DNOPASTIX)
   #message("No PaStiX")
else()
   find_package (Pastix  REQUIRED)
endif(NOT PASTIX)

# Find the dependencies.
find_package (Eigen REQUIRED)
find_package (CppUnit REQUIRED)

if (PASTIX_FOUND)
   message ("-- Found Pastix version ${PASTIX_VERSION}: ${PASTIX_INCLUDE_DIRS}")
endif (PASTIX_FOUND)

if (EIGEN_FOUND)
    message ("-- Found Eigen version ${EIGEN_VERSION} : ${EIGEN_INCLUDE_DIRS}")
endif (EIGEN_FOUND)

if (CPPUNIT_FOUND)
    message ("-- Found CPPUNIT : ${CPPUNIT_INCLUDE_DIRS}")
endif (CPPUNIT_FOUND)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(benchmarks)
add_subdirectory(rigidBody)

