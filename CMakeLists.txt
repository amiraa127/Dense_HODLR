cmake_minimum_required(VERSION 2.8)
project(DenseHODLR)
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Enable test running via 'make test' command
enable_testing()

set(GCC_COVERAGE_COMPILE_FLAGS "-std=c++0x -O3 -Wall")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

# Find the dependencies.
find_package (Eigen REQUIRED)
find_package (CppUnit REQUIRED)

if (EIGEN_FOUND)
    message ("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
endif (EIGEN_FOUND)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(benchmarks)
add_subdirectory(rigidBody)

file(COPY test/data DESTINATION ${DenseHODLR_BINARY_DIR}/test)
file(COPY benchmarks/ACA DESTINATION ${DenseHODLR_BINARY_DIR}/benchmarks)
file(COPY benchmarks/ACA_Blade DESTINATION ${DenseHODLR_BINARY_DIR}/benchmarks)
file(COPY benchmarks/boundaryLR DESTINATION ${DenseHODLR_BINARY_DIR}/benchmarks)
file(COPY rigidBody/input DESTINATION ${DenseHODLR_BINARY_DIR}/rigidBody/)