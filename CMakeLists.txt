cmake_minimum_required(VERSION 2.8)
project(DenseHODLR)

# Enable test running via 'make test' command
enable_testing()

set(GCC_COVERAGE_COMPILE_FLAGS "-std=c++0x -O3 -Wall")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(benchmarks)
add_subdirectory(rigidBody)

file(COPY test/data DESTINATION ${DenseHODLR_BINARY_DIR}/test)
file(COPY benchmarks/ACA DESTINATION ${DenseHODLR_BINARY_DIR}/benchmarks)
file(COPY benchmarks/ACA_Blade DESTINATION ${DenseHODLR_BINARY_DIR}/benchmarks)
file(COPY benchmarks/boundaryLR DESTINATION ${DenseHODLR_BINARY_DIR}/benchmarks)
file(COPY rigidBody/input DESTINATION ${DenseHODLR_BINARY_DIR}/rigidBody/)
