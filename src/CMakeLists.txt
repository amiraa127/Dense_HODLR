# Directories to look for include files
include_directories(${DenseHODLR_SOURCE_DIR}/include)
include_directories(${DenseHODLR_SOURCE_DIR}/IML/include)
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${PASTIX_INCLUDE_DIRS})


set ( src_files "" )

####################################################
# Find all C/C++ files recursively from current dir
####################################################
foreach ( ext IN ITEMS "cpp" "cxx" "cc" "c" )
	file ( GLOB_RECURSE _files "*.${ext}" )
	set ( src_files ${src_files} ${_files} )
endforeach ()

#message ( STATUS "Found: ${src_files}" )
    
# Create a library called HODLRSolver which includes the source files defined
add_library(HODLRSolver ${src_files})
