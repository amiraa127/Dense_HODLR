# - Try to find PaStiX
# Once done this will define
#
#  PASTIX_FOUND        - system has PaStiX
#  PASTIX_INCLUDE_DIRS - include directories for PaStiX
#  PASTIX_LIBRARIES    - libraries for PaStiX
#  PASTIX_VERSION      - the PaStiX version string (MAJOR.MEDIUM.MINOR)

# Check for PaStiX header file
find_path(PASTIX_INCLUDE_DIR pastix.h
  HINTS ${PASTIX_DIR}/include $ENV{PASTIX_DIR}/include ${PASTIX_DIR} $ENV{PASTIX_DIR}
  PATH_SUFFIXES install
  DOC "Directory where the PaStiX header is located"
 )

set(PASTIX_INCLUDE_DIRS ${PASTIX_INCLUDE_DIR})

# Check for PaStiX library
find_library(PASTIX_LIBRARY pastix
  HINTS ${PASTIX_DIR}/lib $ENV{PASTIX_DIR}/lib ${PASTIX_DIR} $ENV{PASTIX_DIR}
  PATH_SUFFIXES install
  DOC "The PaStiX library"
  )

set(PASTIX_LIBRARIES ${PASTIX_LIBRARY})

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT APPLE)
  # Check for rt library
  find_library(RT_LIBRARY rt
    DOC "The RT library"
    )
  # Check for math library
  find_library(M_LIBRARY m
    DOC "The math library"
    )
  set(PASTIX_LIBRARIES ${PASTIX_LIBRARIES} ${RT_LIBRARY} ${M_LIBRARY})
  message(STATUS "PASTIX_LIBRARIES ${PASTIX_LIBRARIES}")
endif()

# Check for hwloc header
find_path(HWLOC_INCLUDE_DIRS hwloc.h
  HINTS ${HWLOC_DIR} $ENV{HWLOC_DIR} ${HWLOC_DIR}/include $ENV{HWLOC_DIR}/include
  DOC "Directory where the hwloc header is located"
 )

# Check for hwloc library
find_library(HWLOC_LIBRARY hwloc
  HINTS ${HWLOC_DIR} $ENV{HWLOC_DIR} ${HWLOC_DIR}/lib $ENV{HWLOC_DIR}/lib
  DOC "The hwloc library"
  )

if (HWLOC_LIBRARY)
  set(PASTIX_LIBRARIES ${PASTIX_LIBRARIES} ${HWLOC_LIBRARY})
endif()

# Add BLAS libs if BLAS has been found
set(CMAKE_LIBRARY_PATH ${BLAS_DIR}/lib $ENV{BLAS_DIR}/lib ${CMAKE_LIBRARY_PATH})
find_package(BLAS)
set(PASTIX_LIBRARIES ${PASTIX_LIBRARIES} ${BLAS_LIBRARIES})

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  find_program(GFORTRAN_EXECUTABLE gfortran)
  if (GFORTRAN_EXECUTABLE)
    execute_process(COMMAND ${GFORTRAN_EXECUTABLE} -print-file-name=libgfortran.so
      OUTPUT_VARIABLE GFORTRAN_LIBRARY
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if (EXISTS "${GFORTRAN_LIBRARY}")
      set(PASTIX_LIBRARIES ${PASTIX_LIBRARIES} ${GFORTRAN_LIBRARY})
    endif()
  endif()
endif()

#Find PaStiX version by looking at pastix.h
if (EXISTS ${PASTIX_INCLUDE_DIR}/pastix.h)
  file(STRINGS ${PASTIX_INCLUDE_DIR}/pastix.h PASTIX_VERSIONS_TMP
    REGEX "^#define PASTIX_[A-Z]+_VERSION[ \t]+[0-9]+$")
  if (PASTIX_VERSIONS_TMP)
    string(REGEX REPLACE ".*#define PASTIX_MAJOR_VERSION[ \t]+([0-9]+).*" "\\1" PASTIX_MAJOR_VERSION ${PASTIX_VERSIONS_TMP})
    string(REGEX REPLACE ".*#define PASTIX_MEDIUM_VERSION[ \t]+([0-9]+).*" "\\1" PASTIX_MEDIUM_VERSION ${PASTIX_VERSIONS_TMP})
    string(REGEX REPLACE ".*#define PASTIX_MINOR_VERSION[ \t]+([0-9]+).*" "\\1" PASTIX_MINOR_VERSION ${PASTIX_VERSIONS_TMP})
    set(PASTIX_VERSION ${PASTIX_MAJOR_VERSION}.${PASTIX_MEDIUM_VERSION}.${PASTIX_MINOR_VERSION} CACHE TYPE STRING)
  else()
    set(PASTIX_VERSION "UNKNOWN")
  endif()
endif()

if (PaStiX_FIND_VERSION)
  # Check if version found is >= required version
  if (NOT "${PASTIX_VERSION}" VERSION_LESS "${PaStiX_FIND_VERSION}")
    set(PASTIX_VERSION_OK TRUE)
  endif()
else()
  # No specific version requested
  set(PASTIX_VERSION_OK TRUE)
endif()

if (PASTIX_VERSION_OK)
   set(PASTIX_FOUND TRUE)
endif()

mark_as_advanced(
  PASTIX_INCLUDE_DIR
  PASTIX_INCLUDE_DIRS
  PASTIX_LIBRARY
  PASTIX_LIBRARIES
  PASTIX_VERSION
  PASTIX_FOUND
  )

