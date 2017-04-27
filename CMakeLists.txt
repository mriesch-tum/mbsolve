cmake_minimum_required(VERSION 3.6)
project(mbsolve)

set(MBSOLVE_VERSION_MAJOR 0)
set(MBSOLVE_VERSION_MINOR 1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

# compiler options
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -mavx2")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -DNDEBUG -g")
endif()

option(WITH_MATLAB "with(out) MATLAB writer" ON)
option(WITH_CUDA "with(out) CUDA solver" ON)
option(WITH_OPENMP "with(out) OPENMP solver" ON)
option(WITH_DOXYGEN "with(out) Doxygen documentation" ON)

# mbsolve library
add_subdirectory(mbsolve-lib)

# MATLAB support?
if(WITH_MATLAB)
  find_package(MATLAB)
  if(MATLAB_FOUND)
    set(MBSOLVE_MATLAB TRUE)
    add_subdirectory(writer-matlab)
  endif(MATLAB_FOUND)
endif(WITH_MATLAB)

# CUDA support?
if(WITH_CUDA)
  find_package(CUDA)
  if(CUDA_FOUND)
    set(MBSOLVE_CUDA TRUE)
    add_subdirectory(solver-cuda)
  endif(CUDA_FOUND)
endif(WITH_CUDA)

# OpenMP support?
if(WITH_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(MBSOLVE_OPENMP TRUE)
    add_subdirectory(solver-openmp)
  endif(OPENMP_FOUND)
endif(WITH_OPENMP)

# Generic support
add_subdirectory(solver-generic)

# executable
add_subdirectory(mbsolve-tool)

# Doxygen documentation
if(WITH_DOXYGEN)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM)
  endif(DOXYGEN_FOUND)
endif(WITH_DOXYGEN)
