cmake_minimum_required(VERSION 2.6)
project(mbsolve)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

option(WITH_MATLAB "with(out) MATLAB writer" ON)
option(WITH_CUDA "with(out) CUDA solver" ON)

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

# Generic support
add_subdirectory(solver-generic)

# executable
add_subdirectory(mbsolve-tool)
