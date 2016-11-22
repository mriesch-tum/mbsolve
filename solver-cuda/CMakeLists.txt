cmake_minimum_required(VERSION 2.6)
project(solver-cuda)

# TODO: it seems the FindCUDA module does not automatically include the header
# directory. As a fix it is provided manually.
#include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
include_directories(../mbsolve-lib/include/ include/)

# TODO: compiler options are not passed by CMake, one has to manually fix it
# with this version
list(APPEND CUDA_NVCC_FLAGS --compiler-options -fPIC -shared)
#list(APPEND CUDA_NVCC_FLAGS -rdc=true)
set(CUDA_SEPARABLE_COMPILATION ON)

cuda_add_library(solver-cuda SHARED src/CUDADensityMatrix.cu
  src/SolverCUDA2lvl.cu include/SolverCUDA2lvl.hpp
  include/CUDADensityMatrix.hpp include/CUDACommon.hpp)

target_include_directories(solver-cuda PUBLIC include/)
target_link_libraries(solver-cuda mbsolve-lib)
