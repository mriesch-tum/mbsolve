cmake_minimum_required(VERSION 3.6)
project(solver-cuda)

# TODO: it seems the FindCUDA module does not automatically include the header
# directory. As a fix it is provided manually.
#include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
#include_directories(../mbsolve-lib/include/ include/)
include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

# TODO: arch as build option?
list(APPEND CUDA_NVCC_FLAGS -arch=sm_35 --compiler-flags)
#list(APPEND CUDA_NVCC_FLAGS -rdc=true)
set(CUDA_SEPARABLE_COMPILATION ON)

cuda_add_library(solver-cuda SHARED
  # src/CUDADensityMatrix.cu include/CUDADensityMatrix.hpp
  src/solver_cuda_2lvl_pc.cu include/solver_cuda_2lvl_pc.hpp
#  src/solver_cuda_2lvl_pc_red.cu include/solver_cuda_2lvl_pc_red.hpp
  include/cuda_common.hpp)

# TODO: how to ensure c++11 properly?
# TODO: use newer cuda version?? could also solve arch issue above
set_property(TARGET solver-cuda PROPERTY CXX_STANDARD 11)

target_include_directories(solver-cuda PUBLIC include/)
target_link_libraries(solver-cuda mbsolve-lib)
