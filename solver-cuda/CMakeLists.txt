cmake_minimum_required(VERSION 3.6)
project(solver-cuda)

# TODO: it seems the FindCUDA module does not automatically include the header
# directory. As a fix it is provided manually.
#include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
include_directories(../mbsolve-lib/include/ include/)
#include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

# TODO: arch as build option?-arch=sm_35;
list(APPEND CUDA_NVCC_FLAGS --std=c++11; -Xcompiler -fPIC)
#list(APPEND CUDA_NVCC_FLAGS -rdc=true)
set(CUDA_SEPARABLE_COMPILATION ON)

cuda_add_library(solver-cuda
  src/common_cuda.cu include/common_cuda.hpp)

set_property(TARGET solver-cuda PROPERTY CXX_STANDARD 11)

target_include_directories(solver-cuda PUBLIC include/)
target_link_libraries(solver-cuda mbsolve-lib)
