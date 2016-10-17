cmake_minimum_required(VERSION 2.6)
project(solver-cuda)

# TODO: it seems the FindCUDA module does not automatically include the header
# directory. As a fix it is provided manually.
#include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
include_directories(../mbsolve-lib/include/ include/)

cuda_add_library(solver-cuda SHARED
  src/SolverCUDA2lvl.cu include/SolverCUDA2lvl.hpp)

target_include_directories(solver-cuda PUBLIC include/)
target_link_libraries(solver-cuda mbsolve-lib)
