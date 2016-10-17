cmake_minimum_required(VERSION 2.6)
project(solver-cuda)

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

cuda_add_library(solver-cuda SHARED src/SolverCUDA.cu)

target_link_libraries(solver-cuda mbsolve-lib)
