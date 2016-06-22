cmake_minimum_required(VERSION 2.6)
project(solver-cuda)

find_package(CUDA REQUIRED)

include_directories(${MBSOLVE_INCLUDE_DIR}
  ${SOLVERCUDA_INCLUDE_DIR})

cuda_add_library(solver-cuda SHARED src/SolverCUDA.cu)
