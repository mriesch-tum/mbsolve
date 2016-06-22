cmake_minimum_required(VERSION 2.6)
project(mbsolve)

set(MBSOLVE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SOLVERCUDA_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/solver-cuda/include)

# solver libraries
add_subdirectory(solver-cuda)

# executable
add_subdirectory(mbsolve-tool)
