cmake_minimum_required(VERSION 2.6)
project(mbsolve)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

set(MBSOLVE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(SOLVERCUDA_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/solver-cuda/include)
set(WRITERMATLAB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/writer-matlab/include)

# writer libraries
add_subdirectory(writer-matlab)

# solver libraries
add_subdirectory(solver-cuda)

# executable
add_subdirectory(mbsolve-tool)
