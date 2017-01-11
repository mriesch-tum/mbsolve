cmake_minimum_required(VERSION 2.6)
project(solver-openmp)

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

add_library(solver-openmp SHARED src/solver_2lvl_pc.cpp src/DensityMatrixData.cpp)

target_compile_options(solver-openmp PRIVATE ${OpenMP_CXX_FLAGS})

target_include_directories(solver-openmp PUBLIC include/)
target_link_libraries(solver-openmp ${OpenMP_CXX_FLAGS})
target_link_libraries(solver-openmp mbsolve-lib)
