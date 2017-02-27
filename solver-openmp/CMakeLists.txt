cmake_minimum_required(VERSION 2.6)
project(solver-openmp)

#string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
#set(CMAKE_CXX_FLAGS_DEBUG "-O3 -g")

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

add_library(solver-openmp SHARED src/solver_2lvl_pc.cpp)

#target_compile_options(solver-openmp PRIVATE "-xCORE-AVX2")
target_compile_options(solver-openmp PRIVATE ${OpenMP_CXX_FLAGS})

target_include_directories(solver-openmp PUBLIC include/)
target_link_libraries(solver-openmp ${OpenMP_CXX_FLAGS})
target_link_libraries(solver-openmp mbsolve-lib)
