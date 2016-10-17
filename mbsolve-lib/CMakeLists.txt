cmake_minimum_required(VERSION 2.6)
project(mbsolve-lib)

include_directories(include/)

add_library(mbsolve-lib SHARED src/Solver.cpp)
target_include_directories(mbsolve-lib PUBLIC include/)
#target_link_libraries(mbsolve-lib solver-generic)
