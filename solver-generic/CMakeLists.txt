cmake_minimum_required(VERSION 3.6)
project(solver-generic)

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

add_library(solver-generic SHARED src/solver_generic.cpp)
set_property(TARGET solver-generic PROPERTY CXX_STANDARD 11)

target_include_directories(solver-generic PUBLIC include/)
target_link_libraries(solver-generic mbsolve-lib)
