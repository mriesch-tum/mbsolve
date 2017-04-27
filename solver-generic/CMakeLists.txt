cmake_minimum_required(VERSION 3.6)
project(solver-generic)

find_package(Boost REQUIRED COMPONENTS program_options system timer)

include_directories(${mbsolve-lib_INCLUDE_DIR} include/ ${Boost_INCLUDE_DIR})

add_library(solver-generic SHARED src/SolverGeneric.cpp)
target_include_directories(solver-generic PUBLIC include/)
target_link_libraries(solver-generic mbsolve-lib)
