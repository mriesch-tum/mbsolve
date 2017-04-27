cmake_minimum_required(VERSION 2.6)
project(mbsolve-lib)

find_package(Boost REQUIRED COMPONENTS program_options system timer)

include_directories(include/ ${Boost_INCLUDE_DIR})

add_library(mbsolve-lib SHARED src/solver.cpp src/writer.cpp)
target_include_directories(mbsolve-lib PUBLIC include/)
