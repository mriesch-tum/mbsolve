cmake_minimum_required(VERSION 3.6)
project(mbsolve-lib)

include_directories(include/)

add_library(mbsolve-lib SHARED src/solver.cpp src/writer.cpp
  src/Quantity.cpp)
set_property(TARGET mbsolve-lib PROPERTY CXX_STANDARD 11)

target_include_directories(mbsolve-lib PUBLIC include/)
