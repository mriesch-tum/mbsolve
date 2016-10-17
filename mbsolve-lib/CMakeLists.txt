cmake_minimum_required(VERSION 2.6)
project(mbsolve-lib)

include_directories(include/)

add_library(mbsolve-lib SHARED
  src/Quantity.cpp src/Solver.cpp src/Writer.cpp
  include/Device.hpp include/Quantity.hpp
  include/Scenario.hpp include/Solver.hpp
  include/Types.hpp include/Writer.hpp)
target_include_directories(mbsolve-lib PUBLIC include/)
