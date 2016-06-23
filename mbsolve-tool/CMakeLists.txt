cmake_minimum_required(VERSION 2.6)
project(mbsolve-tool)

find_package(Boost REQUIRED COMPONENTS program_options system timer)

include_directories(${MBSOLVE_INCLUDE_DIR}
  ${SOLVERCUDA_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR})

add_executable(mbsolve-tool src/mbsolve-tool.cpp)

target_link_libraries(mbsolve-tool ${Boost_LIBRARIES})
target_link_libraries(mbsolve-tool solver-cuda)
