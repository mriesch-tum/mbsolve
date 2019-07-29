cmake_minimum_required(VERSION 3.6)
project(solver-cpu)

string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELWITHDEBINFO
  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

add_library(solver-cpu SHARED
  src/solver_cpu_fdtd.cpp
  src/solver_cpu_fdtd_red.cpp)

set_property(TARGET solver-cpu PROPERTY CXX_STANDARD 11)

if(ARCH STREQUAL "AVX2")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-cpu PRIVATE "-mavx2")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-cpu PRIVATE "-xAVX2")
  else()
  endif()
elseif(ARCH STREQUAL "SSE3")
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-cpu PRIVATE "-msse3")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-cpu PRIVATE "-xSSE3")
  else()
  endif()
else()
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_options(solver-cpu PRIVATE "-march=native")
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    target_compile_options(solver-cpu PRIVATE "-xHost")
  else()
  endif()
endif()

target_compile_options(solver-cpu PRIVATE ${OpenMP_CXX_FLAGS})

target_include_directories(solver-cpu PUBLIC include/)
target_link_libraries(solver-cpu Eigen3::Eigen)
target_link_libraries(solver-cpu ${OpenMP_CXX_FLAGS})
target_link_libraries(solver-cpu mbsolve-lib)

if(MBSOLVE_PYTHON)
  include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
  include_directories(${PYTHON_INCLUDE_PATH})
  set_source_files_properties(swig/pysolvercpu.i PROPERTIES CPLUSPLUS ON)
  swig_add_module(pysolvercpu python swig/pysolvercpu.i)
  set_property(TARGET ${SWIG_MODULE_pysolvercpu_REAL_NAME} PROPERTY
    CXX_STANDARD 11)
  swig_link_libraries(pysolvercpu ${PYTHON_LIBRARIES} mbsolve-lib
    solver-cpu)
endif(MBSOLVE_PYTHON)
