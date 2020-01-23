#
# mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
#
# Copyright (c) 2016, Computational Photonics Group, Technical University of
# Munich.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

cmake_minimum_required(VERSION 3.6)
project(solver-cuda)

# TODO: it seems the FindCUDA module does not automatically include the header
# directory. As a fix it is provided manually.
#include_directories(${mbsolve-lib_INCLUDE_DIR} include/)
include_directories(../mbsolve-lib/include/ include/)
#include_directories(${mbsolve-lib_INCLUDE_DIR} include/)

# TODO: arch as build option?-arch=sm_35;
list(APPEND CUDA_NVCC_FLAGS --std=c++11; -Xcompiler -fPIC)
#list(APPEND CUDA_NVCC_FLAGS -rdc=true)
set(CUDA_SEPARABLE_COMPILATION ON)

cuda_add_library(solver-cuda
  src/common_cuda.cu include/common_cuda.hpp)

set_property(TARGET solver-cuda PROPERTY CXX_STANDARD 11)

target_include_directories(solver-cuda PUBLIC include/)
target_link_libraries(solver-cuda mbsolve-lib)
