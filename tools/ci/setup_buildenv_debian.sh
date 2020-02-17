#!/bin/bash

# packages from distribution repositories
apt update && apt install --no-install-recommends -y ca-certificates \
  doxygen g++ git libhdf5-dev make python3-dev python3-pip \
  python3-setuptools python3-wheel swig wget

# distribution specialization
apt install --no-install-recommends -y lsb-release
codename=`lsb_release -cs`

# Clang compiler
if [ "$codename" == "buster" ]; then
    apt-get install --no-install-recommends -y clang-7 libomp-dev
fi

# clang-format
if [ "$codename" == "buster" ]; then
    apt-get install --no-install-recommends -y clang-format-6.0
fi

# CMake
if [ "$codename" == "jessie" ] || [ "$codename" == "stretch" ]; then
    # we need CMake >= 3.9
    wget https://github.com/Kitware/CMake/releases/download/v3.16.2/cmake-3.16.2-Linux-x86_64.sh
    chmod +x cmake-3.16.2-Linux-x86_64.sh
    ./cmake-3.16.2-Linux-x86_64.sh --skip-license --prefix=/
else
    apt-get install --no-install-recommends -y cmake
fi

# cxxopts library
git clone https://github.com/jarro2783/cxxopts.git && cd cxxopts && \
  cmake -Bbuild -H. -DCXXOPTS_BUILD_EXAMPLES=OFF -DCXXOPTS_BUILD_TESTS=OFF && \
  cmake --build build/ --target install && cd ..

# Eigen3 library
git clone https://gitlab.com/libeigen/eigen.git && cd eigen && \
  git checkout 3.3.7 && cmake -Bbuild -H. -DBUILD_TESTING=OFF && \
  cmake --build build/ --target install && cd ..
