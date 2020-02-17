#!/bin/bash

# packages from distribution repositories
apt update && apt install --no-install-recommends -y g++ libhdf5-dev python3

# distribution specialization
apt install --no-install-recommends -y lsb-release
codename=`lsb_release -cs`

if [ "$codename" == "buster" ]; then
    apt-get install --no-install-recommends -y clang-7 libomp5
fi
