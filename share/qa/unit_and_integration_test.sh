#!/bin/bash -e
#
# Description:
#    This script is used to run unit and integration tests in Travis

# Print ccache version and re-initialize statistics
ccache -V
ccache -z

# Parallel build in Release mode
mkdir ${TRAVIS_BUILD_DIR}/cmake-build-release && cd ${TRAVIS_BUILD_DIR}/cmake-build-release
CC="ccache mpicc" CXX="ccache mpicxx" cmake -DCMAKE_BUILD_TYPE:STRING=Release -DWITH_HDF5CXX:BOOL=ON -DWITH_PYTHON:BOOL=OFF -DBUILD_SHARED_LIBS:BOOL=ON -DWITH_NSOLVER:BOOL=ON -DWITH_NETCDF:STRING=Serial ${TRAVIS_BUILD_DIR}
make -j 4

# Print ccache statistics
ccache -s

# Run tests
make test


