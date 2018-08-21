# Documentation

## Installing Channelflow

Channelflow is hosted on [GitHub](https://github.com/epfl-ecps/channelflow), using the distributed  version control system git.
To obtain a copy of the code clone the current repository running

`git clone https://github.com/epfl-ecps/channelflow`


### Prerequisites

To compile Channelflow:
* CMake (version 3.1 or higher)
* C++ supportive compiler
* FFTW (version 3 or higher)
* Eigen3 (version 3 or higher)

To enable parallelization, to use the NetCDF format (easy visualization with Paraview, VisIt ...) and to enable the parallel I/O:
* MPI
* NetCDF (parallel version 4 or higher)
* HDF5(cxx) (format available for backwards compatibility)

To use Channelflow functions from the Python wrapper:
* boost-python


### Compilation
A Makefile to compile and install Channelflow can be generated using the provided CMake build scripts.
CMake allows two build types: "release" and "debug". The build type decides on compiler options and whether assertions
are compiled. Use "release" for optimal performance (e.g. for production runs) and "debug" for optimal diagnostics
(e.g. for code development).
Out of source builds are recommended.

```
mkdir build
cd build
cmake PATH_TO_SOURCE -DCMAKE_BUILD_TYPE=debug (/release) (configuration options)
make -j
make install
```

Channelflow supports, beneath other standard cmake flags, the following options


|Option                   | Values  | Default   | Description                                                       |
|:------------------------|:--------|:----------|:------------------------------------------------------------------|
|`-DCMAKE_INSTALL_PREFIX` | path    | usr/local | Installation path for make install                                |
|`-DUSE_MPI`              | ON/OFF  | ON        | Enable MPI                                                        |
|`-DWITH_SHARED`          | ON/OFF  | ON        | build shared channelflow and nsolver libraries                    |
|`-DWITH_STATIC`          | ON/OFF  | OFF       | build static libraries (also enables linking to static libraries) |
|`-DWITH_PYTHON`          | ON/OFF  | OFF       | build a python wrapper for flowfields, disabled by default because it requires boost-python |
|`-DWITH_HDF5CXX`         |  ON/OFF | OFF       | enable legacy .h5 file format (using HDF5 C++)                    |


A complete installation, with all features enabled, might look like this:

  `cmake PATH_TO_SOURCE -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=$HOME/usr -DWITH_PYTHON=ON -DWITH_HDF5CXX=ON`


### Tests

The correct functionality of the Channelflow code is tested in two ways:
* **Unit tests** are tests of individual classes or methods.
* **Integration tests** are tests of applications with all underlying code.

Both test types are hooked to the make test command, where they run sequentially.
For verbose output and debugging, run `make test ARGS=“-V”`

To run only the unit tests, run `tests/gtest/runUnitTest` in the CMake build directory.Pass argument <br/>
`–gtest_filter=<foo>`

Integration tests are found in folder tests, and registered in tests/CMakeList.txt. They are compiled into standalong
programs, which are individually listed in tests/CMakeList.txt, from where they are called on make test.



### Running your first simulations


To run a DNS we first need to provide an initial velocity field. This can be done providing a specific file or creating
a random flowfield using the utility `randomfield` in the `tools` folder.

Running the command line

`randomfield -Nx 48 -Ny 81 -Nz 32 -Lx 3 -Lz 2 newfield.nc`

a random flowfield with zero-divergence, Dirichlet boundary conditions and the specified resolution and geometry is created.

Time integration of this is initial condition is performed by the program `simulateflow` in the `programs` folder.
The specific system to be simulated is specified via the program's options.

The following example integrates for 200 advective time units a simple plane Couette system at the default Reynolds number
400 .

`simulateflow -T 200 newfield.nc`

The following example simulates a pressure driven (-dPds) Poiseuille flow (-Uwall 0) at the Reynolds number 1000 using
a 1st order forward-backward euler scheme as initial time stepping algorithm (-is).

`simulateflow -R 1000 -Uwall 0 -dPds -0.002 -is SBDF1 newfield.nc`



