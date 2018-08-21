# Channelflow

Channelflow is a software system for numerical analysis of the incompressible fluid flow in 
channel geometries, written in C++. Channelflow does time integration of plane Couette and 
plane Poiseuille flow using a Fourier x Chebyshev x Fourier spatial discretization, finite 
differencing in time, and multiprocessor parallelism via MPI. Channelflow also computes 
invariant solutions of these flows: equilibria, traveling waves, and periodic orbits. 

The main goals of channelflow are
  * to lower the barrier to entry to numerical research in fluids
  * to enable creation of short, readable, easily-modifiable fluid simulation codes
  * to provide easy access to advancecd algorithms for computing invariant solutions of Navier-Stokes
  
## News

**2018-06-21** Feature freeze in preparation for public release of MPI-parallel channelflow-2.0.
