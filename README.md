# Channelflow Project

Channelflow is a software system for numerical analysis of the incompressible fluid flow in
channel geometries, written in C++.
It consists of two independent libraries:
* **chflow** that integrates Navier-Stokes equations using semi-implicit finite differences in time and spectral
  discretization in space (Fourier x Chebyshev x Fourier).
* **nsolver** to compute invariant solutions (equilibria, traveling waves, and periodic orbits).

Multiprocessor parallelism is implemented via MPI with two dimensional data distribution.
Parallel input/output is available via NetCDF (requires clustered file system hardware).

The main goals of Channelflow are
  * to lower the barrier to entry to numerical research in fluids
  * to enable creation of short, readable, easily-modifiable fluid simulation codes
  * to provide easy access to advanced algorithms for computing invariant solutions of Navier-Stokes

  
## Getting Started
To obtain the source code and compile it follow the [Installation](./INSTALL.md) procedure.


## Features overview
The Channelflow package provides five programs to perform direct numerical simulations and compute invariant solutions.
For detailed information about functionalities and specific options refer to the helpmode of each program running
`<programname> -h (--help)`.

**1. simulateflow**

   Simulateflow performs DNS of several different incompressible flow systems in channel geometry, from the classic plane Couette
   (PCF) and plane Poiseuille (PPF) to Asymptotic boundary layers (ASBL) and mixed systems. This can be done specifying
   the boundary conditions (velocity at the wall, suction velocity, streamwise direction), the forcing strategy (constant
   pressure gradient or fix flux) and eventually including external forces like rotation.
   The numerics depends on the choice of the nonlinear term, time stepping algorithms, variable/fixed time step and dealiasing.


**2. edgetracking**

   Edgetracking follows the edge of chaos between the laminar and the turbulent basins of attraction using bisection in the cross
   flow energy (Ecf) or in the L2norm of the flow field.


**3. findsoln**

   Findsoln finds invariant solutions of Navier-Stokes equations (equilibria, traveling waves and periodic orbits) using
   Newton-Krylov algorithms.
   Newton iterations can be performed using matrix free techniques (gmres, flexiblegmres, bigcstab) or with the full
   Jacobian (eigen) and Newton updates can be optimized specifying linear or hookstep.
   Algorithmic variants are Tuckerman's Stokes preconditioning algorithm  method and multi-shooting.


**4. continuesoln**

   Continuesoln does a parametric (default) or arclength (-al) continuation of an invariant solution based on
   predictor-corrector method (quadratic extrapolation-Newton search).

**5. findeigenvals**

   Findeigenvals computes spectrum of eigenvalues and the most unstable eigenmodes of equilibria, traveling waves, or periodic orbit
   using Arnoldi iteration (perturbative).



### Tools
The tools directory includes the source files of the pre- and post-processing utilities (listed in the table below)
to manipulate, analyse, compare and convert flowfields.
For detailed information about functionalities and specific options refer to the helpmode of each program running
`<programname> -h (--help)`.

| &nbsp;|&nbsp;|
|:---|:---|
|**addfields** | creates a linear combination of flowfields or adds a base flow (from file or creates it specifying DNSflags) to a flowfield |
|**changegrid**| interpolates a given flowfield ont a different grid with specified geometry |
|**diffops**| Applies a differential operation to a given FlowField <br/> (derivatives; gradient; laplacian; curl; divergence; Q criterion; energy operator; pointwise norm; stramwise average)|
|**extrapolatefield**| Quadratic extrapolation of a FlowField u(mu) from three given flowfields and mu parameters|
|**fieldconvert**| Converts flowfield formats in NetCDF, hdf5, binary (in both directions) and from NetCDF, hdf5, binary to vtk or asci.|
|**fieldprops**| Prints informations about a given FlowField <br/> ( geometry; norms; symmetries; mean constraints; mean, spectral, dynamical and wall properties; statistics; mean velocity profile)|
|**findsymmetries**| Finds the symmetries satisfied by a given FlowField|
|**L2op**| Computes the L2 distance or the inner product between two given FlowFields|
|**optphaseshift**| Computes the optimal phase shift between the first input FlowField and all the other input FlowFields|
|**perturbfield**| Random perturbation of a  given FlowField with another field u with zero divergence and Dirichlet BCs|
|**pressure**| Computes the pressure field of a given FlowField with given DNSFlags|
|**randomfield**| Constructs a random field with zero divergence and Dirichlet BCs |
|**symmetrize**| Translates a FlowField to maximize (or minimize) its shift-reflect and/or shift-rotate symmetry|
|**symmetryop**| Apply a symmetry operation to a given FlowField|



## Troubleshooting
Discourse (Coming soon)

## Bugs Report
The Bug report shall be carried out using the [Issues] feature of GitHub.

## Authors
(Coming soon)

## Licence
[Licence](./LICENCE)

