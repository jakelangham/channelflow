# Regression cases

Reference solutions for checking that code changes leave the physics unchanged.
All cases: Re = 400, vs = kappa = 2.5e-3, Lx = 5.51157, Lz = 2.51327 (0.8 pi), y in [-1, 1].
Run serially (1 thread).

    ./run_regression.sh                 # quick tier against ../../build: ~40 s, ~42 MB
    ./run_regression.sh -a              # all tiers (quick + full)
    ./run_regression.sh -b OTHER_BUILD  # against another build
    ./run_regression.sh eq3 settle1     # selected cases, any tier

Tiers: each case script sets `tier=quick` (run after every change) or `tier=full` (slower;
before merges, restructure milestones and production runs). The pseudo-case `unit` (quick)
runs the CTest unit/integration tests that currently pass.

| case | tier | what it checks | cost |
|---|---|---|---|
| unit | quick | CTest: unitTests, tridiag, cheby, helmholtz, poisson, io | <1 s |
| eq1 | quick | Ri = 0 equilibrium, Newton | 4 s |
| eq3 | quick | Ri > 0 travelling wave, Newton + shift | 4 s |
| settle1 | quick | exact settling profile stays steady (analytic) | 15 s |
| traj1 | quick | time stepping at Ri > 0 | 14 s |

The full-resolution originals (`eq1.nc`, `eq3_vk1_sn.nc`, `upo_vk0.56_findsoln/`) are kept
out of git.

## eq1 — equilibrium, Ri = 0

- `eq1_16x35x32.nc`, symmetry `s1.asc` = (1, 1, 1, -1, 0.5, 0), i.e. [u,v,-w](x+Lx/2, y, -z).
- Run with `findsoln_eq1.sh`.
- Provenance: regridded with `changegrid -Nx 16 -Nz 32 -vsok 1` from `eq1.nc`
  (48x35x48, written 2023-10-23 by unpublished revision a26bd36, L2Norm(G) ~ 1e-10),
  then reconverged at commit 061b1aa.
- Reference convergence: L2Norm(G) 4.0e-7 -> 2.6e-10 -> 1.7e-13 in 2 Newton steps
  (22 f^T evaluations), 1 min 19 s, peak RSS 108 MB.
- Truncation: x 5e-5, y 5e-8, z 8e-6 (adequate for regression, not for physics).

## eq3 — travelling wave in x, Ri = 0.0025

- `eq3_16x35x32.nc`, shift `eq3_sigma.asc` (ax = 0.0317574 over T = 20),
  symmetry `sztxz.asc` = (1, 1, 1, -1, 0.5, 0.5), i.e. [u,v,-w](x+Lx/2, y, -z+Lz/2).
- Run with `findsoln_eq3.sh` (`-eqb -xrel`).
- Provenance: `eq3_vk1_sn.nc` (48x35x48, written 2024-06-17, no git revision recorded) is a
  travelling wave at a saddle-node in Ri (Ri = 0.0026442). It is sztxz-symmetric after a
  z-translation of -0.0048046 Lz (anti-symmetric fraction 2e-5 -> 1e-17). Shifted, regridded
  with `changegrid -Nx 16 -Nz 32 -vsok 1`, then converged at commit ad16ec4.
- At the saddle-node Ri the coarse-grid solve converged only linearly (3 steps, 9.7e-4 -> 4.3e-6),
  so the reference was taken at Ri = 0.0025, starting from that iterate.
- Reference convergence at Ri = 0.0025: L2Norm(G) 2.4e-5 -> 3.7e-15 in 9 Newton steps
  (188 f^T evaluations; 6 trust-region-limited steps, then quadratic), 12 min 32 s,
  peak RSS 149 MB.
- Truncation: x 4e-4, y 5e-7, z 5e-6 (adequate for regression, not for physics).

## settle1 — exact settling profile, Ri = 0.0025

- `tools/settlingprofile` writes C(y) = A exp(-r y), r = vs/kappa = 1, A = r / sinh(r)
  (mean C = 1), with zero velocity perturbation. This is an exact steady state: advection by the
  base flow vanishes and buoyancy is balanced hydrostatically. Integrated for T = 10, it must
  not move (strict tolerance 1e-12); observed 8e-14 (velocity, absolute) and 2e-14 (scalar).
- Sensitivity: integrating the same profile with vs = 2.0e-3 instead of 2.5e-3 moves the scalar by
  a relative 1.4e-2 in T = 10.
- The scalar is stored relative to the background -y of the current code (u[3] = C + y).

## traj1 — time stepping, Ri = 0.0025

- Starts from `eq1_16x35x32.nc` (steady only at Ri = 0) and integrates to T = 10 at Ri = 0.0025
  with fixed dt = 0.02 (sbdf3, sbdf1 start-up). Reference: `traj1_T10.nc`. Run with
  `simulate_traj1.sh`; about 16 s, 42 MB.
- Over T = 10 the velocity moves 3.7% (relative) from the initial state, the scalar 0.5%.
- Time-stepping error at dt = 0.02, from runs at dt = 0.02, 0.01, 0.005: relative 3.8e-7
  (velocity), 6.4e-8 (scalar). Observed order is 2, not 3: the first-order sbdf1 start-up step
  limits it. (The higher-order start-ups use NSE::linear, which ignores the scalar.)
- Tolerances: strict 1e-10 (unchanged discretisation); loose 4e-6 velocity, 6e-7 scalar
  (10x the time-stepping error) for changes to the time discretisation.

## Effect of the buoyancy tau-correction fix (commit 9ab93c7)

Running these cases against the code before any fixes (d870ccf):

| case | quantity | pre-fix vs fixed |
|---|---|---|
| eq1 (Ri = 0) | G0 | 1.68e-13 vs 1.68e-13 (control, unaffected) |
| eq3 (Ri = 0.0025) | G0 | 4.8e-10 vs 3.9e-15 |
| eq3 | reconverged solution | rel L2 distance 1.1e-7 (velocity), 1.2e-8 (scalar); ax differs by 2.1e-8 |
| traj1 | u(T = 20), before traj1 was shortened to T = 10 | rel L2 distance 3.0e-9 (velocity), 5.1e-11 (scalar) |

The error scales with Ri and with the size of the top Chebyshev coefficients of the scalar,
so it is smaller still at production resolution. Old solutions at Ri > 0 reconverge in about one
Newton step under the fixed code.
