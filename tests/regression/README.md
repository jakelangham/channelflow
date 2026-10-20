# Regression cases

Reference solutions for checking that code changes leave the physics unchanged.
All cases: Re = 400, vs = kappa = 2.5e-3, Lx = 5.51157, Lz = 2.51327 (0.8 pi), y in [-1, 1].
Run serially (1 thread).

    ./run_regression.sh                 # quick tier against ../../build: ~40 s, ~42 MB
    ./run_regression.sh -a              # all tiers (quick + full)
    ./run_regression.sh -a -m MPI_BUILD # ... including the MPI case
    ./run_regression.sh -b OTHER_BUILD  # against another build (from 2026-10-18 on: the case
                                        # scripts use -vsett, which older builds do not accept)
    ./run_regression.sh eq3 settle1     # selected cases, any tier

Tiers: each case script sets `tier=quick` (run after every change) or `tier=full` (slower;
before merges, restructure milestones and production runs). The pseudo-case `unit` (quick)
runs the CTest unit/integration tests that currently pass.

| case | tier | what it checks | cost |
|---|---|---|---|
| unit | quick | CTest: unitTests, tridiag, cheby, laminar, helmholtz, tausolver, poisson, io, vector2field | 1 s |
| eq1 | quick | Ri = 0 equilibrium, Newton | 4 s |
| eq3 | quick | Ri > 0 travelling wave, Newton + shift | 4 s |
| settle1 | quick | exact settling profile stays steady (analytic), vs/kappa = 1 | 15 s |
| settle2 | quick | the same at vs/kappa = 14 (bottom boundary layer), Ny = 49 | 15 s |
| decay1 | quick | exact decaying scalar mode (m = 1, n = 1, vs/kappa = 1) | 14 s |
| traj1 | quick | time stepping at Ri > 0 | 14 s |
| rotation, nonlinearity | quick | unsupported options with the scalar stop with an error | 1 s |
| vsett1 | full | continuation in vs keeps the no-flux BC consistent | 8 min |
| kappa1 | full | continuation in kappa (no-flux BC and diffusion follow) | 8 min |
| ri1 | full | continuation in Ri from 0 (eq1 becomes a travelling wave, -xrel) | 15 min |
| decay2 | full | exact decaying scalar mode (m = 2, n = 2, vs/kappa = 4) | 14 s |
| steady | full | exact steady profile under every scheme and both constraints (7 XFAIL) | 45 s |
| onset1 | full | stability onset of the laminar state vs Langham & Hogg (2026) App. B | 10 min |
| upo1 | full | periodic orbit (T ~ 87.8), Newton + period | 1 min |
| mpi | full | eq1, eq3, traj1 on 2 MPI processes, grids 2x1 and 1x2, match serial (needs `-m`) | 30 s |

Every case also checks that its output field keeps the NetCDF layout of the reference
(dimensions, `Component_0`..`Component_3`, geometry attributes), since external scripts read it.

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

## settle2 — exact settling profile at vs/kappa = 14

- As settle1, with vs = 0.035 (vs/kappa = 14, the largest value in Langham & Hogg 2026) and Ny = 49,
  which resolves the exponential profile to ~4e-25 (its Chebyshev coefficients scale like I_n(14)).
  Exercises the no-flux BC rows at large alpha = 14 and explicit settling near its time-step limit
  (vs dt / dy_min ~ 0.33 at dt = 0.02). Observed: moves by 3e-14 (velocity) and 2.6e-13 (scalar).

## decay1, decay2 — exact decaying scalar modes, Ri = 0

- With zero velocity perturbation and Ri = 0, a spanwise scalar perturbation
  eps phi_m(y) cos(2 pi n z / Lz) exp(sigma t) on the base profile is an exact solution (the Couette
  base flow does not advect z-dependence). Substituting phi = exp(-r y/2) psi, r = vs/kappa, the
  no-flux problem becomes psi'' = -mu^2 psi, psi' + (r/2) psi = 0 at y = +-1, whose eigenvalues
  satisfy sin(mu) cos(mu) = 0: mu = m pi/2 and sigma = -kappa (k^2 + r^2/4 + mu^2) (m >= 1), or
  phi = exp(-r y), sigma = -kappa k^2 (m = 0). `tools/settlingprofile -eps -m -n -t` writes it.
- The only tests of the non-mean scalar modes and the Robin BC against an exact solution (settle1/2
  test the mean mode). Errors are second order in dt: decay1 2.9e-8 at dt = 0.02, 7.2e-9 at 0.01;
  decay2 (faster decay, x0.38 over T = 10) 2.2e-7. Tolerances are 10x these. The m = 0 mode was also
  checked (6.7e-9) but is not a case.

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

## vsett1 — continuation in the settling velocity

- `continue_vsett1.sh` continues eq1 (Ri = 0) in `vsett` for two steps (vs = 2.5e-3 to 3.1e-3),
  then evaluates the last solution with a fresh `findsoln -Nn 0` at that vs, which builds the
  no-flux BC from scratch. Residual must be below 1e-11; observed 2e-15.
- Sensitivity: the same evaluation at the previous step's vs gives ~1e-3, i.e. a frozen BC would
  be detected. (Before commit 2179c53, cfDSI kept the starting vs/kappa in the BC throughout.)
- At Ri = 0 the velocity does not depend on vs, so only the scalar changes along this branch.

## steady — every time-stepping scheme and both mean constraints

- From upstream's dnsParabolaTest/dnsZeroTest: the exact steady settling profile (vs/kappa = 1,
  Ri = 0.0025, Couette flow) must stay put for T = 2 under each scheme:startup:constraint
  combination listed in `schemes_steady.sh`.
- Passing (1e-14 to 3e-14): CNFE1, SBDF1-4 with first-order start-up, for -mc gradp.
- Expected failures (XFAIL), recorded in the script's xfail list and reported as XPASS once fixed:
  CNAB2, CNRK2, SMRK2, and the SMRK2/CNRK2 start-ups, use NSE::linear, which ignores the scalar
  (the scalar drifts by ~0.1, or ~3e-3 from the start-up alone; the velocity stays exact);
  -mc bulkv stops with an error (not supported with the scalar yet).

## rotation, nonlinearity — unsupported options

- `error_rotation.sh` (-rot) and `error_nonlinearity.sh` (-nl conv) must stop with a message rather
  than run: rotation with the scalar was silently wrong, and only the rotational form of the
  nonlinear term handles the scalar.

## onset1 — linear-stability onset against the paper

- Langham & Hogg (2026), Appendix B: at Re = 400, Sc = 1, vs/kappa = 1, in this box, the
  streamwise-invariant (0, 1) mode of the laminar state becomes unstable at Ri_b = -9.08e-4
  (computed there with Chebfun, independently of channelflow). `eigen_onset1.sh` runs
  findeigenvals (T = 10, 40 Arnoldi iterations, 16x35x32) on the laminar state at Ri_b = -8.5e-4
  and -9.5e-4 and interpolates where the leading growth rate crosses zero.
- The leading mode is a degenerate pair (cos and sin in z). Arnoldi converges one copy to a residual
  of ~1e-7 and the other only to ~1e-3, and the poorly converged copy varies between runs, so the case
  takes the largest growth rate among eigenvalues with residual < 1e-5.
- Observed: sigma = -7.090e-4 and +5.014e-4, onset -9.0858e-4, 0.063% from the paper (tolerance
  0.2%). Checks the coupled velocity-scalar linear operator, buoyancy, the scalar BCs and the
  linearisation/Arnoldi path against an external result.

## upo1 — periodic orbit, Ri = 0, vs/kappa = 0.56

- `upo1_24x35x40.nc`, period T = 87.84095659533736, symmetry `sztxz.asc`, box Lx = 1.75 pi,
  Lz = 1.2 pi (different from the other cases). vs = 1.4e-3, kappa = 2.5e-3. Run with
  `findsoln_upo1.sh` (`-orb`); about 50 s per f^T.
- Provenance: the user's `upo_vk0.56_findsoln/ubest.nc` (60x61x60, T = 87.8902, L2Norm(G) = 9.9e-12,
  written 2026-06-23 at ef3d6d6). 16x35x32 was too coarse (initial L2Norm(G) = 1.8e-2). At 32x35x48
  it converged to 2.0e-11 in 5 Newton steps (T = 87.91901, ~95 s per f^T). Regridding that orbit to
  24x35x40 and reconverging: 1.7e-3 -> 2.2e-4 -> 1.8e-6 -> 6.0e-9 -> 2.5e-10 -> 1.2e-11 in 5 steps
  (1 h, 157 MB), T = 87.84096 (0.09% below the 32x35x48 period), ~46 s per f^T.
- Truncation at 24x35x40: x 3e-4, y 3e-7, z 2e-4.

## mpi — parallel decomposition

- Re-runs eq1, eq3 and traj1 with `mpirun -np 2`, on process grids 2x1 and 1x2, using an MPI build (`-m DIR` or `MPI_BUILD`), with
  the same references and tolerances as the serial cases. Observed: residuals agree with serial to
  the last digits (e.g. 1.67888e-13 vs 1.67907e-13), distances <= 2e-13.
- `RR_MPI_FINDSOLN` / `RR_MPI_SIMULATEFLOW` override the MPI build's programs (as `RR_FINDSOLN` etc.
  do for the serial cases).

## kappa1, ri1 — continuation in kappa and Ri

- Built like vsett1 (`continue_kappa1.sh`, `continue_ri1.sh`): two continuation steps from eq1, then
  a fresh `findsoln -Nn 0` at the last parameter value, read at full precision from the search
  directory's mu.asc (MuD.asc rounds to 8 digits, which alone gave a residual of 4e-11).
- kappa1: kappa = 2.5e-3 to 2.97e-3, fresh residual 1.1e-15.
- ri1: Ri = 0 to 6.28e-4 with -xrel: eq1 becomes a travelling wave as Ri moves away from 0 (the
  shift-rotate symmetry is broken; Langham & Hogg 2026, section 4); fresh residual 6.6e-15. Ri
  continuation is how most of the paper's results were obtained.

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
