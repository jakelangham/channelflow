# Regression cases

Reference solutions for checking that code changes leave the physics unchanged.
All cases: Re = 400, vs = kappa = 2.5e-3, Lx = 5.51157, Lz = 2.51327 (0.8 pi), y in [-1, 1].
Run serially (1 thread).

The full-resolution originals (`eq1.nc`, `eq3_vk1_sn.nc`) are kept out of git.

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
