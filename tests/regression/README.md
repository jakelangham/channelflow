# Regression cases

Reference solutions for checking that code changes leave the physics unchanged.
All cases: Re = 400, vs = kappa = 2.5e-3, Lx = 5.51157, Lz = 2.51327 (0.8 pi), y in [-1, 1].
Run serially (1 thread).

The full-resolution originals (`eq1.nc`, `eq3_vk1_sn.nc`) are kept out of git.

## eq1 — equilibrium, Ri = 0

- `eq1_16x35x32.nc`, symmetry `s1.asc` = (1, 1, 1, -1, 0.5, 0), i.e. [u,v,-w](x+Lx/2, y, -z).
- Run with `findsoln.sh`.
- Provenance: regridded with `changegrid -Nx 16 -Nz 32 -vsok 1` from `eq1.nc`
  (48x35x48, written 2023-10-23 by unpublished revision a26bd36, L2Norm(G) ~ 1e-10),
  then reconverged at commit 061b1aa.
- Reference convergence: L2Norm(G) 4.0e-7 -> 2.6e-10 -> 1.7e-13 in 2 Newton steps
  (22 f^T evaluations), 1 min 19 s, peak RSS 108 MB.
- Truncation: x 5e-5, y 5e-8, z 8e-6 (adequate for regression, not for physics).
