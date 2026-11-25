"""Linear stability of stratified plane Couette flow, in the settling module's conventions.

Base flow U = y (walls at y = -1, 1 moving at -1, 1), gravity along g = y or z, density perturbation c
about the background -g, buoyancy -Ri c e_g, diffusivity 1/(Re Sc):

    lam u = -i kx y u - v - i kx p + Lap u / Re
    lam v = -i kx y v - Dp + Lap v / Re          (- Ri c if g = y)
    lam w = -i kx y w - i kz p + Lap w / Re      (- Ri c if g = z)
        0 = i kx u + Dv + i kz w
    lam c = -i kx y c + u_g + Lap c / (Re Sc)

no slip; c = 0 (diri) or Dc = 0 (neum) at the walls; Sc = inf drops diffusion and the c BC.
Facchini et al. (2018) have U = -y (x -> -x), Ri = 1/Fr^2, and omega = i lam: Im(omega) = Re(lam).

Chebyshev collocation on Gauss-Lobatto points, generalised eigenproblem A q = lam B q.
Usage: python3 linstab.py Re Ri Sc kx kz g bc [N]
"""
import sys
import numpy as np
import scipy.linalg as sla


def cheb(N):
    # Trefethen's differentiation matrix on x_j = cos(pi j/N), j = 0..N
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2, np.ones(N - 1), 2]) * (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1 / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


def eigs(Re, Ri, Sc, kx, kz, g='z', bc='diri', N=80):
    D, y = cheb(N)
    D2 = D @ D
    n = N + 1
    I = np.eye(n)
    Z = np.zeros((n, n))
    k2 = kx ** 2 + kz ** 2
    Lap = D2 - k2 * I
    adv = -1j * kx * np.diag(y)
    inv_pe = 0.0 if np.isinf(Sc) else 1.0 / (Re * Sc)
    # unknowns q = [u, v, w, c, p]
    A = np.zeros((5 * n, 5 * n), dtype=complex)
    B = np.zeros((5 * n, 5 * n), dtype=complex)
    blk = lambda i, j: (slice(i * n, (i + 1) * n), slice(j * n, (j + 1) * n))
    gi = 1 if g == 'y' else 2
    # momentum
    A[blk(0, 0)] = adv + Lap / Re
    A[blk(0, 1)] = -I                      # -U' v
    A[blk(0, 4)] = -1j * kx * I
    A[blk(1, 1)] = adv + Lap / Re
    A[blk(1, 4)] = -D
    A[blk(2, 2)] = adv + Lap / Re
    A[blk(2, 4)] = -1j * kz * I
    A[blk(gi, 3)] = -Ri * I
    for i in range(3):
        B[blk(i, i)] = I
    # continuity
    A[blk(3, 0)] = 1j * kx * I
    A[blk(3, 1)] = D
    A[blk(3, 2)] = 1j * kz * I
    # scalar (equation rows 4)
    A[blk(4, 3)] = adv + inv_pe * Lap
    A[blk(4, gi)] = I
    B[blk(4, 3)] = I
    # Rows: momentum 0..2, continuity 3, scalar 4. Velocity BCs replace the momentum rows at the walls.
    for i in range(3):
        for j in (0, N):
            r = i * n + j
            A[r, :] = 0
            B[r, :] = 0
            A[r, i * n + j] = 1
    if inv_pe > 0:
        for j in (0, N):
            r = 4 * n + j
            A[r, :] = 0
            B[r, :] = 0
            if bc == 'diri':
                A[r, 3 * n + j] = 1
            else:
                A[r, 3 * n:4 * n] = D[j]
    lam = sla.eig(A, B, right=False)
    lam = lam[np.isfinite(lam)]
    lam = lam[np.abs(lam) < 1e3]
    return lam[np.argsort(-lam.real)]


if __name__ == '__main__':
    Re, Ri, Sc, kx, kz = (float(a) for a in sys.argv[1:6])
    g, bc = sys.argv[6], sys.argv[7]
    N = int(sys.argv[8]) if len(sys.argv) > 8 else 80
    for NN in (N, int(1.5 * N)):
        lam = eigs(Re, Ri, Sc, kx, kz, g, bc, NN)
        print('N =', NN, ' leading eigenvalues:', ' '.join('%.8f%+.8fi' % (l.real, l.imag) for l in lam[:4]))
