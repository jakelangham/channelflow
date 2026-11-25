/**
 * O(N) solver for the settling scalar's Helmholtz problem in Chebyshev-tau form,
 *
 *   nu u''(y) - lambda u(y) = f(y) on [a, b],
 *
 * with the boundary conditions of BoundaryCond (Diri, DiriRobin, NeumRobin, NoFlux, Neum) in rows 0 and 1
 * of the tau system and, for NoFlux with normalise_mass, mean(u) = mass in row 2. It solves the
 * same linear system as DenseHelmholtzSolver (which stays as its reference in the unit tests) but
 * without the dense N x N factorisation.
 *
 * Method: rows n >= 2 of the system (Canuto & Hussaini eqn 5.1.24) couple only u_{n-2}, u_n,
 * u_{n+2}, so the even and the odd coefficients each form a tridiagonal chain. The leading
 * coefficient of each chain (u_0, u_1, and u_2 too when row 2 is the mass constraint) is left
 * free: each chain's solution is a particular solution plus the free coefficients times
 * homogeneous solutions, which do not depend on f and are computed once. The BC rows (and the mass
 * row) then give a 2x2 (or 3x3) system for the free coefficients. The chains are eliminated from
 * the highest mode down (UL), as in upstream's BandedTridiag.
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#ifndef SETTLING_SCALARSOLVER_H
#define SETTLING_SCALARSOLVER_H

#include <vector>
#include "channelflow/bandedtridiag.h"
#include "channelflow/chebyshev.h"
#include "modules/settling/boundarycond.h"

namespace chflow {

class ScalarSolver {
   public:
    ScalarSolver();
    ScalarSolver(int numberModes, BoundaryCond bc, Real a, Real b, Real lambda, Real nu = 1.0,
                 bool normalise_mass = false);

    // g0 and g1 are the right-hand sides of the BC rows, mass that of the mass row (if any)
    void solve(ChebyCoeff& u, const ChebyCoeff& f, Real g0, Real g1, Real mass_loading) const;

    Real lambda() const { return lambda_; }

   private:
    // One parity chain: unknowns x_k = u_{2k + parity}, k = 0..K-1; equations k = nfree..K-1:
    //   lo_[k] x_{k-1} + d_[k] x_k + up_[k] x_{k+1} = r_k   (up_[K-1] = 0)
    // with x_0..x_{nfree-1} free. UL-factorised from k = K-1 down to k = nfree.
    struct Chain {
        int K = 0;
        int nfree = 0;
        std::vector<Real> lo, d, up;  // coefficients, then (d, lo) overwritten by the UL factors
        void factor();
        // x_k for k >= nfree from r_k (k >= nfree) and the free values x_0..x_{nfree-1} (in x)
        void solve(std::vector<Real>& x, const std::vector<Real>& r) const;
    };

    int N_;           // Canuto & Hussaini's N, section 5.1.2
    int nModes_;      // N_+1
    int nEvenModes_;  // number of even modes
    int nOddModes_;   // number of odd modes
    Real a_;
    Real b_;
    Real lambda_;
    Real nu_;
    bool normalise_mass_;
    inline int beta(int n) const { return (n > N_ - 2) ? 0 : 1; }
    inline int c(int n) const { return (n == 0 || n == N_) ? 2 : 1; }

    BandedTridiag Be_;  // RHS Tridiag for even modes
    BandedTridiag Bo_;  // RHS Tridiag for odd modes
    Chain even_, odd_;

    // the constraint rows (BC rows 0 and 1, and the mass row) over all nModes_ coefficients
    std::vector<std::vector<Real>> crow_;
    // free unknowns, as coefficient indices (0, 1 and, with the mass row, 2), their homogeneous
    // solutions over all coefficients, and the LU factors of M = crow_ . homogeneous solutions
    std::vector<int> free_;
    std::vector<std::vector<Real>> hom_;
    std::vector<Real> M_;  // nc x nc, row major, LU with partial pivoting
    std::vector<int> piv_;
};

}  // namespace chflow
#endif  // SETTLING_SCALARSOLVER_H
