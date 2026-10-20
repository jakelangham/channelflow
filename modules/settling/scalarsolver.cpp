/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#include "modules/settling/scalarsolver.h"

#include <cmath>
#include "cfbasics/cfbasics.h"

using namespace std;

namespace chflow {

void ScalarSolver::Chain::factor() {
    // UL elimination from the last row up: row k becomes d_k x_k + lo_k x_{k-1} = r'_k with
    //   d_k <- d_k - up_k lo_{k+1} / d_{k+1}
    for (int k = K - 2; k >= nfree; --k)
        d[k] -= up[k] * lo[k + 1] / d[k + 1];
}

void ScalarSolver::Chain::solve(vector<Real>& x, const vector<Real>& r) const {
    // upward sweep: r'_k = r_k - up_k r'_{k+1} / d_{k+1}; then downward: x_k = (r'_k - lo_k x_{k-1}) / d_k
    vector<Real> rr(K, 0.0);
    for (int k = K - 1; k >= nfree; --k)
        rr[k] = r[k] - (k + 1 < K ? up[k] * rr[k + 1] / d[k + 1] : 0.0);
    for (int k = nfree; k < K; ++k)
        x[k] = (rr[k] - lo[k] * x[k - 1]) / d[k];
}

ScalarSolver::ScalarSolver()
    : N_(0), nModes_(0), nEvenModes_(0), nOddModes_(0), a_(0), b_(0), lambda_(0), nu_(0), normalise_mass_(false) {}

ScalarSolver::ScalarSolver(int numberModes, BoundaryCond bc, Real a, Real b, Real lambda, Real nu, bool normalise_mass)
    : N_(numberModes - 1),
      nModes_(numberModes),
      nEvenModes_(N_ / 2 + 1),
      nOddModes_(N_ / 2),
      a_(a),
      b_(b),
      lambda_(lambda),
      nu_(nu),
      // Total mass is only conserved (and the constraint row only built) for no-flux BCs
      normalise_mass_(normalise_mass && bc.type_ == ScalarBCType::NoFlux),
      Be_(nEvenModes_),
      Bo_(nOddModes_) {
    assert(nModes_ % 2 == 1);
    assert(nModes_ > 4);

    const Real nuscaled = nu / square((b_ - a_) / 2);

    // Tau rows n >= 2 (Canuto & Hussaini eqn 5.1.24), as DenseHelmholtzSolver: the left-hand side
    // of each parity chain, and the right-hand side matrices Be_, Bo_ (row 0 of each is unused)
    even_.K = nEvenModes_;
    odd_.K = nOddModes_;
    even_.nfree = normalise_mass_ ? 2 : 1;
    odd_.nfree = 1;
    for (Chain* ch : {&even_, &odd_}) {
        ch->lo.assign(ch->K, 0.0);
        ch->d.assign(ch->K, 0.0);
        ch->up.assign(ch->K, 0.0);
    }
    Be_.diag(0) = 1.0;
    Bo_.diag(0) = 1.0;
    for (int parity = 0; parity < 2; ++parity) {
        Chain& ch = parity == 0 ? even_ : odd_;
        BandedTridiag& B = parity == 0 ? Be_ : Bo_;
        for (int k = 1; k < ch.K; ++k) {
            const int n = 2 * k + parity;
            const int nn = n * n;
            ch.lo[k] = -(c(n - 2) * lambda) / (4 * n * (n - 1));
            ch.d[k] = (nuscaled + (beta(n) * lambda) / (2 * (nn - 1)));
            if (beta(n + 2) && k + 1 < ch.K)
                ch.up[k] = -lambda / (4 * n * (n + 1));

            B.lodiag(k) = ((Real)c(n - 2)) / (4 * n * (n - 1));
            B.diag(k) = -((Real)beta(n)) / (2 * (nn - 1));
            if (beta(n + 2) && k + 1 < ch.K)
                B.updiag(k) = 1.0 / (4 * n * (n + 1));
        }
    }

    // Constraint rows, as rows 0, 1 (and 2) of DenseHelmholtzSolver's matrix (walls at y = -1, 1)
    vector<Real> r0(nModes_, 0.0), r1(nModes_, 0.0);
    for (int i = 0; i < nModes_; ++i) {
        const Real ii = Real(i) * i;
        const Real sgn = (i % 2 == 0) ? 1.0 : -1.0;  // T_i(-1)
        switch (bc.type_) {
            case ScalarBCType::Diri:
                r0[i] = sgn;  // u(-1)
                r1[i] = 1.0;  // u(1)
                break;
            case ScalarBCType::DiriRobin:
                r0[i] = sgn;                // u(-1)
                r1[i] = ii + bc.alpha_;     // u'(1) + alpha u(1)
                break;
            case ScalarBCType::NeumRobin:
                r0[i] = -sgn * ii;          // u'(-1)
                r1[i] = ii + bc.alpha_;     // u'(1) + alpha u(1)
                break;
            case ScalarBCType::NoFlux:
                r0[i] = sgn * (bc.alpha_ - ii);  // u'(-1) + alpha u(-1)
                r1[i] = ii + bc.alpha_;          // u'(1) + alpha u(1)
                break;
        }
    }
    crow_ = {r0, r1};
    free_ = {0, 1};
    if (normalise_mass_) {
        vector<Real> r2(nModes_, 0.0);
        for (int i = 0; i < nModes_; i += 2)
            r2[i] = 1.0 / (1.0 - Real(i) * i);  // mean(u)
        crow_.push_back(r2);
        free_.push_back(2);
    }

    even_.factor();
    odd_.factor();

    // Homogeneous solutions: free unknown j = 1, the other free unknowns 0, chain right-hand sides 0
    const int nc = free_.size();
    hom_.assign(nc, vector<Real>(nModes_, 0.0));
    for (int j = 0; j < nc; ++j) {
        const int parity = free_[j] % 2;
        Chain& ch = parity == 0 ? even_ : odd_;
        vector<Real> x(ch.K, 0.0), r(ch.K, 0.0);
        x[free_[j] / 2] = 1.0;
        ch.solve(x, r);
        for (int k = 0; k < ch.K; ++k)
            hom_[j][2 * k + parity] = x[k];
    }

    // M_{cj} = crow_c . hom_j, LU-factorised with partial pivoting
    M_.assign(nc * nc, 0.0);
    for (int ci = 0; ci < nc; ++ci)
        for (int j = 0; j < nc; ++j) {
            Real s = 0.0;
            for (int i = 0; i < nModes_; ++i)
                s += crow_[ci][i] * hom_[j][i];
            M_[ci * nc + j] = s;
        }
    piv_.assign(nc, 0);
    for (int k = 0; k < nc; ++k) {
        int p = k;
        for (int i = k + 1; i < nc; ++i)
            if (abs(M_[i * nc + k]) > abs(M_[p * nc + k]))
                p = i;
        piv_[k] = p;
        if (p != k)
            for (int j = 0; j < nc; ++j)
                swap(M_[k * nc + j], M_[p * nc + j]);
        if (M_[k * nc + k] == 0.0)
            cferror("ScalarSolver : singular constraint system (check the BC and lambda)");
        for (int i = k + 1; i < nc; ++i) {
            M_[i * nc + k] /= M_[k * nc + k];
            for (int j = k + 1; j < nc; ++j)
                M_[i * nc + j] -= M_[i * nc + k] * M_[k * nc + j];
        }
    }
}

void ScalarSolver::solve(ChebyCoeff& u, const ChebyCoeff& f, Real g0, Real g1, Real mass_loading) const {
    assert(f.state() == Spectral);
    // right-hand sides of the tau rows, as DenseHelmholtzSolver
    Be_.multiplyStrided(f, u, 0, 2);
    Bo_.multiplyStrided(f, u, 1, 2);

    // particular solution: chains with the free unknowns 0
    vector<Real> p(nModes_, 0.0);
    for (int parity = 0; parity < 2; ++parity) {
        const Chain& ch = parity == 0 ? even_ : odd_;
        vector<Real> x(ch.K, 0.0), r(ch.K, 0.0);
        for (int k = ch.nfree; k < ch.K; ++k)
            r[k] = u[2 * k + parity];
        ch.solve(x, r);
        for (int k = 0; k < ch.K; ++k)
            p[2 * k + parity] = x[k];
    }

    // free unknowns from the constraint rows: M c = g - crow . p
    const int nc = free_.size();
    Real g[3] = {g0, g1, mass_loading};
    for (int ci = 0; ci < nc; ++ci)
        for (int i = 0; i < nModes_; ++i)
            g[ci] -= crow_[ci][i] * p[i];
    // all row interchanges first (the factorisation swapped whole rows, multipliers included), then
    // L and U
    for (int k = 0; k < nc; ++k)
        swap(g[k], g[piv_[k]]);
    for (int k = 0; k < nc; ++k)
        for (int i = k + 1; i < nc; ++i)
            g[i] -= M_[i * nc + k] * g[k];
    for (int k = nc - 1; k >= 0; --k) {
        for (int j = k + 1; j < nc; ++j)
            g[k] -= M_[k * nc + j] * g[j];
        g[k] /= M_[k * nc + k];
    }

    for (int i = 0; i < nModes_; ++i) {
        Real ui = p[i];
        for (int j = 0; j < nc; ++j)
            ui += g[j] * hom_[j][i];
        u[i] = ui;
    }
    u.setState(Spectral);
}

}  // namespace chflow
