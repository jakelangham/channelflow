/**
 * Navier-Stokes equations coupled to dilute settling suspensions (Langham & Hogg 2026)
 *
 *   du/dt + u.grad u = -grad p + nu lap u - sum_i Ri_i c_i e_y,   div u = 0,
 *   dc_i/dt + (u - vs_i e_y).grad c_i = kappa_i lap c_i,
 *   no slip for u; no flux kappa_i dc_i/dy + vs_i c_i = 0 at y = -1, 1; mean(c_i) = 1,
 *
 * on fields {u, c_0, ..., c_{n-1}, p}. Each concentration is stored relative to the background -y
 * (the stored field is c + y), the convention of the existing 4-component flowfield files.
 *
 * The scalars are solved first in each Fourier mode (implicit diffusion, with settling, advection
 * and the background term explicit), then the velocity with the upstream TauSolver and the
 * buoyancy -sum_i Ri_i c_i^{n+1} added to the y-forcing (semi-implicit buoyancy).
 *
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#ifndef SETTLINGNSE_H
#define SETTLINGNSE_H

#include <vector>
#include "modules/settling/boundarycond.h"
#include "modules/settling/densehelmholtz.h"
#include "channelflow/diffops.h"
#include "channelflow/dnsflags.h"
#include "channelflow/flowfield.h"
#include "channelflow/nse.h"
#include "channelflow/tausolver.h"

namespace chflow {

/** Parameters of one settling scalar */
struct SettlingScalar {
    Real Ri;     // buoyancy coefficient (bulk Richardson number)
    Real vs;     // settling velocity, in -y
    Real kappa;  // diffusivity

    // No-flux BC kappa c' + vs c = 0 for the stored field c + y: c' + alpha c = 1 + alpha y,
    // alpha = vs/kappa, i.e. g = 1 -/+ alpha at y = -1, 1 for the mean mode, 0 for the others.
    Real alpha() const { return vs / kappa; }
    BoundaryCond bc(int Mx, int Mz) const { return BoundaryCond(ScalarBCType::NoFlux, Mx, Mz, 1.0 - alpha(), 1.0 + alpha(), alpha()); }
};

/** The settling scalars defined by DNSFlags (currently one: Ri, vs, kappa) */
std::vector<SettlingScalar> settlingScalars(const DNSFlags& flags);

/** Split a 4-component field (velocity + scalar in component 3) and recombine it */
void splitVelocityScalar(const FlowField& u4, FlowField& u, FlowField& c);
void combineVelocityScalar(const FlowField& u, const FlowField& c, FlowField& u4);

/** fc = (u - vs e_y).grad c - (v - vs): advection and settling of the stored scalar c + y, including
 * the background -y. u must include the base flow. */
void scalarAdvection(const FlowField& u, const FlowField& c, Real vs, FlowField& fc, FlowField& tmp);

class SettlingNSE : public NSE {
   public:
    SettlingNSE(const std::vector<FlowField>& fields, const DNSFlags& flags);
    SettlingNSE(const SettlingNSE&) = delete;
    SettlingNSE& operator=(const SettlingNSE&) = delete;
    virtual ~SettlingNSE();

    void nonlinear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields) override;
    void linear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields) override;
    void solve(std::vector<FlowField>& outfields, const std::vector<FlowField>& infields, const int i = 0) override;
    void reset_lambda(const std::vector<Real> lambda_t) override;
    std::vector<FlowField> createRHS(const std::vector<FlowField>& fields) const override;
    std::vector<cfarray<FieldSymmetry>> createSymmVec() const override;

    int nscalars() const { return scalars_.size(); }
    const std::vector<SettlingScalar>& scalars() const { return scalars_; }

   protected:
    std::vector<SettlingScalar> scalars_;
    // scalarsolver_[i][j][mx][mz]: scalar i, substep j, Fourier mode (mx, mz)
    std::vector<DenseHelmholtzSolver***> scalarsolver_;
    std::vector<ComplexChebyCoeff> ck_;  // scalar profiles of the current mode
    ComplexChebyCoeff Rck_;              // scalar right-hand side of the current mode

    void freeSolvers();
};

}  // namespace chflow
#endif  // SETTLINGNSE_H
