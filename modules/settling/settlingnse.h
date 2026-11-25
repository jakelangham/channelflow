/**
 * Navier-Stokes equations coupled to a scalar: a dilute settling suspension (Langham & Hogg 2026)
 *
 *   du/dt + u.grad u = -grad p + nu lap u - sum_i Ri_i c_i e_y,   div u = 0,
 *   dc_i/dt + (u - vs_i e_y).grad c_i = kappa_i lap c_i,
 *   no slip for u; no flux kappa_i dc_i/dy + vs_i c_i = 0 at y = -1, 1; mean(c_i) = 1,
 *
 * or a linear stable stratification along g = y (wall-normal) or z (spanwise) (DNSFlags::scalarmodel)
 *
 *   du/dt + u.grad u = -grad p + nu lap u - Ri c e_g,   div u = 0,
 *   dc/dt + u.grad c - u_g = (nu/Sc) lap c,
 *   no slip for u; c = 0 (Dirichlet) or dc/dy = 0 (Neumann) at y = -1, 1,
 *
 * where c is the density perturbation about the background -g (Ri = N^2 h^2/U^2 = 1/Fr^2), e.g.
 * Facchini et al. (2018, JFM 853) for g = z with x -> -x (their U = -y). Fields {u, c_0, ..., c_{n-1}, p}.
 * A settling concentration is likewise stored relative to the background -y (the stored field is
 * c + y), the convention of the existing 4-component flowfield files.
 *
 * The scalars are solved first in each Fourier mode (implicit diffusion, with settling, advection
 * and the background term explicit), then the velocity with the upstream TauSolver and the
 * buoyancy -sum_i Ri_i c_i^{n+1} added to the g-forcing (semi-implicit buoyancy).
 *
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#ifndef SETTLINGNSE_H
#define SETTLINGNSE_H

#include <functional>
#include <vector>
#include "modules/settling/boundarycond.h"
#include "modules/settling/scalarsolver.h"
#include "channelflow/diffops.h"
#include "channelflow/dnsflags.h"
#include "channelflow/flowfield.h"
#include "channelflow/nse.h"
#include "channelflow/tausolver.h"

namespace chflow {

/** Parameters of one scalar: a settling suspension or a linear stratification (DNSFlags::scalarmodel) */
struct SettlingScalar {
    Real Ri;     // buoyancy coefficient (bulk Richardson number)
    Real vs;     // settling velocity, in -y (0 for a stratification)
    Real kappa;  // diffusivity (nu/Sc for a stratification)
    ScalarModel model;
    ScalarWallBC wallbc;  // stratification only

    SettlingScalar(Real Ri_, Real vs_, Real kappa_, ScalarModel model_ = SettlingModel,
                   ScalarWallBC wallbc_ = DirichletWalls)
        : Ri(Ri_), vs(vs_), kappa(kappa_), model(model_), wallbc(wallbc_) {}

    bool stratified() const { return model != SettlingModel; }
    // velocity component along gravity: the buoyancy acts on it and it advects the background
    int gravity() const { return model == StratifiedZ ? 2 : 1; }
    // a stratification's density perturbation changes sign with the reflection of the gravity axis (as
    // u_g does); a settling concentration does not change sign (-1: no axis)
    int parityAxis() const { return stratified() ? gravity() : -1; }

    // No-flux BC kappa c' + vs c = 0 for the stored field c + y: c' + alpha c = 1 + alpha y,
    // alpha = vs/kappa, i.e. g = 1 -/+ alpha at y = -1, 1 for the mean mode, 0 for the others.
    // A stratification's perturbation has homogeneous BCs in every mode.
    Real alpha() const { return vs / kappa; }
    BoundaryCond bc(int Mx, int Mz) const {
        if (!stratified())
            return BoundaryCond(ScalarBCType::NoFlux, Mx, Mz, 1.0 - alpha(), 1.0 + alpha(), alpha());
        return BoundaryCond(wallbc == NeumannWalls ? ScalarBCType::Neum : ScalarBCType::Diri, Mx, Mz, 0.0, 0.0);
    }
    // the BC values of the mean mode (kx = kz = 0)
    Real meanga() const { return stratified() ? 0.0 : 1.0 - alpha(); }
    Real meangb() const { return stratified() ? 0.0 : 1.0 + alpha(); }
    // the settling scalar's mean is fixed at 1 (its mass)
    bool normaliseMass() const { return !stratified(); }
    // with Neumann walls and gravity along y, c + const is also a solution (the constant is taken up by
    // the hydrostatic pressure) and the mean of c is conserved: Newton fixes it at 0. Along z the
    // mean is not neutral: -Ri mean(c) e_z drives a mean spanwise flow, which advects the background.
    bool pinMean() const { return model == StratifiedY && wallbc == NeumannWalls; }
};

/** The scalars defined by DNSFlags (currently one: Ri, vs, kappa, or for a stratification Ri, nu/Sc) */
std::vector<SettlingScalar> settlingScalars(const DNSFlags& flags);

/** Split a 4-component field (velocity + scalar in component 3) and recombine it */
void splitVelocityScalar(const FlowField& u4, FlowField& u, FlowField& c);
void combineVelocityScalar(const FlowField& u, const FlowField& c, FlowField& u4);

/** fc = (u - vs e_y).grad c - (u_g - vs): advection and settling of the stored scalar, including the
 * background -g (g = 1 for a settling scalar, stored as c + y; vs = 0 for a stratification). u must
 * include the base flow. */
void scalarAdvection(const FlowField& u, const FlowField& c, Real vs, int g, FlowField& fc, FlowField& tmp);

/** The symmetries usym as they act on a 1-component scalar that changes sign with axis (as the velocity
 * component along it): s times the sign of that axis; axis < 0 leaves them as they are. */
cfarray<FieldSymmetry> scalarSymmetries(const cfarray<FieldSymmetry>& usym, int axis);

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

    // Prescribed velocity (Ri = 0 only): step only the scalars; the velocity and pressure fields are
    // left as they are, for the caller to set to u(t) before each step. Multistep schemes only.
    void setPrescribedVelocity(bool prescribed);
    bool prescribedVelocity() const { return prescribedVelocity_; }

    // With a prescribed velocity, its RHS is not assembled (NSE::prescribed). Optionally a provider sets
    // it at the start of each step (NSE::prestep): provider(u) must store u(t) of the step about to be
    // taken in u (spectral, on this grid), so the caller can advance several steps per DNS::advance.
    void setVelocityProvider(std::function<void(FlowField&)> provider) { provider_ = std::move(provider); }
    bool prescribed(int i) const override { return prescribedVelocity_ && i == 0; }
    void prestep(std::vector<FlowField>& fields) override {
        if (prescribedVelocity_ && provider_)
            provider_(fields[0]);
    }

   protected:
    std::vector<SettlingScalar> scalars_;
    // scalarsolver_[i][j][mx][mz]: scalar i, substep j, Fourier mode (mx, mz)
    std::vector<ScalarSolver***> scalarsolver_;
    std::vector<ComplexChebyCoeff> ck_;  // scalar profiles of the current mode
    ComplexChebyCoeff Rck_;              // scalar right-hand side of the current mode

    FlowField utot_;  // work field: total velocity u + Ubase, physical in rotationalNonlinear
    bool prescribedVelocity_ = false;
    std::function<void(FlowField&)> provider_;

    void rotationalNonlinear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields);
    void genericNonlinear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields);
    void freeSolvers();
};

}  // namespace chflow
#endif  // SETTLINGNSE_H
