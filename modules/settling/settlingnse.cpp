/**
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include "modules/settling/settlingnse.h"

using namespace std;

namespace chflow {

vector<SettlingScalar> settlingScalars(const DNSFlags& flags) { return {SettlingScalar{flags.Ri, flags.vs, flags.kappa}}; }

void splitVelocityScalar(const FlowField& u4, FlowField& u, FlowField& c) {
    if (u4.Nd() != 4)
        cferror("splitVelocityScalar : expected a 4-component field (velocity and scalar), got " +
                i2s(u4.Nd()) + " components");
    u = FlowField(u4.Nx(), u4.Ny(), u4.Nz(), 3, u4.Lx(), u4.Lz(), u4.a(), u4.b(), u4.cfmpi(), u4.xzstate(),
                  u4.ystate());
    c = FlowField(u4.Nx(), u4.Ny(), u4.Nz(), 1, u4.Lx(), u4.Lz(), u4.a(), u4.b(), u4.cfmpi(), u4.xzstate(),
                  u4.ystate());
    u.setPadded(u4.padded());
    c.setPadded(u4.padded());
    u.copySubfields(u4, {0, 1, 2}, {0, 1, 2});
    c.copySubfields(u4, {3}, {0});
}

void combineVelocityScalar(const FlowField& u, const FlowField& c, FlowField& u4) {
    assert(u.Nd() == 3 && c.Nd() == 1);
    if (!u4.geomCongruent(u) || u4.Nd() != 4)
        u4.resize(u.Nx(), u.Ny(), u.Nz(), 4, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    u4.setPadded(u.padded());
    u4.copySubfields(u, {0, 1, 2}, {0, 1, 2});
    u4.copySubfields(c, {0}, {3});
}

void scalarAdvection(const FlowField& u_, const FlowField& c_, Real vs, FlowField& fc, FlowField& tmp) {
    FlowField& u = const_cast<FlowField&>(u_);
    FlowField& c = const_cast<FlowField&>(c_);
    FlowField& gradc = tmp;
    assert(u.Nd() >= 3 && c.Nd() == 1);

    if (!c.geomCongruent(fc) || fc.Nd() != 1)
        fc.resize(c.Nx(), c.Ny(), c.Nz(), 1, c.Lx(), c.Lz(), c.a(), c.b(), c.cfmpi());
    if (!c.geomCongruent(gradc) || gradc.Nd() < 3)
        gradc.resize(c.Nx(), c.Ny(), c.Nz(), 3, c.Lx(), c.Lz(), c.a(), c.b(), c.cfmpi());

    const fieldstate uxzstate = u.xzstate(), uystate = u.ystate();
    const fieldstate cxzstate = c.xzstate(), cystate = c.ystate();

    c.makeSpectral();
    grad(c, gradc);
    gradc.makePhysical();
    u.makePhysical();
    fc.setState(Physical, Physical);

    const lint Nz = u.Nz();
    const lint nxlocmin = u.nxlocmin();
    const lint nxlocmax = u.nxlocmin() + u.Nxloc();
    const lint nylocmin = u.nylocmin();
    const lint nylocmax = u.nylocmax();
    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (lint nz = 0; nz < Nz; ++nz) {
                Real f = 0.0;
                for (int j = 0; j < 3; ++j)
                    f += u(nx, ny, nz, j) * gradc(nx, ny, nz, j);
                f -= u(nx, ny, nz, 1) - vs;          // advection of the background -y, with settling
                f -= vs * gradc(nx, ny, nz, 1);       // settling of the perturbation
                fc(nx, ny, nz, 0) = f;
            }

    fc.makeSpectral();
    u.makeState(uxzstate, uystate);
    c.makeState(cxzstate, cystate);
}

SettlingNSE::SettlingNSE(const vector<FlowField>& fields, const DNSFlags& flags)
    : NSE(fields, flags), scalars_(settlingScalars(flags)), Rck_(Nyd_, a_, b_, Spectral) {
    if (fields.size() != scalars_.size() + 2)
        cferror("SettlingNSE : expected fields {u, c_0, ..., c_n-1, p} with " + i2s(scalars_.size()) +
                " scalar(s), got " + i2s(fields.size()) + " fields");
    if (fields[0].Nd() != 3)
        cferror("SettlingNSE : the velocity field must have 3 components");
    for (uint i = 0; i < scalars_.size(); ++i)
        if (fields[1 + i].Nd() != 1)
            cferror("SettlingNSE : scalar fields must have 1 component");
    if (Nyd_ % 2 == 0)
        cferror("SettlingNSE : the scalar solver needs an odd number of Chebyshev modes");
    if (flags.nonlinearity == LinearAboutProfile)
        cferror("SettlingNSE : the linearised nonlinearity is not implemented for the settling scalar");
    ck_.assign(scalars_.size(), ComplexChebyCoeff(Nyd_, a_, b_, Spectral));
    scalarsolver_.assign(scalars_.size(), nullptr);
}

SettlingNSE::~SettlingNSE() { freeSolvers(); }

void SettlingNSE::freeSolvers() {
    for (auto& solver : scalarsolver_) {
        if (!solver)
            continue;
        for (uint j = 0; j < lambda_t_.size(); ++j) {
            for (int mx = 0; mx < Mxloc_; ++mx)
                delete[] solver[j][mx];
            delete[] solver[j];
        }
        delete[] solver;
        solver = nullptr;
    }
}

void SettlingNSE::setPrescribedVelocity(bool prescribed) {
    if (prescribed) {
        for (const SettlingScalar& s : scalars_)
            if (s.Ri != 0.0)
                cferror("SettlingNSE : a prescribed velocity needs Ri == 0 (a passive scalar)");
        if (flags_.nonlinearity != Rotational)
            cferror("SettlingNSE : a prescribed velocity needs the rotational nonlinearity");
    }
    prescribedVelocity_ = prescribed;
}

// Project a spectral field onto real fields: the kz = 0 modes kx and -kx become complex conjugates and
// the kx = 0 one real, as a transform to physical space and back would leave them. Returns false,
// changing nothing, if this process does not hold every kx (x distributed over MPI processes); that is
// the same on every process, so the caller's fallback stays collective.
static bool makeRealField(FlowField& f) {
    assert(f.xzstate() == Spectral);
    if (f.Mxloc() != f.Mx())
        return false;
    if (f.mzlocmin() != 0)
        return true;  // no kz = 0 modes on this process
    for (int mx = 0; mx < f.Mx(); ++mx) {
        const int kx = f.kx(mx);
        if (kx < 0)
            continue;
        const int mxm = f.mx(-kx);
        const bool self = (kx == 0 || mxm == mx);  // kx = 0, or the Nyquist mode of an even Nx
        for (int i = 0; i < f.Nd(); ++i)
            for (int my = 0; my < f.My(); ++my) {
                if (self) {
                    f.cmplx(mx, my, 0, i) = Complex(Re(f.cmplx(mx, my, 0, i)), 0.0);
                } else {
                    const Complex a = 0.5 * (f.cmplx(mx, my, 0, i) + conj(f.cmplx(mxm, my, 0, i)));
                    f.cmplx(mx, my, 0, i) = a;
                    f.cmplx(mxm, my, 0, i) = conj(a);
                }
            }
    }
    return true;
}

void SettlingNSE::nonlinear(const vector<FlowField>& infields, vector<FlowField>& outfields) {
    // Keep the fields real. The solve updates every Fourier mode on its own, so the kz = 0 modes kx and
    // -kx can drift by round-off from the conjugate symmetry of a real field. That part is invisible in
    // physical space, but at high vs/kappa the scalar's grows: from t ~ 300 at vs/kappa = 10 (Ri = 0,
    // upo1 and its 60x61x60 parent; not at 0.56). The upstream nonlinear terms, and the original settling
    // code's densityAdvection, transformed the fields in place each step, which removes it; projecting
    // directly does the same far more cheaply (the round trip cost 38% at 24x35x40), and the round trip
    // remains the fallback with x distributed over MPI processes. With a prescribed velocity, the caller
    // supplies a velocity that has been through this.
    for (int i = prescribedVelocity_ ? 1 : 0; i <= nscalars(); ++i) {
        FlowField& f = const_cast<FlowField&>(infields[i]);
        if (!makeRealField(f)) {
            f.makePhysical();
            f.makeSpectral();
        }
    }

    if (flags_.nonlinearity == Rotational)
        rotationalNonlinear(infields, outfields);
    else
        genericNonlinear(infields, outfields);

    if (flags_.dealias_xz())
        for (int i = 0; i <= nscalars(); ++i)
            outfields[i].zeroPaddedModes();
}

// Rotational form in one pass: the total velocity u + Ubase (with suction) goes to physical space once,
// for both (curl u) x u and the scalars' advection, and is never transformed back (it is a copy).
// The same arithmetic as navierstokesNL + scalarAdvection, with 13 instead of 22 transforms per scalar.
void SettlingNSE::rotationalNonlinear(const vector<FlowField>& infields, vector<FlowField>& outfields) {
    const FlowField& u0 = infields[0];
    assert(u0.xzstate() == Spectral && u0.ystate() == Spectral);
    utot_ = u0;
    if (utot_.taskid() == utot_.task_coeff(0, 0)) {
        for (int ny = 0; ny < utot_.Ny(); ++ny) {
            utot_.cmplx(0, ny, 0, 0) += Complex(Ubase_(ny), 0.0);
            utot_.cmplx(0, ny, 0, 2) += Complex(Wbase_(ny), 0.0);
        }
        utot_.cmplx(0, 0, 0, 1) -= Complex(flags_.Vsuck, 0.);
    }

    const lint Nz = utot_.Nz();
    const lint nxlocmin = utot_.nxlocmin();
    const lint nxlocmax = utot_.nxlocmin() + utot_.Nxloc();
    const lint nylocmin = utot_.nylocmin();
    const lint nylocmax = utot_.nylocmax();

    // velocity: f = (curl u) x u, minus the rotation term, as rotationalNL and navierstokesNL
    FlowField& f = outfields[0];
    if (prescribedVelocity_) {
        utot_.makePhysical();
        f.setToZero();
    } else {
        FlowField& vort = tmp_;
        curl(utot_, vort);
        utot_.makePhysical();
        vort.makePhysical();
        cross(vort, utot_, f, Physical);
    }
    if (flags_.rotation != 0.0 && !prescribedVelocity_)
        for (lint ny = nylocmin; ny < nylocmax; ++ny)
            for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                for (lint nz = 0; nz < Nz; ++nz) {
                    f(nx, ny, nz, 0) -= flags_.rotation * utot_(nx, ny, nz, 1);
                    f(nx, ny, nz, 1) += flags_.rotation * utot_(nx, ny, nz, 0);
                }
    if (!prescribedVelocity_)
        f.makeSpectral();

    // scalars: as scalarAdvection, with the physical total velocity
    for (int i = 0; i < nscalars(); ++i) {
        FlowField& c = const_cast<FlowField&>(infields[1 + i]);
        FlowField& fc = outfields[1 + i];
        const Real vs = scalars_[i].vs;
        assert(c.xzstate() == Spectral && c.ystate() == Spectral);
        if (!c.geomCongruent(fc) || fc.Nd() != 1)
            fc.resize(c.Nx(), c.Ny(), c.Nz(), 1, c.Lx(), c.Lz(), c.a(), c.b(), c.cfmpi());
        FlowField& gradc = tmp_;
        grad(c, gradc);
        gradc.makePhysical();
        fc.setState(Physical, Physical);
        for (lint ny = nylocmin; ny < nylocmax; ++ny)
            for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
                for (lint nz = 0; nz < Nz; ++nz) {
                    Real fcn = 0.0;
                    for (int j = 0; j < 3; ++j)
                        fcn += utot_(nx, ny, nz, j) * gradc(nx, ny, nz, j);
                    fcn -= utot_(nx, ny, nz, 1) - vs;     // advection of the background -y, with settling
                    fcn -= vs * gradc(nx, ny, nz, 1);      // settling of the perturbation
                    fc(nx, ny, nz, 0) = fcn;
                }
        fc.makeSpectral();
    }
}

// Any nonlinear form: the upstream navierstokesNL for the velocity, then scalarAdvection
void SettlingNSE::genericNonlinear(const vector<FlowField>& infields, vector<FlowField>& outfields) {
    // Velocity: the upstream nonlinear term (any method; base flow and suction added internally)
    navierstokesNL(infields[0], Ubase_, Wbase_, outfields[0], tmp_, flags_);

    // Scalars: advection by the total velocity u + Ubase (with suction), settling, background
    FlowField& u = const_cast<FlowField&>(infields[0]);
    const bool owner = (u.taskid() == u.task_coeff(0, 0));
    for (int ny = 0; ny < u.Ny(); ++ny)
        if (owner) {
            u.cmplx(0, ny, 0, 0) += Complex(Ubase_(ny), 0.0);
            u.cmplx(0, ny, 0, 2) += Complex(Wbase_(ny), 0.0);
        }
    if (owner)
        u.cmplx(0, 0, 0, 1) -= Complex(flags_.Vsuck, 0.);

    for (int i = 0; i < nscalars(); ++i)
        scalarAdvection(u, infields[1 + i], scalars_[i].vs, outfields[1 + i], tmp_);

    for (int ny = 0; ny < u.Ny(); ++ny)
        if (owner) {
            u.cmplx(0, ny, 0, 0) -= Complex(Ubase_(ny), 0.0);
            u.cmplx(0, ny, 0, 2) -= Complex(Wbase_(ny), 0.0);
        }
    if (owner)
        u.cmplx(0, 0, 0, 1) += Complex(flags_.Vsuck, 0.);
}

void SettlingNSE::linear(const vector<FlowField>& infields, vector<FlowField>& outfields) {
    // Linear terms L for fields {u, c_i} (pressure is the last infield):
    //   L u = nu u'' - k^2 nu u - grad p + C,   L c_i = kappa_i (c_i'' - k^2 c_i)
    assert(infields.size() == outfields.size() + 1);
    const FlowField& press = infields.back();
    const int kxmax = infields[0].kxmax();
    const int kzmax = infields[0].kzmax();

    for (lint mx = mxlocmin_; mx < mxlocmin_ + Mxloc_; ++mx) {
        const int kx = infields[0].kx(mx);
        for (lint mz = mzlocmin_; mz < mzlocmin_ + Mzloc_; ++mz) {
            const int kz = infields[0].kz(mz);
            if ((kx == kxmax || kz == kzmax) || (flags_.dealias_xz() && isAliasedMode(kx, kz)))
                break;

            for (int ny = 0; ny < Nyd_; ++ny) {
                uk_.set(ny, flags_.nu * infields[0].cmplx(mx, ny, mz, 0));
                vk_.set(ny, flags_.nu * infields[0].cmplx(mx, ny, mz, 1));
                wk_.set(ny, flags_.nu * infields[0].cmplx(mx, ny, mz, 2));
                Pk_.set(ny, press.cmplx(mx, ny, mz, 0));
            }
            diff2(uk_, Ruk_, Pyk_);
            diff2(vk_, Rvk_, Pyk_);
            diff2(wk_, Rwk_, Pyk_);
            diff(Pk_, Pyk_);

            const Real kappa2 = 4 * pi * pi * (square(kx / Lx_) + square(kz / Lz_));
            const Complex Dx = infields[0].Dx(mx);
            const Complex Dz = infields[0].Dz(mz);
            for (int ny = 0; ny < Nyd_; ++ny) {
                outfields[0].cmplx(mx, ny, mz, 0) = Ruk_[ny] - kappa2 * uk_[ny] - Dx * Pk_[ny];
                outfields[0].cmplx(mx, ny, mz, 1) = Rvk_[ny] - kappa2 * vk_[ny] - Pyk_[ny];
                outfields[0].cmplx(mx, ny, mz, 2) = Rwk_[ny] - kappa2 * wk_[ny] - Dz * Pk_[ny];
            }

            if (kx == 0 && kz == 0) {
                if (Ubaseyy_.length() > 0)
                    for (int ny = 0; ny < My_; ++ny)
                        outfields[0].cmplx(mx, ny, mz, 0) += Complex(flags_.nu * Ubaseyy_[ny], 0);
                if (Wbaseyy_.length() > 0)
                    for (int ny = 0; ny < My_; ++ny)
                        outfields[0].cmplx(mx, ny, mz, 2) += Complex(flags_.nu * Wbaseyy_[ny], 0);
                if (flags_.constraint == PressureGradient) {
                    outfields[0].cmplx(mx, 0, mz, 0) -= Complex(dPdxRef_, 0);
                    outfields[0].cmplx(mx, 0, mz, 2) -= Complex(dPdzRef_, 0);
                } else {
                    Real Ly = b_ - a_;
                    diff(uk_, Ruk_);
                    diff(wk_, Rwk_);
                    Real dPdxAct = Re(Ruk_.eval_b() - Ruk_.eval_a()) / Ly;
                    Real dPdzAct = Re(Rwk_.eval_b() - Rwk_.eval_a()) / Ly;
                    ChebyCoeff Ubasey = diff(Ubase_);
                    ChebyCoeff Wbasey = diff(Wbase_);
                    if (Ubase_.length() != 0)
                        dPdxAct += flags_.nu * (Ubasey.eval_b() - Ubasey.eval_a()) / Ly;
                    if (Wbase_.length() != 0)
                        dPdzAct += flags_.nu * (Wbasey.eval_b() - Wbasey.eval_a()) / Ly;
                    outfields[0].cmplx(mx, 0, mz, 0) -= Complex(dPdxAct, 0);
                    outfields[0].cmplx(mx, 0, mz, 2) -= Complex(dPdzAct, 0);
                }
            }

            for (int i = 0; i < nscalars(); ++i) {
                ComplexChebyCoeff& ck = ck_[i];
                for (int ny = 0; ny < Nyd_; ++ny)
                    ck.set(ny, scalars_[i].kappa * infields[1 + i].cmplx(mx, ny, mz, 0));
                diff2(ck, Rck_, Pyk_);
                for (int ny = 0; ny < Nyd_; ++ny)
                    outfields[1 + i].cmplx(mx, ny, mz, 0) = Rck_[ny] - kappa2 * ck[ny];
            }
        }
    }
}

void SettlingNSE::solve(vector<FlowField>& outfields, const vector<FlowField>& rhs, const int s) {
    // outfields {u, c_i, p}, rhs {Ru, Rc_i}
    assert(outfields.size() == rhs.size() + 1);
    FlowField& press = outfields.back();
    const int kxmax = outfields[0].kxmax();
    const int kzmax = outfields[0].kzmax();

    for (lint ix = 0; ix < Mxloc_; ++ix) {
        const lint mx = ix + mxlocmin_;
        const int kx = outfields[0].kx(mx);
        for (lint iz = 0; iz < Mzloc_; ++iz) {
            const lint mz = iz + mzlocmin_;
            const int kz = outfields[0].kz(mz);
            if ((kx == kxmax || kz == kzmax) || (flags_.dealias_xz() && isAliasedMode(kx, kz)))
                break;
            const bool mean = (kx == 0 && kz == 0);

            for (int ny = 0; ny < Nyd_; ++ny) {
                Ruk_.set(ny, rhs[0].cmplx(mx, ny, mz, 0));
                Rvk_.set(ny, rhs[0].cmplx(mx, ny, mz, 1));
                Rwk_.set(ny, rhs[0].cmplx(mx, ny, mz, 2));
            }

            // Scalars first, so that the buoyancy uses c^{n+1}
            for (int i = 0; i < nscalars(); ++i) {
                for (int ny = 0; ny < Nyd_; ++ny)
                    Rck_.set(ny, -rhs[1 + i].cmplx(mx, ny, mz, 0));
                const Real ga = mean ? 1.0 - scalars_[i].alpha() : 0.0;
                const Real gb = mean ? 1.0 + scalars_[i].alpha() : 0.0;
                scalarsolver_[i][s][ix][iz].solve(ck_[i].re, Rck_.re, ga, gb, 1.0);
                scalarsolver_[i][s][ix][iz].solve(ck_[i].im, Rck_.im, 0.0, 0.0, 0.0);
                for (int ny = 0; ny < Nyd_; ++ny)
                    Rvk_.sub(ny, scalars_[i].Ri * ck_[i][ny]);
            }

            if (prescribedVelocity_) {
                const bool realmode = mean || (outfields[0].Nx() % 2 == 0 && kx == kxmax && kz == 0) ||
                                      (outfields[0].Nz() % 2 == 0 && kz == kzmax && kx == 0) ||
                                      (outfields[0].Nx() % 2 == 0 && outfields[0].Nz() % 2 == 0 &&
                                       kx == kxmax && kz == kzmax);
                for (int i = 0; i < nscalars(); ++i)
                    for (int ny = 0; ny < Nyd_; ++ny)
                        outfields[1 + i].cmplx(mx, ny, mz, 0) = realmode ? Complex(Re(ck_[i][ny]), 0.0) : ck_[i][ny];
                continue;
            }

            if (!mean)
                tausolver_[s][ix][iz].solve(uk_, vk_, wk_, Pk_, Ruk_, Rvk_, Rwk_);
            else {
                if (Ubaseyy_.length() > 0)
                    for (int ny = 0; ny < My_; ++ny)
                        Ruk_.re[ny] += flags_.nu * Ubaseyy_[ny];
                if (Wbaseyy_.length() > 0)
                    for (int ny = 0; ny < My_; ++ny)
                        Rwk_.re[ny] += flags_.nu * Wbaseyy_[ny];
                if (flags_.constraint == PressureGradient) {
                    Ruk_.re[0] -= dPdxRef_;
                    Rwk_.re[0] -= dPdzRef_;
                    tausolver_[s][ix][iz].solve(uk_, vk_, wk_, Pk_, Ruk_, Rvk_, Rwk_);
                } else {  // bulk velocity: dPdx is a free variable fixed by UbulkBase + mean(u) = UbulkRef
                    tausolver_[s][ix][iz].solve(uk_, vk_, wk_, Pk_, dPdxAct_, dPdzAct_, Ruk_, Rvk_, Rwk_,
                                                UbulkRef_ - UbulkBase_, WbulkRef_ - WbulkBase_);
                }
            }

            // The (0,0) mode, and the kxmax/kzmax modes for even Nx/Nz, must be real (FFTW symmetry)
            const bool realmode = mean || (outfields[0].Nx() % 2 == 0 && kx == kxmax && kz == 0) ||
                                  (outfields[0].Nz() % 2 == 0 && kz == kzmax && kx == 0) ||
                                  (outfields[0].Nx() % 2 == 0 && outfields[0].Nz() % 2 == 0 && kx == kxmax &&
                                   kz == kzmax);
            auto put = [realmode](Complex z) { return realmode ? Complex(Re(z), 0.0) : z; };
            for (int ny = 0; ny < Nyd_; ++ny) {
                outfields[0].cmplx(mx, ny, mz, 0) = put(uk_[ny]);
                outfields[0].cmplx(mx, ny, mz, 1) = put(vk_[ny]);
                outfields[0].cmplx(mx, ny, mz, 2) = put(wk_[ny]);
                for (int i = 0; i < nscalars(); ++i)
                    outfields[1 + i].cmplx(mx, ny, mz, 0) = put(ck_[i][ny]);
                press.cmplx(mx, ny, mz, 0) = put(Pk_[ny]);
            }
        }
    }
}

void SettlingNSE::reset_lambda(const vector<Real> lambda_t) {
    if (tausolver_ != 0 && lambda_t.size() != lambda_t_.size())
        cferror("SettlingNSE::reset_lambda : the number of substeps cannot change");
    lambda_t_ = lambda_t;
    const uint nsub = lambda_t.size();

    if (tausolver_ == 0) {
        tausolver_ = new TauSolver**[nsub];
        for (uint j = 0; j < nsub; ++j) {
            tausolver_[j] = new TauSolver*[Mxloc_];
            for (int mx = 0; mx < Mxloc_; ++mx)
                tausolver_[j][mx] = new TauSolver[Mzloc_];
        }
    }
    for (auto& solver : scalarsolver_)
        if (!solver) {
            solver = new ScalarSolver**[nsub];
            for (uint j = 0; j < nsub; ++j) {
                solver[j] = new ScalarSolver*[Mxloc_];
                for (int mx = 0; mx < Mxloc_; ++mx)
                    solver[j][mx] = new ScalarSolver[Mzloc_];
            }
        }

    for (uint j = 0; j < nsub; ++j)
        for (int mx = 0; mx < Mxloc_; ++mx) {
            const int kx = kxloc_[mx];
            for (int mz = 0; mz < Mzloc_; ++mz) {
                const int kz = kzloc_[mz];
                if ((kx == kxmax_ && kz == kzmax_) || (flags_.dealias_xz() && isAliasedMode(kx, kz)))
                    continue;
                const Real k2 = 4 * square(pi) * (square(kx / Lx_) + square(kz / Lz_));
                tausolver_[j][mx][mz] =
                    TauSolver(kx, kz, Lx_, Lz_, a_, b_, lambda_t[j] + flags_.nu * k2, flags_.nu, Nyd_,
                              flags_.taucorrection);
                for (int i = 0; i < nscalars(); ++i) {
                    const SettlingScalar& sc = scalars_[i];
                    scalarsolver_[i][j][mx][mz] =
                        ScalarSolver(Nyd_, sc.bc(Mxloc_, Mzloc_), a_, b_, lambda_t[j] + sc.kappa * k2,
                                             sc.kappa, kx == 0 && kz == 0);
                }
            }
        }
}

vector<FlowField> SettlingNSE::createRHS(const vector<FlowField>& fields) const {
    return vector<FlowField>(fields.begin(), fields.begin() + 1 + nscalars());
}

vector<cfarray<FieldSymmetry>> SettlingNSE::createSymmVec() const {
    // The same symmetries act on the velocity and on each scalar (a 1-component field is unsigned
    // except for the overall sign s); none on the pressure.
    cfarray<FieldSymmetry> usym = SymmetryList(1);
    if (flags_.symmetries.length() > 0)
        usym = flags_.symmetries;
    vector<cfarray<FieldSymmetry>> symm(1 + nscalars(), usym);
    symm.push_back(SymmetryList(1));
    return symm;
}

}  // namespace chflow
