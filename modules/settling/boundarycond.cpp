/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#include "modules/settling/boundarycond.h"
#include "channelflow/utilfuncs.h"

using namespace std;

namespace chflow {

BoundaryCond::BoundaryCond()
    : type_(ScalarBCType::Diri), Mx_(0), Mz_(0), ga_(0), gb_(0), alpha_(0) {
}

BoundaryCond::BoundaryCond(ScalarBCType type, int Mx, int Mz, Real ga, Real gb)
    : type_(type), Mx_(Mx), Mz_(Mz), ga_(ga), gb_(gb), alpha_(0),
    gaxz_re_(Mx_ * Mz_), gaxz_im_(Mx_ * Mz_) { 
    gaxz_re_[0] = ga;
}

BoundaryCond::BoundaryCond(ScalarBCType type, int Mx, int Mz, Real ga, Real gb, Real alpha)
    : type_(type), Mx_(Mx), Mz_(Mz), ga_(ga), gb_(gb), alpha_(alpha),
    gaxz_re_(Mx_ * Mz_), gaxz_im_(Mx_ * Mz_) { 
    gaxz_re_[0] = ga;
}

// Set Fourier components of ga - the bottom wall boundary condition
void BoundaryCond::setgaxz(Vector& g_re, Vector& g_im, int Mx)
{
    gaxz_re_ = g_re;
    gaxz_im_ = g_im;
    Mx_ = Mx;
}

Complex BoundaryCond::gaxz(int mx, int mz) const
{
    return Complex(gaxz_re_[mx + mz * Mx_], gaxz_im_[mx + mz * Mx_]);
}

void fixBC(ChebyCoeff& f, BoundaryCond bc) {
    if (bc.type_ == ScalarBCType::Diri)
        fixDiri(f);
    else if (bc.type_ == ScalarBCType::DiriRobin)
        fixDiriRobin(f, 0.0, 0.0, bc.alpha_); // JL important to retain alpha here!
    else if (bc.type_ == ScalarBCType::NeumRobin)
        fixNeumRobin(f, 0.0, 0.0, bc.alpha_); 
    else if (bc.type_ == ScalarBCType::NoFlux)
        fixNoFlux(f, 0.0, 0.0, bc.alpha_); 
}

void fixBC(ChebyCoeff& f, Real ga, Real gb, BoundaryCond bc) {
    if (bc.type_ == ScalarBCType::Diri)
        fixDiri(f, ga, gb);
    else if (bc.type_ == ScalarBCType::DiriRobin)
        fixDiriRobin(f, ga, gb, bc.alpha_);
    else if (bc.type_ == ScalarBCType::NoFlux)
        fixNoFlux(f, ga, gb, bc.alpha_);
}

// JL cater for nonzero Dirichlet boundary values ga, gb
void fixDiri(ChebyCoeff& f, Real ga, Real gb) {
    Real fa = f.eval_a();
    Real fb = f.eval_b();
    Real mean = 0.5 * (fb + fa);
    Real slop = 0.5 * (fb - fa);
    f[0] -= mean - 0.5 * (gb + ga);
    f[1] -= slop - 0.5 * (gb - ga);
}

// JL mixed condition with Dirichlet at ya and Robin at yb defined by
// df/dy + alpha * f = gb
void fixDiriRobin(ChebyCoeff& f, Real ga, Real gb, Real alpha) {
    Real fa = f.eval_a();
    Real fb = f.eval_b();
    Real sb = f.slope_b();

    f[0] -= (sb + (1.0 + alpha) * fa + alpha * fb - (gb + (1.0 + alpha) * ga)) / (1.0 + 2.0 * alpha);
    f[1] -= (sb + alpha * fb - alpha * fa - gb + alpha * ga) / (1.0 + 2.0 * alpha);
}

// JL mixed condition with Neumann  at ya and Robin at yb
void fixNeumRobin(ChebyCoeff& f, Real ga, Real gb, Real alpha) {
    Real fb = f.eval_b();
    Real sa = f.slope_a();
    Real sb = f.slope_b();

    f[0] -= (sb - (1.0 + alpha) * sa + alpha * fb - gb + (1.0 + alpha) * ga) / alpha;
    f[1] -= -ga + sa;
}

// JL no flux condition with Robin bcs top and bottom
void fixNoFlux(ChebyCoeff& f, Real ga, Real gb, Real alpha) {
    Real fa = f.eval_a();
    Real fb = f.eval_b();
    Real sa = f.slope_a();
    Real sb = f.slope_b();

    Real tmp = (alpha + 1) * ga + (alpha - 1) * gb;
    for (int n = 2u; n < f.N(); n++)
        tmp -= (alpha - 1) * (alpha + n * n) * f[n];
    for (int n = 2u; n < f.N(); n += 2)
        tmp -= (alpha + 1) * (alpha - n * n) * f[n];
    for (int n = 3u; n < f.N(); n += 2)
        tmp += (alpha + 1) * (alpha - n * n) * f[n];
    tmp /= 2.0 * alpha * alpha;
    //f[0] -= (alpha - 1) * (alpha * fb + sb - gb) / 2.0 / alpha + (alpha + 1) * (alpha * fa + sa - ga) / 2.0 / alpha;
    f[0] = tmp;
    f[1] -= (sb + alpha * fb - sa - alpha * fa + ga - gb) / (2.0 * alpha);
}

void fixNoFluxNormaliseMass(ChebyCoeff& f, Real ga, Real gb, Real alpha) {
    // probably this can be simplified
    Real ap1_o_2asqr = (alpha + 1) / (2 * alpha * alpha);
    Real am1_o_2asqr = (alpha - 1) / (2 * alpha * alpha);
    Real tmp = 1.0 - ga * ap1_o_2asqr - gb * am1_o_2asqr;
    for (int n = 4u; n < f.N(); n += 2)
        tmp -= f[n] / (1 - n * n) - ap1_o_2asqr * (alpha - n * n) * f[n];
    for (int n = 3u; n < f.N(); n += 2)
        tmp -= ap1_o_2asqr * (alpha - n * n) * f[n];
    for (int n = 3u; n < f.N(); n ++)
        tmp += am1_o_2asqr * (alpha + n * n) * f[n];
    f[2] = tmp / (-ap1_o_2asqr * (alpha - 4) - am1_o_2asqr * (alpha + 4) - 1.0 / 3.0);

    tmp = 1.0;
    for (int n = 2u; n < f.N(); n += 2)
        tmp -= f[n] / (1 - n * n);
    f[0] = tmp;
    //f[0] -= (alpha - 1) * (alpha * fb + sb - gb) / 2.0 / alpha + (alpha + 1) * (alpha * fa + sa - ga) / 2.0 / alpha;
    // note computing these *after* fixing f[2]
    Real fa = f.eval_a();
    Real fb = f.eval_b();
    Real sa = f.slope_a();
    Real sb = f.slope_b();
    f[1] -= (sb + alpha * fb - sa - alpha * fa + ga - gb) / (2.0 * alpha);
}

void fixBC(ComplexChebyCoeff& f, BoundaryCond bc) {
    fixBC(f.re, bc);
    fixBC(f.im, bc);
}

void fixBC(ComplexChebyCoeff& f, Complex ga, Real gb, BoundaryCond bc) {
    fixBC(f.re, ga.real(), gb, bc);
    fixBC(f.im, ga.imag(), gb, bc);
}

}  // namespace chflow
