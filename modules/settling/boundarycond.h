/**
 * Boundary conditions of a scalar field in y (Dirichlet, Dirichlet/Robin, Neumann/Robin, no flux)
 * and the functions that impose them on Chebyshev expansions.
 *
 * Robin conditions are f' + alpha f = g. No flux, kappa c' + vs c = 0, is a Robin condition at both
 * walls with alpha = vs/kappa (for the stored field c + y, g = 1 -/+ alpha for the mean mode).
 * The settling module uses NoFlux; the other types are kept for later use and are covered by the
 * module's unit tests only.
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */

#ifndef SETTLING_BOUNDARYCOND_H
#define SETTLING_BOUNDARYCOND_H

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"

namespace chflow {

enum class ScalarBCType { Diri, DiriRobin, NeumRobin, NoFlux };

class BoundaryCond {
   public:
    BoundaryCond();
    BoundaryCond(ScalarBCType type, int Mx, int Mz, Real ga, Real gb);
    BoundaryCond(ScalarBCType type, int Mx, int Mz, Real ga, Real gb, Real alpha);
    void setgaxz(Vector& g_re, Vector& g_im, int Mx);
    Complex gaxz(int mx, int mz) const;

    ScalarBCType type_;
    int Mx_;
    int Mz_;
    Real ga_;
    Real gb_;
    Real alpha_;

    // these are required in the case of nonconstant ga
    Vector gaxz_re_;
    Vector gaxz_im_;
};

// Impose the homogeneous (first two) or inhomogeneous (ga, gb) form of bc on f, by changing its
// first two Chebyshev coefficients
void fixBC(ChebyCoeff& f, BoundaryCond bc);
void fixBC(ComplexChebyCoeff& f, BoundaryCond bc);
void fixBC(ChebyCoeff& f, Real ga, Real gb, BoundaryCond bc);
void fixBC(ComplexChebyCoeff& f, Complex ga, Real gb, BoundaryCond bc);
void fixDiri(ChebyCoeff& f, Real ga, Real gb);
void fixDiriRobin(ChebyCoeff& f, Real ga = 0, Real gb = 0, Real alpha = 0);
void fixNeumRobin(ChebyCoeff& f, Real ga = 0, Real gb = 0, Real alpha = 0);
void fixNoFlux(ChebyCoeff& f, Real ga = 0, Real gb = 0, Real alpha = 0);
// fixNoFlux and mean(f) = 1, by changing coefficients 0, 1 and 2
void fixNoFluxNormaliseMass(ChebyCoeff& f, Real ga = 0, Real gb = 0, Real alpha = 0);

}  // namespace chflow

#endif
