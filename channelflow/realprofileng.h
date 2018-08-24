/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#ifndef CHANNELFLOW_REALPROFILE_NG
#define CHANNELFLOW_REALPROFILE_NG

#include <vector>
#include "channelflow/chebyshev.h"
#include "channelflow/symmetry.h"

namespace channelflow {

/**
 * Represents a vector field on the xz-periodic domain [0, Lx] x [a,b] x [0,Lz]
 * of the form :       Phi_u(y) g_jx{alpha x} g_-jz{gamma z} xhat +
 *                             Phi_v(y) g_-jx{alpha x} g_-jz{gamma z} yhat +
 *                              Phi_w(y) g_-jx{alpha x} g_jz{gamma z} zhat
 *
 * Where               [ cos(j * x) for  j >= 0
 *            g_j(x) = {
 *                         [ sin(-j * x) for J < 0
 */
class RealProfileNG {
   public:
    RealProfileNG(const int jx, const int jz, const int Nd, const int Ny, const cfbasics::Real Lx,
                  const cfbasics::Real Lz, const cfbasics::Real a, const cfbasics::Real b,
                  const cfbasics::fieldstate state = cfbasics::Spectral);

    RealProfileNG(const std::vector<ChebyCoeff> u, const int jx, const int jz, const cfbasics::Real Lx,
                  const cfbasics::Real Lz);

    RealProfileNG();

    RealProfileNG(const RealProfileNG&);

    inline int jx() const;
    inline int jz() const;
    inline int Nd() const;
    inline int Ny() const;
    inline cfbasics::Real Lx() const;
    inline cfbasics::Real Lz() const;
    inline cfbasics::Real a() const;
    inline cfbasics::Real b() const;
    inline cfbasics::fieldstate state() const;
    inline void setJx(int jx);
    inline void setJz(int jz);
    RealProfileNG& operator=(const RealProfileNG&);
    RealProfileNG& operator*=(const cfbasics::Real c);
    RealProfileNG& operator+=(const RealProfileNG& e);
    RealProfileNG& operator-=(const RealProfileNG& e);
    RealProfileNG& operator*=(const FieldSymmetry& s);
    inline const ChebyCoeff& operator[](int i) const;
    inline ChebyCoeff& operator[](int i);

    // Can be added/subtracted to e (true if congruent and on same geometry)
    bool compatible(const RealProfileNG& e) const;
    // True if on some geometry and in same state
    bool congruent(const RealProfileNG& e) const;

    void makeSpectral();                           // if Physical, transform to Spectral
    void makePhysical();                           // if Spectral, transform to Physical
    void makeState(const cfbasics::fieldstate s);  // if state != s, transform to state s

    void makeSpectral(const ChebyTransform& t);
    void makePhysical(const ChebyTransform& t);
    void makeState(const cfbasics::fieldstate s, const ChebyTransform& t);

    // When converting to a FlowField
    // Gives the appropriate normalization factor for the (+kx,kz) fourier mode
    cfbasics::Complex normalization_p(const int d) const;

    // Gives the appropriate normalization factor for the (-kx,kz) fourier mode
    cfbasics::Complex normalization_m(const int d) const;

   private:
    cfbasics::fieldstate state_;
    int jx_;
    int jz_;
    int Nd_;
    int Ny_;
    cfbasics::Real Lx_;
    cfbasics::Real Lz_;
    cfbasics::Real a_;
    cfbasics::Real b_;

   public:
    std::vector<ChebyCoeff> u_;
};

inline const ChebyCoeff& RealProfileNG::operator[](int i) const {
    assert(i >= 0 && i < Nd());
    return u_[i];
}
inline ChebyCoeff& RealProfileNG::operator[](int i) {
    assert(i >= 0 && i < Nd());
    return u_[i];
}

inline int RealProfileNG::jx() const { return jx_; }
inline int RealProfileNG::jz() const { return jz_; }
inline int RealProfileNG::Nd() const { return Nd_; }
inline int RealProfileNG::Ny() const { return Ny_; }
inline cfbasics::Real RealProfileNG::Lx() const { return Lx_; }
inline cfbasics::Real RealProfileNG::Lz() const { return Lz_; }
inline cfbasics::Real RealProfileNG::a() const { return a_; }
inline cfbasics::Real RealProfileNG::b() const { return b_; }
inline cfbasics::fieldstate RealProfileNG::state() const { return state_; }
inline void RealProfileNG::setJx(int jx) { jx_ = jx; }
inline void RealProfileNG::setJz(int jz) { jz_ = jz; }

cfbasics::Real L2InnerProduct(const RealProfileNG& e1, const RealProfileNG& e2, const bool normalize = true);
cfbasics::Real L2Norm2(const RealProfileNG& e, bool normalize = true);
inline cfbasics::Real L2Norm(const RealProfileNG& e, bool normalize = true) { return sqrt(L2Norm2(e, normalize)); }

std::vector<RealProfileNG> realBasisNG(const int Ny, const int kxmax, const int kzmax, const cfbasics::Real Lx,
                                       const cfbasics::Real Lz, const cfbasics::Real a, const cfbasics::Real b);

void orthonormalize(std::vector<RealProfileNG>& basis);

// Remove all elements from basis which are not symmetric under all members of s, to a given tolerence
void selectSymmetries(std::vector<RealProfileNG>& basis, const std::vector<FieldSymmetry>& s,
                      const cfbasics::Real tolerance = 1e-13);
}  // namespace channelflow

#endif
