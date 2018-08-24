/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#ifndef CHANNELFLOW_REALPROFILE_H
#define CHANNELFLOW_REALPROFILE_H
#include <vector>
#include "cfbasics/mathdefs.h"
#include "channelflow/basisfunc.h"

namespace channelflow {

enum Sign { Minus = -1, Plus = 1 };
std::ostream& operator<<(std::ostream&, Sign s);
std::istream& operator>>(std::istream&, Sign& s);

// RealProfile is used primarily to represent real-valued "Fourier profiles",
// Due to the complex representation of the Fourier modes of FlowFields,
// a RealProfile is interally expressed as a linear combination of a
// ComplexProfile and its complex conjugate --this is the form needed
// to calculate things like inner products against velocity fields.

// A RealProfile f can be in one of two canonical states
// f = (psi + psi*)
// f = (psi - psi*)/i
// where psi is a complex-valued vector BasisFunction, and with psi having
// kx,kz==0,0;   (in which case f == psi+psi* always).
// kz=0, kx>0; or
// kz>0, kx positive, negative, or zero
// and, usually, psi normalized so that L2Norm(f) == 1.

class RealProfile {
   public:
    RealProfile();
    RealProfile(const std::string& filebase);

    RealProfile(int Nd, int Ny, int kx, int kz, cfbasics::Real Lx, cfbasics::Real Lz, cfbasics::Real a,
                cfbasics::Real b, Sign sign, cfbasics::fieldstate state = cfbasics::Spectral);

    RealProfile(int Ny, const RealProfile& f);

    RealProfile(const ComplexChebyCoeff& u, const ComplexChebyCoeff& v, const ComplexChebyCoeff& w, int kx, int kz,
                cfbasics::Real Lx, cfbasics::Real Lz, Sign s);

    // RealProfile == (psi + psi^*)   == 2 Re psi, sign == Plus
    //             == (psi - psi^*)/i == 2 Im psi, sign == Minus
    RealProfile(const BasisFunc& phi, Sign s);

    void save(const std::string& filebase, cfbasics::fieldstate s = cfbasics::Physical) const;
    void binaryDump(std::ostream& os) const;
    void binaryLoad(std::istream& is);

    void canonicalize();  // swap signs and conjugate to put in canonical form
    void randomize(cfbasics::Real magn, cfbasics::Real decay, BC aBC, BC bBC);
    void interpolate(const RealProfile& f);
    void reflect(const RealProfile& f);
    void fill(const RealProfile& f);

    inline int Nd() const;
    inline int Ny() const;
    inline int kx() const;
    inline int kz() const;
    inline cfbasics::Real Lx() const;
    inline cfbasics::Real Lz() const;
    inline cfbasics::Real a() const;
    inline cfbasics::Real b() const;
    inline Sign sign() const;
    inline cfbasics::fieldstate state() const;

    void reconfig(const RealProfile& f);
    void resize(int Ny, int Nd);
    void setBounds(cfbasics::Real Lx, cfbasics::Real Lz, cfbasics::Real a, cfbasics::Real b);
    void setkxkzSign(int kx, int kz, Sign s);

    void setState(cfbasics::fieldstate s);
    void setSign(Sign s);
    void setToZero();

    void chebyfft();
    void ichebyfft();
    void makeSpectral();
    void makePhysical();
    void makeState(cfbasics::fieldstate s);

    void chebyfft(const ChebyTransform& t);
    void ichebyfft(const ChebyTransform& t);
    void makeSpectral(const ChebyTransform& t);
    void makePhysical(const ChebyTransform& t);
    void makeState(cfbasics::fieldstate s, const ChebyTransform& t);

    // return the ith component as a new RealProfile
    RealProfile operator[](int i) const;

    bool geomCongruent(const RealProfile& f) const;
    bool congruent(const RealProfile& f) const;      // geomCongruent & kx,kz equality
    bool interoperable(const RealProfile& f) const;  // congruent && same state

    RealProfile& operator*=(const FieldSymmetry& s);
    RealProfile& operator*=(const RealProfile& g);
    RealProfile& operator+=(const RealProfile& g);
    RealProfile& operator-=(const RealProfile& g);
    RealProfile& operator*=(cfbasics::Real c);
    RealProfile& operator*=(cfbasics::Complex c);

    BasisFunc psi;  // The psi from which (psi+psi*) or (psi-psi*)/i is formed

   private:
    Sign sign_;  // Minus, Plus
};

std::vector<RealProfile> realBasisKxKz(int Ny, int kx, int kz, cfbasics::Real Lx, cfbasics::Real Lz, cfbasics::Real a,
                                       cfbasics::Real b, const BasisFlags& flags);

std::vector<RealProfile> realBasis(int Ny, int kxmax, int kzmax, cfbasics::Real Lx, cfbasics::Real Lz, cfbasics::Real a,
                                   cfbasics::Real b, const BasisFlags& flags);

void orthonormalize(std::vector<RealProfile>& basis);

void checkBasis(const std::vector<RealProfile>& e, const BasisFlags& flags, bool orthogcheck = false);

bool operator==(const RealProfile& f, const RealProfile& g);
bool operator!=(const RealProfile& f, const RealProfile& g);

cfbasics::Real L2Norm(const RealProfile& f, bool normalize = true);
cfbasics::Real L2Norm2(const RealProfile& f, bool normalize = true);
cfbasics::Real L2Dist(const RealProfile& f, const RealProfile& g, bool normalize = true);
cfbasics::Real L2Dist2(const RealProfile& f, const RealProfile& g, bool normalize = true);
cfbasics::Real L2InnerProduct(const RealProfile& f, const RealProfile& g, bool normalize = true);

cfbasics::Real chebyNorm(const RealProfile& f, bool normalize = true);
cfbasics::Real chebyNorm2(const RealProfile& f, bool normalize = true);
cfbasics::Real chebyDist(const RealProfile& f, const RealProfile& g, bool normalize = true);
cfbasics::Real chebyDist2(const RealProfile& f, const RealProfile& g, bool normalize = true);
cfbasics::Real chebyInnerProduct(const RealProfile& f, const RealProfile& g, bool normalize = true);

cfbasics::Real norm(const RealProfile& f, NormType n, bool normalize = true);
cfbasics::Real norm2(const RealProfile& f, NormType n, bool normalize = true);
cfbasics::Real dist(const RealProfile& f, const RealProfile& g, NormType n, bool nrmlz = true);
cfbasics::Real dist2(const RealProfile& f, const RealProfile& g, NormType n, bool nrmlz = true);
cfbasics::Real innerProduct(const RealProfile& f, const RealProfile& g, NormType n, bool normalize = true);

cfbasics::Real bcNorm(const RealProfile& f, bool normalize = true);
cfbasics::Real bcDist(const RealProfile& f, const RealProfile& g, bool normalize = true);
cfbasics::Real bcNorm2(const RealProfile& f, bool normalize = true);
cfbasics::Real bcDist2(const RealProfile& f, const RealProfile& g, bool normalize = true);

cfbasics::Real divNorm(const RealProfile& f, bool normalize = true);
cfbasics::Real divDist(const RealProfile& f, const RealProfile& g, bool normalize = true);
cfbasics::Real divNorm2(const RealProfile& f, bool normalize = true);
cfbasics::Real divDist2(const RealProfile& f, const RealProfile& g, bool normalize = true);

RealProfile xdiff(const RealProfile& f);
RealProfile ydiff(const RealProfile& f);
RealProfile zdiff(const RealProfile& f);
RealProfile grad(const RealProfile& f);
RealProfile curl(const RealProfile& f);
RealProfile lapl(const RealProfile& f);
RealProfile div(const RealProfile& f);

// Binary RealProfile operations produce TWO RealProfiles!
// So we can't define binary ops like this. Use void versions below.
//
// RealProfile dotgrad(const RealProfile& f, const RealProfile& g);
// RealProfile cross(const RealProfile& f, const RealProfile& g);
// ChebyCoeff  dot(const RealProfile& f, const RealProfile& g);

void xdiff(const RealProfile& f, RealProfile& fx);
void ydiff(const RealProfile& f, RealProfile& fy);
void zdiff(const RealProfile& f, RealProfile& fz);
void grad(const RealProfile& f, RealProfile& gradf);
void curl(const RealProfile& f, RealProfile& curlf);
void lapl(const RealProfile& f, RealProfile& laplf);
void div(const RealProfile& f, RealProfile& divf);

void cross(const RealProfile& f, const RealProfile& g, RealProfile& rtn1, RealProfile& rtn2);

void dot(const RealProfile& f, const RealProfile& g, RealProfile& rtn1, RealProfile& rtn2);

void dotgrad(const RealProfile& f, const RealProfile& g, RealProfile& rtn1, RealProfile& rtn2);

// For a multiplicative operation on binary operation, f x g,
// where  f = a+a* or (a-a*)/i and
//        g = b+b* or (b-b*)/i,
// let rtn1.psi = a x b
// let rtn2.psi = a x b*
// Then multiplicativeOpSetSigns modifies the Sign members of rtn1 and rtn2
// and possibly mutlipies them by -1 in order to set the correct form of
// rtn1 = a x b +- (a x b)*
// rtn2 = (+- a x b* +- (a x b*)*/i
// so that
// f x g = rtn1 + rtn2

void multiplySetSigns(Sign fsign, Sign gsign, RealProfile& rtn1, RealProfile& rtn2);

// An expensive convenience.

inline int RealProfile::Nd() const { return psi.Nd(); }
inline int RealProfile::Ny() const { return psi.Ny(); }
inline int RealProfile::kx() const { return psi.kx(); }
inline int RealProfile::kz() const { return psi.kz(); }
inline cfbasics::Real RealProfile::Lx() const { return psi.Lx(); }
inline cfbasics::Real RealProfile::Lz() const { return psi.Lz(); }
inline cfbasics::Real RealProfile::a() const { return psi.a(); }
inline cfbasics::Real RealProfile::b() const { return psi.b(); }
inline Sign RealProfile::sign() const { return sign_; }
inline cfbasics::fieldstate RealProfile::state() const { return psi.state(); }

inline cfbasics::Real L2IP(const RealProfile& f, const RealProfile& g, bool normalize = true) {
    return L2InnerProduct(f, g, normalize);
}

}  // namespace channelflow
#endif
