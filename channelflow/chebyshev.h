/**
 * Real- and Complex-valued Chebyshev expansion classes.
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#ifndef CHANNELFLOW_CHEBYSHEV_H
#define CHANNELFLOW_CHEBYSHEV_H

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/cfmpi.h"

#include <fftw3.h>
#ifdef HAVE_MPI
#include <fftw3-mpi.h>
#endif

#include <memory>
#include <type_traits>

namespace channelflow {

enum BC { Free, Diri };
enum NormType { Uniform, Chebyshev };  // Uniform weighting in y or 1/sqrt(1-y^2)

void fftw_loadwisdom(const char* filename = 0);  // defaults to ~/.fftw_wisdom
void fftw_savewisdom(const char* filename = 0);

cfbasics::Vector chebypoints(int N, cfbasics::Real a, cfbasics::Real b);  // N   pts

cfbasics::Real chebyIP(int m, int n);  // Integral_(-1)^1 Tm(y) Tn(y) dy
inline int cheby_c(int n) { return (n > 0) ? 1 : ((n == 0) ? 2 : 0); }

cfbasics::Real legendre(int n, cfbasics::Real x);   // Value of Legendre polynomial P_n(x)
cfbasics::Real chebyshev(int n, cfbasics::Real x);  // Value of Chebyshev polynomial T_n(x)

void gaussLegendreQuadrature(int N, cfbasics::Real a, cfbasics::Real b, cfbasics::Vector& x, cfbasics::Vector& w);

class ChebyTransform;

class ChebyCoeff : public cfbasics::Vector {
   public:
    ChebyCoeff();  /// Null constructor
    ChebyCoeff(int N, cfbasics::Real a, cfbasics::Real b, cfbasics::fieldstate s = cfbasics::Spectral);
    ChebyCoeff(const Vector& v, cfbasics::Real a, cfbasics::Real b, cfbasics::fieldstate s = cfbasics::Spectral);
    ChebyCoeff(int N, const ChebyCoeff& g);   // copy first N elems.
    ChebyCoeff(const std::string& filebase);  // read ascii from file
    ~ChebyCoeff();

    void save(const std::string& filebase, cfbasics::fieldstate s = cfbasics::Physical) const;
    void binaryDump(std::ostream& os) const;
    void binaryLoad(std::istream& is);

    void reconfig(const ChebyCoeff& f);
    void randomize(cfbasics::Real magn, cfbasics::Real smoothness, BC aBC, BC bBC);
    void setBounds(cfbasics::Real a, cfbasics::Real b);
    void setState(cfbasics::fieldstate s);
    void setToZero();
    void fill(const ChebyCoeff& g);
    void interpolate(const ChebyCoeff& g);                  // *this is on subdomain of g.
    void reflect(const ChebyCoeff& g, cfbasics::parity p);  // reflect g about u.a()

    cfbasics::Real eval_a() const;
    cfbasics::Real eval_b() const;
    cfbasics::Real eval(cfbasics::Real x) const;
    ChebyCoeff eval(const Vector& x) const;
    void eval(const Vector& x, ChebyCoeff& g) const;
    cfbasics::Real slope_a() const;
    cfbasics::Real slope_b() const;

    inline cfbasics::Real a() const;
    inline cfbasics::Real b() const;
    inline cfbasics::Real L() const;
    inline int N() const;
    inline int numModes() const;
    inline cfbasics::fieldstate state() const;
    cfbasics::Real mean() const;
    // Real energy() const;

    ChebyCoeff& operator*=(cfbasics::Real c);
    ChebyCoeff& operator+=(const ChebyCoeff& g);
    ChebyCoeff& operator-=(const ChebyCoeff& g);
    ChebyCoeff& operator*=(const ChebyCoeff& g);  // dottimes, only for Physical

    // These transforms make temp ChebyTransforms thus are more expensive.
    void chebyfft();                         // transform from Physcial to Spectral
    void ichebyfft();                        // transform from Spectral to Physcial
    void makeSpectral();                     // if Physical, transform to Spectral
    void makePhysical();                     // if Spectral, transform to Physical
    void makeState(cfbasics::fieldstate s);  // if state != s, transform to state s

    void chebyfft(const ChebyTransform& t);
    void ichebyfft(const ChebyTransform& t);
    void makeSpectral(const ChebyTransform& t);
    void makePhysical(const ChebyTransform& t);
    void makeState(cfbasics::fieldstate s, const ChebyTransform& t);

    bool congruent(const ChebyCoeff& g) const;

    friend void swap(ChebyCoeff& f, ChebyCoeff& g);

   private:
    cfbasics::Real a_;            // lower bound of domain
    cfbasics::Real b_;            // upper bound of domain
    cfbasics::fieldstate state_;  // indicates Physical or Spectral state of object
    friend class ChebyTransform;
};

class ComplexChebyCoeff {
   public:
    ComplexChebyCoeff();
    ComplexChebyCoeff(int N, cfbasics::Real a, cfbasics::Real b, cfbasics::fieldstate s);
    ComplexChebyCoeff(int N, const ComplexChebyCoeff& f);  // copy first N elems.
    ComplexChebyCoeff(const ChebyCoeff& re, const ChebyCoeff& im);
    ComplexChebyCoeff(const std::string& filename);  // read ascii from file

    void reconfig(const ComplexChebyCoeff& f);
    void resize(int N);
    void randomize(cfbasics::Real magn, cfbasics::Real smoothness, BC aBC, BC bBC);
    void setToZero();
    void setBounds(cfbasics::Real a, cfbasics::Real b);
    void setState(cfbasics::fieldstate s);
    void fill(const ComplexChebyCoeff& g);                         // set this[i]=v[i], rest=0
    void interpolate(const ComplexChebyCoeff& g);                  // set this[i]=v[i], rest=0
    void reflect(const ComplexChebyCoeff& g, cfbasics::parity p);  // reflect v about u.a()

    cfbasics::Complex eval_a() const;
    cfbasics::Complex eval_b() const;
    cfbasics::Complex eval(cfbasics::Real x) const;
    cfbasics::Complex slope_a() const;
    cfbasics::Complex slope_b() const;

    cfbasics::Complex mean() const;

    inline cfbasics::Real a() const;
    inline cfbasics::Real b() const;
    inline cfbasics::Real L() const;
    inline int N() const;
    inline int length() const;
    inline int numModes() const;
    inline cfbasics::fieldstate state() const;

    inline cfbasics::Complex operator[](int n) const;
    inline void set(int n, cfbasics::Complex c);
    inline void add(int n, cfbasics::Complex c);
    inline void sub(int n, cfbasics::Complex c);

    ComplexChebyCoeff& operator+=(const ComplexChebyCoeff& f);
    ComplexChebyCoeff& operator-=(const ComplexChebyCoeff& f);
    ComplexChebyCoeff& operator*=(cfbasics::Real c);
    ComplexChebyCoeff& operator*=(cfbasics::Complex c);
    ComplexChebyCoeff& operator*=(const ComplexChebyCoeff& c);  // dottimes

    void conjugate();  // destructive
    void save(const std::string& filebase, cfbasics::fieldstate s = cfbasics::Physical) const;
    void binaryDump(std::ostream& os) const;
    void binaryLoad(std::istream& is);

    bool congruent(const ComplexChebyCoeff& g) const;

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

    friend void swap(ComplexChebyCoeff& f, ComplexChebyCoeff& g);

    ChebyCoeff re;
    ChebyCoeff im;
};

class ChebyTransform {
   public:
    ChebyTransform(int N, uint fftw_flags = FFTW_ESTIMATE);
    ChebyTransform(const ChebyTransform& t);             // unimplemented
    ChebyTransform& operator=(const ChebyTransform& t);  // unimplemented

    inline int N() const;
    inline int length() const;

    friend class ChebyCoeff;
    friend class ComplexChebyCoeff;

   private:
    int N_;
    uint flags_;
    std::unique_ptr<std::remove_pointer<fftw_plan>::type, void (*)(fftw_plan)> cosfftw_plan_;
};

ChebyCoeff operator*(cfbasics::Real c, const ChebyCoeff& g);
ChebyCoeff operator+(const ChebyCoeff& f, const ChebyCoeff& g);
ChebyCoeff operator-(const ChebyCoeff& f, const ChebyCoeff& g);
bool operator==(const ChebyCoeff& f, const ChebyCoeff& g);
bool operator!=(const ChebyCoeff& f, const ChebyCoeff& g);

void diff(const ChebyCoeff& f, ChebyCoeff& df);
void diff2(const ChebyCoeff& f, ChebyCoeff& d2f);
void diff2(const ChebyCoeff& f, ChebyCoeff& d2f, ChebyCoeff& tmp);
void diff(const ChebyCoeff& f, ChebyCoeff& df, int n);
ChebyCoeff diff(const ChebyCoeff& f);
ChebyCoeff diff2(const ChebyCoeff& f);
ChebyCoeff diff(const ChebyCoeff& f, int n);

// Integrate sets the arbitrary const of integration so that mean(f)==0.
void integrate(const ChebyCoeff& df, ChebyCoeff& f);
ChebyCoeff integrate(const ChebyCoeff& df);

void legendre(int n, ChebyCoeff& u, ChebyTransform& trans, bool normalize = false);

//   normalize  ?   false                  :  true
// L2Norm2(f)   =      Int_a^b f^2 dy            1/(b-a) Int_a^b f^2 dy
// L2Dist2(f,g) =      Int_a^b (f-g)^2 dy        1/(b-a) Int_a^b (f-g)^2 dy
// L2Norm(f)    = sqrt(Int_a^b f^2 dy)      sqrt(1/(b-a) Int_a^b f^2 dy)
// L2Dist(f,g)  = sqrt(Int_a^b (f-g)^2 dy)  sqrt(1/(b-a) Int_a^b (f-g)^2 dy
// L2IP...(f,g) =      Int_a^b f g dy            1/(b-a) Int_a^b f g dy
//

cfbasics::Real L2Norm2(const ChebyCoeff& f, bool normalize = true);
cfbasics::Real L2Dist2(const ChebyCoeff& f, const ChebyCoeff& g, bool normalize = true);
cfbasics::Real L2Norm(const ChebyCoeff& f, bool normalize = true);
cfbasics::Real L2Dist(const ChebyCoeff& f, const ChebyCoeff& g, bool normalize = true);
cfbasics::Real L2InnerProduct(const ChebyCoeff& f, const ChebyCoeff& g, bool normalize = true);

cfbasics::Real chebyNorm2(const ChebyCoeff& f, bool normalize = true);
cfbasics::Real chebyDist2(const ChebyCoeff& f, const ChebyCoeff& g, bool normalize = true);
cfbasics::Real chebyNorm(const ChebyCoeff& f, bool normalize = true);
cfbasics::Real chebyDist(const ChebyCoeff& f, const ChebyCoeff& g, bool normalize = true);
cfbasics::Real chebyInnerProduct(const ChebyCoeff& f, const ChebyCoeff& g, bool normalize = true);

cfbasics::Real norm2(const ChebyCoeff& f, NormType n, bool normalize = true);

cfbasics::Real norm(const ChebyCoeff& f, NormType n, bool normalize = true);
cfbasics::Real dist(const ChebyCoeff& f, const ChebyCoeff& g, NormType n, bool normalize = true);
cfbasics::Real innerProduct(const ChebyCoeff& f, const ChebyCoeff& g, NormType n, bool normalize = true);

cfbasics::Real L1Norm(const ChebyCoeff& f, bool normalize = true);
cfbasics::Real L1Dist(const ChebyCoeff& f, const ChebyCoeff& g, bool normalize = true);

cfbasics::Real LinfNorm(const ChebyCoeff& f);
cfbasics::Real LinfDist(const ChebyCoeff& f, const ChebyCoeff& g);

ComplexChebyCoeff operator*(cfbasics::Real c, const ComplexChebyCoeff& g);
ComplexChebyCoeff operator+(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g);
ComplexChebyCoeff operator-(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g);
bool operator==(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g);
bool operator!=(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g);

void diff(const ComplexChebyCoeff& f, ComplexChebyCoeff& df);
void diff2(const ComplexChebyCoeff& f, ComplexChebyCoeff& d2f);
void diff2(const ComplexChebyCoeff& f, ComplexChebyCoeff& d2f, ComplexChebyCoeff& tmp);
void diff(const ComplexChebyCoeff& f, ComplexChebyCoeff& d2f, int n);

ComplexChebyCoeff diff(const ComplexChebyCoeff& f);
ComplexChebyCoeff diff2(const ComplexChebyCoeff& f);
ComplexChebyCoeff diff(const ComplexChebyCoeff& f, int n);

void integrate(const ComplexChebyCoeff& df, ComplexChebyCoeff& f);
ComplexChebyCoeff integrate(const ComplexChebyCoeff& df);

std::ostream& operator<<(std::ostream& os, const ComplexChebyCoeff& f);

inline ChebyCoeff& Re(ComplexChebyCoeff& f);
inline ChebyCoeff& Im(ComplexChebyCoeff& f);
inline const ChebyCoeff& Re(const ComplexChebyCoeff& f);
inline const ChebyCoeff& Im(const ComplexChebyCoeff& f);

// L1Norms are simple sums of |u_n| with no normalization whatsoever.
// Real L1Norm(const ComplexChebyCoeff& f);
// Real L1Dist(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g);

//   normalize  ?           false                   true (c= 1/(b-a))
// L2Norm2(f)   =      Int_a^b f f* dy         1/(b-a) Int_a^b f f* dy
// L2Dist2(f,g) =      Int_a^b (f-g)(f-g)* dy  1/(b-a) Int_a^b (f-g)(f-g)* dy
// L2Norm(f)    = sqrt(Int_a^b f f* dy)         sqrt(c Int_a^b f f* dy)
// L2Dist(f,g)  = sqrt(Int_a^b (f-g)(f-g)* dy)  sqrt(c Int_a^b (f-g)(f-g)* dy
// L2IP...(f,g) =      Int_a^b f g dy          1/(b-a) Int_a^b f g* dy
//
cfbasics::Real L2Norm2(const ComplexChebyCoeff& f, bool normalize = true);
cfbasics::Real L2Dist2(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, bool normalize = true);
cfbasics::Real L2Norm(const ComplexChebyCoeff& f, bool normalize = true);
cfbasics::Real L2Dist(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, bool normalize = true);
cfbasics::Complex L2InnerProduct(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, bool normalize = true);

cfbasics::Real chebyNorm2(const ComplexChebyCoeff& f, bool normalize = true);
cfbasics::Real chebyDist2(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, bool normalize = true);
cfbasics::Real chebyNorm(const ComplexChebyCoeff& f, bool normalize = true);
cfbasics::Real chebyDist(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, bool normalize = true);
cfbasics::Complex chebyInnerProduct(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, bool normalize = true);

cfbasics::Real norm2(const ComplexChebyCoeff& f, NormType n, bool normalize = true);

cfbasics::Real norm(const ComplexChebyCoeff& f, NormType n, bool normalize = true);
cfbasics::Real dist(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, NormType n, bool normalize = true);
cfbasics::Complex innerProduct(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, NormType n,
                               bool normalize = true);

cfbasics::Real L1Norm(const ComplexChebyCoeff& f, bool normalize = true);
cfbasics::Real L1Dist(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g, bool normalize = true);

cfbasics::Real LinfNorm(const ComplexChebyCoeff& f);
cfbasics::Real LinfDist(const ComplexChebyCoeff& f, const ComplexChebyCoeff& g);

inline cfbasics::Real ChebyCoeff::a() const { return a_; }
inline cfbasics::Real ChebyCoeff::b() const { return b_; }
inline cfbasics::Real ChebyCoeff::L() const { return b_ - a_; }
inline int ChebyCoeff::N() const { return data_.size(); }
inline int ChebyCoeff::numModes() const { return data_.size(); }
inline cfbasics::fieldstate ChebyCoeff::state() const { return state_; }

inline cfbasics::Real ComplexChebyCoeff::a() const { return re.a(); }
inline cfbasics::Real ComplexChebyCoeff::b() const { return re.b(); }
inline cfbasics::Real ComplexChebyCoeff::L() const { return re.L(); }
inline int ComplexChebyCoeff::N() const { return re.N(); }
inline int ComplexChebyCoeff::length() const { return re.length(); }
inline int ComplexChebyCoeff::numModes() const { return re.numModes(); }
inline cfbasics::fieldstate ComplexChebyCoeff::state() const { return re.state(); }

inline ChebyCoeff& Re(ComplexChebyCoeff& f) { return f.re; }
inline ChebyCoeff& Im(ComplexChebyCoeff& f) { return f.im; }
inline const ChebyCoeff& Re(const ComplexChebyCoeff& f) { return f.re; }
inline const ChebyCoeff& Im(const ComplexChebyCoeff& f) { return f.im; }

inline cfbasics::Complex ComplexChebyCoeff::operator[](int i) const { return re[i] + cfbasics::I * im[i]; }
inline void ComplexChebyCoeff::set(int i, cfbasics::Complex c) {
    re[i] = cfbasics::Re(c);
    im[i] = cfbasics::Im(c);
}
inline void ComplexChebyCoeff::add(int i, cfbasics::Complex c) {
    re[i] += cfbasics::Re(c);
    im[i] += cfbasics::Im(c);
}
inline void ComplexChebyCoeff::sub(int i, cfbasics::Complex c) {
    re[i] -= cfbasics::Re(c);
    im[i] -= cfbasics::Im(c);
}

inline int ChebyTransform::N() const { return N_; }
inline int ChebyTransform::length() const { return N_; }

}  // namespace channelflow

#endif /* CHEBYSHEV_H */
