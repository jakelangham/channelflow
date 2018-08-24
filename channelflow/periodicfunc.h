/**
 * Fourier representation of periodic functions
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#ifndef CHANNELFLOW_PERIODICFUNC_H
#define CHANNELFLOW_PERIODICFUNC_H

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"

#include <fftw3.h>

#include <iostream>
#include <memory>

namespace channelflow {

cfbasics::Vector periodicpoints(int N, cfbasics::Real L);  // N+1 pts {0, L/N, 2L/N, ..., L}

// A class for representing 1d periodic functions on the interval [0, L]
class PeriodicFunc {
   public:
    PeriodicFunc();
    PeriodicFunc(const PeriodicFunc& f);
    PeriodicFunc(uint N, cfbasics::Real L, cfbasics::fieldstate s = cfbasics::Spectral, uint flag = FFTW_ESTIMATE);
    PeriodicFunc(const std::string& filebase);  // read ascii from file
    PeriodicFunc& operator=(const PeriodicFunc& f);

    inline cfbasics::Real operator()(uint n) const;  // gridpoint value f(x_n)
    inline cfbasics::Real& operator()(uint n);       // gridpoint value f(x_n)
    inline cfbasics::Complex cmplx(uint k) const;    // fourier coefficient f_k
    inline cfbasics::Complex& cmplx(uint k);         // fourier coefficient f_k
    inline cfbasics::Real operator[](uint n) const;  // nth data cfarray value
    inline cfbasics::Real& operator[](uint n);       // nth data cfarray value
    inline cfbasics::Real x(uint n) const;           // position of nth gridpt x_n

    void save(const std::string& filebase, cfbasics::fieldstate s = cfbasics::Physical) const;

    void reconfig(const PeriodicFunc& f);
    void randomize(cfbasics::Real magn, cfbasics::Real decay);
    void setLength(cfbasics::Real L);
    void setState(cfbasics::fieldstate s);
    void setToZero();
    void resize(uint N, cfbasics::Real L);

    cfbasics::Real operator()(cfbasics::Real x) const;
    cfbasics::Real eval(cfbasics::Real x) const;

    // void eval(const Vector& x, PeriodicFunc& g) const;
    // PeriodicFunc eval(const Vector& x) const;
    // void binaryDump(std::ostream& os) const;
    // void binaryLoad(std::istream& is);
    // void fill(const PeriodicFunc& g);
    // void interpolate(const PeriodicFunc& g); // *this is on subdomain of g.
    // void reflect(const PeriodicFunc& g, parity p); // reflect g about u.a()
    // inline uint Ngridpts() const; // same as N

    inline cfbasics::Real L() const;
    inline uint N() const;
    inline uint Npad() const;    // 2*(N/2+1)
    inline uint Nmodes() const;  // N/2+1
    inline uint kmax() const;    // N/2 (note: N/2 mode for N even is set to 0)

    inline cfbasics::fieldstate state() const;
    cfbasics::Real mean() const;

    PeriodicFunc& operator*=(cfbasics::Real c);
    PeriodicFunc& operator+=(const PeriodicFunc& g);
    PeriodicFunc& operator-=(const PeriodicFunc& g);
    PeriodicFunc& operator*=(const PeriodicFunc& g);  // dottimes, only for Physical

    void fft();                              // transform from Physical to Spectral
    void ifft();                             // transform from Spectral to Physcial
    void makeSpectral();                     // if Physical, transform to Spectral
    void makePhysical();                     // if Spectral, transform to Physical
    void makeState(cfbasics::fieldstate s);  // if state != s, transform to state s

    bool congruent(const PeriodicFunc& g) const;

    friend void swap(PeriodicFunc& f, PeriodicFunc& g);

   private:
    uint N_ = 0;                                          // number of gridpoints
    cfbasics::Real L_ = 0;                                // upper bound of domain
    std::unique_ptr<void, void (*)(void*)> data_handle_;  // Manages the lifetime of FFTW data buffer
    cfbasics::Real* rdata_ = nullptr;                     // cfarray for Fourier coeffs or physical data
    cfbasics::Complex* cdata_ = nullptr;                  // Complex alias for rdata_
    cfbasics::fieldstate state_;                          // indicates Physical or Spectral state of object

    void fftw_initialize(uint fftw_flags = FFTW_ESTIMATE);
    uint fftw_flags_;
    std::unique_ptr<std::remove_pointer<fftw_plan>::type, void (*)(fftw_plan)> forward_plan_;
    std::unique_ptr<std::remove_pointer<fftw_plan>::type, void (*)(fftw_plan)> inverse_plan_;
};

PeriodicFunc operator*(cfbasics::Real c, const PeriodicFunc& g);
PeriodicFunc operator+(const PeriodicFunc& f, const PeriodicFunc& g);
PeriodicFunc operator-(const PeriodicFunc& f, const PeriodicFunc& g);
bool operator==(const PeriodicFunc& f, const PeriodicFunc& g);
bool operator!=(const PeriodicFunc& f, const PeriodicFunc& g);

void diff(const PeriodicFunc& f, PeriodicFunc& df);
void diff2(const PeriodicFunc& f, PeriodicFunc& d2f);
void diff2(const PeriodicFunc& f, PeriodicFunc& d2f, PeriodicFunc& tmp);
void diff(const PeriodicFunc& f, PeriodicFunc& df, uint n);
PeriodicFunc diff(const PeriodicFunc& f);
PeriodicFunc diff2(const PeriodicFunc& f);
PeriodicFunc diff(const PeriodicFunc& f, uint n);

// Integrate sets the arbitrary const of integration so that mean(u)==0.
void integrate(const PeriodicFunc& df, PeriodicFunc& f);
PeriodicFunc integrate(const PeriodicFunc& df);

// Find a zero of f(x) via Newton search.
cfbasics::Real newtonSearch(const PeriodicFunc& f, cfbasics::Real xguess, int Nmax = 10, cfbasics::Real eps = 1e-12);

//   normalize  ?   false                  :  true
// L2Norm2(f)   =      Int_a^b f^2 dy            1/(b-a) Int_a^b f^2 dy
// L2Dist2(f,g) =      Int_a^b (f-g)^2 dy        1/(b-a) Int_a^b (f-g)^2 dy
// L2Norm(f)    = sqrt(Int_a^b f^2 dy)      sqrt(1/(b-a) Int_a^b f^2 dy)
// L2Dist(f,g)  = sqrt(Int_a^b (f-g)^2 dy)  sqrt(1/(b-a) Int_a^b (f-g)^2 dy
// L2IP...(f,g) =      Int_a^b f g dy            1/(b-a) Int_a^b f g dy
//
// L2IP is L2 Inner Product

cfbasics::Real L2Norm(const PeriodicFunc& f, bool normalize = true);
cfbasics::Real L2Norm2(const PeriodicFunc& f, bool normalize = true);
cfbasics::Real L2Dist(const PeriodicFunc& f, const PeriodicFunc& g, bool normalize = true);
cfbasics::Real L2Dist2(const PeriodicFunc& f, const PeriodicFunc& g, bool normalize = true);
cfbasics::Real L2IP(const PeriodicFunc& f, const PeriodicFunc& g, bool normalize = true);

// Vector zeros(const PeriodicFunc& f, Real epsilon=1e-14);

std::ostream& operator<<(std::ostream& os, const PeriodicFunc& f);

inline cfbasics::Real PeriodicFunc::L() const { return L_; }
inline uint PeriodicFunc::N() const { return N_; }
inline uint PeriodicFunc::Npad() const { return 2 * (N_ / 2 + 1); }
inline uint PeriodicFunc::Nmodes() const { return (N_ / 2 + 1); }
inline uint PeriodicFunc::kmax() const { return N_ / 2; }
inline cfbasics::fieldstate PeriodicFunc::state() const { return state_; }

inline cfbasics::Real PeriodicFunc::operator()(uint n) const {
    assert(state_ == cfbasics::Physical);
    assert(n < N_);
    return rdata_[n];
}

inline cfbasics::Real& PeriodicFunc::operator()(uint n) {
    assert(state_ == cfbasics::Physical);
    assert(n < N_);
    return rdata_[n];
}

inline cfbasics::Complex PeriodicFunc::cmplx(uint k) const {
    assert(state_ == cfbasics::Spectral);
    assert(k < Nmodes());
    return cdata_[k];
}

inline cfbasics::Complex& PeriodicFunc::cmplx(uint k) {
    assert(state_ == cfbasics::Spectral);
    assert(k < Nmodes());
    return cdata_[k];
}

inline cfbasics::Real PeriodicFunc::operator[](uint n) const {
    assert(n < Npad());
    return rdata_[n];
}

inline cfbasics::Real& PeriodicFunc::operator[](uint n) {
    assert(n < Npad());
    return rdata_[n];
}

inline cfbasics::Real PeriodicFunc::x(uint n) const { return n * (L_ / N_); }

}  // namespace channelflow

#endif
