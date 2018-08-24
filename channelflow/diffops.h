/**
 * Differential operators for FlowFields
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#ifndef CHANNELFLOW_DIFFOPS_H
#define CHANNELFLOW_DIFFOPS_H

#include <vector>

#include "cfbasics/cfbasics.h"
#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"
#include "channelflow/flowfield.h"

namespace channelflow {

// innerProduct(f,g) == Int f g*  dx dy dz   (/(Lx Ly Lz) if normalize==true)
// L2Norm2(f)  == Int ||f||^2  dx dy dz      (/(Lx Ly Lz) if normalize==true)
// bcNorm2(f)  == Int ||f||^2  dx dz at y=a,b   (/(Lx Lz) if normalize==true)
// divNorm2(f) == L2Norm2(div(u))
//
// L2Norm(f)    == sqrt(L2Norm2(f))
// bcNorm(f)    == sqrt(bcNorm2(f))
// L2Dist(f,g)  == sqrt(L2Dist2(f,g))
// L2Dist2(f,g) == L2Norm2(f-g)
// etc.

cfbasics::Real L1Norm(const FlowField& f, bool normalize = true);
cfbasics::Real L1Dist(const FlowField& f, const FlowField& g, bool normalize = true);
cfbasics::Real LinfNorm(const FlowField& f);
cfbasics::Real LinfDist(const FlowField& f, const FlowField& g);

cfbasics::Real L2Norm(const FlowField& f, bool normalize = true);
cfbasics::Real L2Norm2(const FlowField& f, bool normalize = true);
cfbasics::Real L2Dist(const FlowField& f, const FlowField& g, bool normalize = true);
cfbasics::Real L2Dist2(const FlowField& f, const FlowField& g, bool normalize = true);

cfbasics::Real chebyNorm(const FlowField& f, bool normalize = true);
cfbasics::Real chebyNorm2(const FlowField& f, bool normalize = true);
cfbasics::Real chebyDist(const FlowField& f, const FlowField& g, bool normalize = true);
cfbasics::Real chebyDist2(const FlowField& f, const FlowField& g, bool normalize = true);

cfbasics::Real bcNorm(const FlowField& f, bool normalize = true);
cfbasics::Real bcNorm2(const FlowField& f, bool normalize = true);
cfbasics::Real bcDist(const FlowField& f, const FlowField& g, bool normalize = true);
cfbasics::Real bcDist2(const FlowField& f, const FlowField& g, bool normalize = true);

cfbasics::Real divNorm(const FlowField& f, bool normalize = true);
cfbasics::Real divNorm2(const FlowField& f, bool normalize = true);
cfbasics::Real divDist(const FlowField& f, const FlowField& g, bool normalize = true);
cfbasics::Real divDist2(const FlowField& f, const FlowField& g, bool normalize = true);

// Restrict sums in norm computation to |kx|<=kxmax, |kz|<=kzmax,
cfbasics::Real L2Norm(const FlowField& f, int kxmax, int kzmax, bool normalize = true);
cfbasics::Real L2Norm2(const FlowField& f, int kxmax, int kzmax, bool normalize = true);
cfbasics::Real L2Dist(const FlowField& f, const FlowField& g, int kxmax, int kzmax, bool normalize = true);
cfbasics::Real L2Dist2(const FlowField& f, const FlowField& g, int kxmax, int kzmax, bool normalize = true);
cfbasics::Real L2InnerProduct(const FlowField& f, const FlowField& g, int kxmax, int kzmax, bool normalize = true);

// Vector  vectorL1Dist(const FlowField& f, const FlowField& v);
// Complex innerProduct(const FlowField& f, const FlowField& v);
cfbasics::Real L2InnerProduct(const FlowField& f, const FlowField& g, bool normalize = true);

cfbasics::Real dissipation(const FlowField& f, bool normalize = true);     // 1/(LxLyLz) int_V |curl u|^2 dV
cfbasics::Real wallshear(const FlowField& f, bool normalize = true);       // |wallshearLower| + |wallshearUpper|
cfbasics::Real wallshearLower(const FlowField& f, bool normalize = true);  // 1/(2LxLz)  int_lower_wall du/dy dx dz
cfbasics::Real wallshearUpper(const FlowField& f, bool normalize = true);  // 1/(2LxLz)  int_upper_wall du/dy dx dz

cfbasics::Real L2Norm2_3d(const FlowField& f, bool normalize = true);  // L2Norm2 of all kx!=0 modes
cfbasics::Real L2Norm3d(const FlowField& f, bool normalize = true);    // L2Norm  of all kx!=0 modes

cfbasics::cfarray<cfbasics::Real> truncerr(const FlowField& f);

cfbasics::Real min_x_L2Dist(const FlowField& u0, const FlowField& u1, cfbasics::Real tol = 1e-3);

// ====================================================================
// Note: these were previously defaulted to normalize=false 2006-03-02
cfbasics::Complex L2InnerProduct(const FlowField& f, const BasisFunc& phi, bool normalize = true);

cfbasics::Real L2InnerProduct(const FlowField& f, const RealProfile& phi, bool normalize = true);

inline cfbasics::Real L2IP(const FlowField& f, const FlowField& g, bool normalize = true) {
    return L2InnerProduct(f, g, normalize);
}

inline cfbasics::Real L2IP(const FlowField& f, const FlowField& g, int kxmax, int kzmax, bool normalize = true) {
    return L2InnerProduct(f, g, kxmax, kzmax, normalize);
}

inline cfbasics::Complex L2IP(const FlowField& f, const BasisFunc& phi, bool normalize = true) {
    return L2InnerProduct(f, phi, normalize);
}

inline cfbasics::Real L2IP(const FlowField& f, const RealProfile& phi, bool normalize = true) {
    return L2InnerProduct(f, phi, normalize);
}

// ====================================================================
cfbasics::Real L2InnerProduct(const RealProfileNG& e, const FlowField& u, bool normalize = true);

inline cfbasics::Real L2InnerProduct(const FlowField& u, const RealProfileNG& e, bool normalize = true) {
    return L2InnerProduct(e, u, normalize);
}

inline cfbasics::Real L2IP(const RealProfileNG& e, const FlowField& u, bool normalize = true) {
    return L2InnerProduct(e, u, normalize);
}
inline cfbasics::Real L2IP(const FlowField& u, const RealProfileNG& e, bool normalize = true) {
    return L2InnerProduct(e, u, normalize);
}

// The following functions assume an orthonormal basis.
void field2coeff(const std::vector<RealProfile>& basis, const FlowField& u, cfbasics::cfarray<cfbasics::Real>& a);
void coeff2field(const std::vector<RealProfile>& basis, const cfbasics::cfarray<cfbasics::Real>& a, FlowField& u);

void field2coeff(const std::vector<BasisFunc>& basis, const FlowField& u, cfbasics::cfarray<cfbasics::Complex>& a);
void coeff2field(const std::vector<BasisFunc>& basis, const cfbasics::cfarray<cfbasics::Complex>& a, FlowField& u);

void field2coeff(const std::vector<RealProfileNG>& basis, const FlowField& u, std::vector<cfbasics::Real>& a);
void coeff2field(const std::vector<RealProfileNG>& basis, const std::vector<cfbasics::Real>& a, FlowField& u);

void swap(FlowField& f, FlowField& v);  // exchange data of two flow fields.

// In general, these fucntions take input fields in any state and return
// output fields in Spectral,Spectral leaving inputs unchanged.
// But it's most efficient to send input fields as Spectral,Spectral.
// Input fields are logically constant, but they may be modified and
// unmodified (transforms and inverses) during the course of the computation.

//  lapl(f,laplf) computes laplf(nx,ny,nz,i) = d^2 f(nx,ny,nz,i)/dx_j^2
//  grad(f,gradf) computes
//                gradf(nx,ny,nz,i)        = df(nx,ny,nz,0)/dx_i for 1d f
//                gradf(nx,ny,nz,i3j(i,j)) = df(nx,ny,nz,i)/dx_j for 3d f

// Warning: these functions are convenient but inefficient due to
// construction and copying of returned FlowFields. Prefer the
// void-returning versions below for repeated calculations.
FlowField xdiff(const FlowField& f, int n = 1);              // d^nf/dx^n
FlowField ydiff(const FlowField& f, int n = 1);              // d^nf/dy^n
FlowField zdiff(const FlowField& f, int n = 1);              // d^nf/dz^n
FlowField diff(const FlowField& f, int i, int n);            // d^nf/dx[i]^n
FlowField diff(const FlowField& f, int nx, int ny, int nz);  // d^(nx+ny+nz)/(dx^nx dy^ny dz^nz) f

FlowField grad(const FlowField& f);   // grad f
FlowField lapl(const FlowField& f);   // lapl f
FlowField curl(const FlowField& f);   // del cross f
FlowField curl(const FlowField& f);   // del cross f
FlowField norm(const FlowField& f);   // ||f||
FlowField norm2(const FlowField& f);  // ||f||^2
FlowField div(const FlowField& f);    // del dot f

FlowField cross(const FlowField& f, const FlowField& g);
FlowField outer(const FlowField& f, const FlowField& g);  // fi gj ei ej
FlowField dot(const FlowField& f, const FlowField& g);
FlowField energy(const FlowField& u);
FlowField energy(const FlowField& u, ChebyCoeff& U);

// The following are some fluid dynamics measures of vorticity
// Q = 1/2 (||A||^2 - ||S||^2) pointwise
//     with A,S are antisymm/symm parts of grad u
// swirling = imaginary part of eig(grad v) pointwise
FlowField Qcriterion(const FlowField& u);
FlowField xyavg(FlowField& u);

// These functions are preferred
void xdiff(const FlowField& f, FlowField& dfdx, int n = 1);
void ydiff(const FlowField& f, FlowField& dfdy, int n = 1);
void zdiff(const FlowField& f, FlowField& dfdz, int n = 1);
void diff(const FlowField& f, FlowField& df, int i, int n);
void diff(const FlowField& f, FlowField& df, int nx, int ny, int nz);

void ydiffOld(const FlowField& f, FlowField& dfdy, int n = 1);
void ydiffOnce(const FlowField& f, FlowField& dfdy);

void grad(const FlowField& f, FlowField& grad_f);    // grad f, 1d->3d or 3d->9d
void d_dx(const FlowField& f, FlowField& ddx_f);     // SlicesHack d_dz f, 3d
void d_dz(const FlowField& f, FlowField& ddz_f);     // SlicesHack d_dz f, 3d
void lapl(const FlowField& f, FlowField& lapl_f);    // lapl f
void curl(const FlowField& f, FlowField& curl_f);    // del cross f
void norm(const FlowField& f, FlowField& norm_f);    // ||f||
void norm2(const FlowField& f, FlowField& norm2_f);  // ||f||^2
void div(const FlowField& f, FlowField& divf, const cfbasics::fieldstate finalstate = cfbasics::Spectral);  // del dot f

void cross(const FlowField& f, const FlowField& g, FlowField& f_cross_g,
           const cfbasics::fieldstate finalstate = cfbasics::Spectral);
void outer(const FlowField& f, const FlowField& g, FlowField& fg);  // fi gj ei ej
void dot(const FlowField& f, const FlowField& g, FlowField& f_dot_g);
void energy(const FlowField& u, FlowField& e);
void energy(const FlowField& u, const ChebyCoeff& U, FlowField& e);
void Qcriterion(const FlowField& u, FlowField& Q);

// Different ways to calculate nonlinear term of Navier-Stokes equation.
// Let u' = u + U ex
// Rotational:  f = (curl u)  cross u + U du/dx + v dUdy ex, tmp = curl u
// Rotational2: f = (curl u') cross u'                       tmp = curl u'
// Convection:  f = u' dot grad u',                          tmp = grad u'
// Divergence:  f = div (u' u'),                             tmp = u' u'
// Skewsymm  :  f = 1/2 [u' dot (grad u') + div (u' u')],    tmp = grad u'
// Alternating:  alternate btwn divergence and convection
// Alternating_: alternate btwn convection and divergence
// Linearized:  f = u dot grad (U ex) + (U ex) dot grad u
//                = v dU/dy ex + U du/dx
//
// LinearizedAboutField:
//            utot = u + U = (u0+U) + (u-u0) = ubase + du
//              f = (curl ubase) cross ubase
//                  (curl ubase) cross du
//                  (curl du)    cross ubase
//                = (curl utot) cross (utot) - (curl du) cross (du)

// Note: these functions return f in Spectral, Spectral and u,U in entry states
// The functions are most efficient if u,U enter in Spectral, Spectral.
// If you want dealiased nonlinearities, call u.setAliasedModesToZero()
// beforehand and possibly f.setAliasedModesToZero() afterwards.

// Rotational f = (curl u) x u
// Convection f = u dot grad u
// Divergence f = div (u u)
// Skew-Symm  f = 1/2 [u dot grad u + div (u u)]
// Linearized f = dU/dy v ex + du/dx U
// Linearized f = (curl ubase)   x u + (curl u) x  ubase

void rotationalNL(const FlowField& u, FlowField& f, FlowField& tmp,
                  const cfbasics::fieldstate finalstate = cfbasics::Spectral);
void convectionNL(const FlowField& u, FlowField& f, FlowField& tmp,
                  const cfbasics::fieldstate finalstate = cfbasics::Spectral);
void divergenceNL(const FlowField& u, FlowField& f, FlowField& tmp,
                  const cfbasics::fieldstate finalstate = cfbasics::Spectral);
void skewsymmetricNL(const FlowField& u, FlowField& f, FlowField& tmp,
                     const cfbasics::fieldstate finalstate = cfbasics::Spectral);
void linearizedNL(const FlowField& u, const ChebyCoeff& U, const ChebyCoeff& Wbase, FlowField& f,
                  const cfbasics::fieldstate finalstate = cfbasics::Spectral);

// Adjoint terms
// u+ is the adjoint velocity
// u  is the direct velocity
// U  is the base flow

// adjointTerms f = + (u + U) grad u+ - u+ grad^T (u + U)		tmp = grad u	tmpadj = grad^T(u + U)

void adjointTerms(const FlowField& u, FlowField& u_dir, FlowField& f, FlowField& tmp, FlowField& tmpadj,
                  const cfbasics::fieldstate finalstate = cfbasics::Spectral);

// Perturbation terms
// u+ is the adjoint velocity
// u  is the direct velocity
// U  is the base flow

// perturbationTermsLin f = + u grad (u_ecs + U) + (u_ecs + U) grad u 		tmp = grad u	tmppert = grad(u_ecs +
// U)

void perturbationTermsLin(const FlowField& u, FlowField& u_nlin, FlowField& f, FlowField& tmp, FlowField& tmppert,
                          const cfbasics::fieldstate finalstate = cfbasics::Spectral);

// perturbationTermsNLin f = + u grad u + u grad (u_ecs + U) + (u_ecs + U) grad u 		tmp = grad u	tmppert
// = grad(u_ecs + U)

void perturbationTermsNLin(const FlowField& u, FlowField& u_nlin, FlowField& f, FlowField& tmp, FlowField& tmppert,
                           const cfbasics::fieldstate finalstate = cfbasics::Spectral);

void linearAboutFieldNL(const FlowField& u, const FlowField& ubase, const ChebyCoeff& U, FlowField& f, FlowField& tmp,
                        FlowField& tmp2, const cfbasics::fieldstate finalstate = cfbasics::Spectral);

void linearizedNL(const FlowField& u_, const FlowField& ubtot, const FlowField& grad_ubtot, FlowField& f,
                  FlowField& tmp, const cfbasics::fieldstate finalstate = cfbasics::Spectral);

void assignOrrSommField(FlowField& u, FlowField& P, cfbasics::Real t, cfbasics::Real Reynolds, cfbasics::Complex omega,
                        const ComplexChebyCoeff& ueig, const ComplexChebyCoeff& veig, const ComplexChebyCoeff& peig);

void dotgrad(const FlowField& u, const FlowField& v, FlowField& u_dotgrad_v, FlowField& tmp);
void dotgradScalar(const FlowField& u, const FlowField& s, FlowField& u_dotgrad_s, FlowField& tmp);
FlowField dotgrad(const FlowField& u, const FlowField& v, FlowField& tmp);

void randomUprofile(ComplexChebyCoeff& u, cfbasics::Real mag, cfbasics::Real spectralDecay);
void randomVprofile(ComplexChebyCoeff& v, cfbasics::Real mag, cfbasics::Real spectralDecay);
void randomProfile(ComplexChebyCoeff& u, ComplexChebyCoeff& v, ComplexChebyCoeff& w, int kx, int kz, cfbasics::Real Lx,
                   cfbasics::Real Lz, cfbasics::Real mag, cfbasics::Real spectralDecay);

void chebyUprofile(ComplexChebyCoeff& u, int n, cfbasics::Real decay);
void chebyVprofile(ComplexChebyCoeff& v, int n, cfbasics::Real decay);
void chebyProfile(ComplexChebyCoeff& u, ComplexChebyCoeff& v, ComplexChebyCoeff& w, int un, int vn, int kx, int kz,
                  cfbasics::Real Lx, cfbasics::Real Lz, cfbasics::Real decay);

FlowField extractRolls(const FlowField& u);

cfbasics::Real Ecf(const FlowField& u);
cfbasics::Real L2Norm_uvw(const FlowField& u, const bool ux, const bool uy, const bool uz);

cfbasics::Real getdPdx(const FlowField& u, cfbasics::Real nu);
cfbasics::Real getdPdz(const FlowField& u, cfbasics::Real nu);
cfbasics::Real getUbulk(const FlowField& u);
cfbasics::Real getWbulk(const FlowField& u);

std::string fieldstats_t(const FlowField& u, cfbasics::Real t);  // return some information about u
std::string fieldstatsheader_t(const std::string tname = "t");   // header for fieldstats
std::string fieldstats(const FlowField& u);                      // return some information about u
std::string fieldstatsheader();                                  // header for fieldstats

}  // namespace channelflow
#endif
