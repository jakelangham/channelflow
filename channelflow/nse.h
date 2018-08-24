/**
 * System class of Navier-Stokes Equation for standard channel flows
 *
 * Encapsulated code from DNSAlgo classes (FS)
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#ifndef CHANNELFLOW_NSE_H
#define CHANNELFLOW_NSE_H

#include "channelflow/diffops.h"
#include "channelflow/dnsflags.h"
#include "channelflow/flowfield.h"
#include "channelflow/tausolver.h"

namespace channelflow {

void navierstokesNL(const channelflow::FlowField& u, channelflow::ChebyCoeff Ubase, channelflow::ChebyCoeff Wbase,
                    channelflow::FlowField& f, channelflow::FlowField& tmp, channelflow::DNSFlags& flags);

class NSE {
   public:
    NSE();
    NSE(const NSE& nse);
    NSE(const std::vector<FlowField>& fields, const DNSFlags& flags);
    NSE(const std::vector<FlowField>& fields, const std::vector<ChebyCoeff>& base, const DNSFlags& flags);
    virtual ~NSE();

    virtual void nonlinear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields);
    virtual void linear(const std::vector<FlowField>& infields, std::vector<FlowField>& outfields);
    // calls a tausolver for each Fourier mode
    virtual void solve(std::vector<FlowField>& outfields, const std::vector<FlowField>& infields, const int i = 0);

    // redefines the tausolver objects with new time-stepping constant (allocates memory for tausolver at first use)
    virtual void reset_lambda(const std::vector<cfbasics::Real> lambda_t);

    // vector of RHS is smaller than of fields because of missing pressure equation
    virtual std::vector<FlowField> createRHS(const std::vector<FlowField>& fields) const;

    // returns vector of symmetries confining the vector of fields to a subspace
    virtual std::vector<cfbasics::cfarray<FieldSymmetry>> createSymmVec() const;

    inline int taskid() const;

    void reset_gradp(cfbasics::Real dPdx, cfbasics::Real dPdz);    // change dPdx and enforce const dPdx
    void reset_bulkv(cfbasics::Real Ubulk, cfbasics::Real Wbulk);  // change Ubulk and enforce const Ubulk

    // getter functions
    int Nx() const;
    int Ny() const;
    int Nz() const;

    cfbasics::Real Lx() const;
    cfbasics::Real Lz() const;
    cfbasics::Real a() const;
    cfbasics::Real b() const;
    cfbasics::Real nu() const;

    cfbasics::Real dPdx() const;  // the mean pressure gradient at the current time
    cfbasics::Real dPdz() const;
    cfbasics::Real Ubulk() const;  // the actual bulk velocity at the current time
    cfbasics::Real Wbulk() const;
    cfbasics::Real dPdxRef() const;  // the mean press grad enforced during integration
    cfbasics::Real dPdzRef() const;
    cfbasics::Real UbulkRef() const;  // the bulk velocity enforced during integ.
    cfbasics::Real WbulkRef() const;

    virtual const ChebyCoeff& Ubase() const;
    virtual const ChebyCoeff& Wbase() const;

   protected:
    std::vector<cfbasics::Real> lambda_t_;
    TauSolver*** tausolver_;  // 3d cfarray of tausolvers, indexed by [i][mx][mz]

    DNSFlags flags_;  // User-defined integration parameters
    int taskid_;

    // Spatial parameter members
    cfbasics::lint nxlocmin_;
    cfbasics::lint Nxloc_;
    cfbasics::lint nylocmin_;
    cfbasics::lint nylocmax_;
    int Nz_;  // number of Z gridpoints
    cfbasics::lint mxlocmin_;
    cfbasics::lint Mxloc_;
    cfbasics::lint My_;
    cfbasics::lint mzlocmin_;
    cfbasics::lint Mzloc_;
    int Nyd_;      // number of dealiased Chebyshev T(y) modes
    int kxd_max_;  // maximum value of kx among dealiased modes
    int kzd_max_;  // maximum value of kz among dealiased modes
    cfbasics::Real Lx_;
    cfbasics::Real Lz_;
    cfbasics::Real a_;
    cfbasics::Real b_;
    int kxmax_;
    int kzmax_;
    std::vector<int> kxloc_;
    std::vector<int> kzloc_;

    // Base flow members
    cfbasics::Real dPdxRef_;    // Enforced mean pressure gradient (0.0 if unused).
    cfbasics::Real dPdxAct_;    // Actual   mean pressure gradient at previous timestep.
    cfbasics::Real dPdzRef_;    //
    cfbasics::Real dPdzAct_;    //
    cfbasics::Real UbulkRef_;   // Enforced total bulk velocity (0.0 if unused).
    cfbasics::Real UbulkAct_;   // Actual total bulk velocity bulk obtained.
    cfbasics::Real UbulkBase_;  // Bulk velocity of Ubase
    cfbasics::Real WbulkRef_;
    cfbasics::Real WbulkAct_;
    cfbasics::Real WbulkBase_;

    ChebyCoeff Ubase_;    // baseflow physical
    ChebyCoeff Ubaseyy_;  // baseflow'' physical
    ChebyCoeff Wbase_;    // baseflow physical
    ChebyCoeff Wbaseyy_;  // baseflow'' physical

    // Memspace members
    FlowField tmp_;  // tmp space for nonlinearity calculation
    // These variables are used as temp storage when solving indpt tau problems.
    ComplexChebyCoeff uk_;  // profile of u_{kx,kz} (y) at t = n dt
    ComplexChebyCoeff vk_;
    ComplexChebyCoeff wk_;
    ComplexChebyCoeff Pk_;   // profile of P_{kx,kz} (y)
    ComplexChebyCoeff Pyk_;  // profile of dP_{kx,kz}/dy (y)
    ComplexChebyCoeff Ruk_;
    ComplexChebyCoeff Rvk_;
    ComplexChebyCoeff Rwk_;

    int kxmaxDealiased() const;
    int kzmaxDealiased() const;
    bool isAliasedMode(int kx, int kz) const;

    void init(FlowField& u);  // common constructor code

   private:
    void createCFBaseFlow();                    // method called only at construction
    void initCFConstraint(const FlowField& u);  // method called only at construction
};

inline int NSE::taskid() const { return taskid_; }

/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Given a baseflow, fluctation, modified pressure triple (U0,u0,q0) and
// a new baseflow U1, compute new fluctuation u1 and modified pressure q1.
//     utot == U0 + u0 == U1 + u1
// pressure == q0 - 1/2 || u0 ||^2 == q1 - 1/2 || u1 ||^2
void changeBaseFlow(const ChebyCoeff& U0, const FlowField& u0, const FlowField& q0, const ChebyCoeff& U1, FlowField& u1,
                    FlowField& q1);

// Construct laminar flow profile for given flow parameters.
// [a,b]   == y position of [lower, upper] walls
// [ua,ub] == in-plane speed of [lower, upper] walls
// constraint == is mean pressure gradient fixed, or mean (bulk) velocity?
// dPdx, Ubulk == value of fixed pressure gradient or fixed Ubulk velocity
// Vsuck == suction velocity at walls (asymptotic suction boundary layer)
ChebyCoeff laminarProfile(cfbasics::Real nu, MeanConstraint constraint, cfbasics::Real dPdx, cfbasics::Real Ubulk,
                          cfbasics::Real Vsuck, cfbasics::Real a, cfbasics::Real b, cfbasics::Real ua,
                          cfbasics::Real ub, int Ny);

ChebyCoeff laminarProfile(const DNSFlags& flags, cfbasics::Real a, cfbasics::Real b, int Ny);

// Return viscosity nu == Uh/Reynolds, where U is determined by VelocityScale.
// For BulkScale or ParabolicScale, U is determined from either dPdx or Ubulk,
// depending on whether MeanConstraint is PressureGradient or BulkVelocity
cfbasics::Real viscosity(cfbasics::Real Reynolds, VelocityScale vscale, MeanConstraint constraint, cfbasics::Real dPdx,
                         cfbasics::Real Ubulk, cfbasics::Real Uwall, cfbasics::Real h);

}  // namespace channelflow
#endif
