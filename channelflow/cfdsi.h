/**
 * Channelflow Dynamical System Interface (DSI)
 *
 * Interface to use NSolver's dynamical system methods
 * with the Channelflow Navier-Stokes solver.
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#ifndef CFDSI_H
#define CFDSI_H

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "cfbasics/cfvector.h"
#include "channelflow/cfmpi.h"
#include "channelflow/chebyshev.h"
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/periodicfunc.h"
#include "channelflow/symmetry.h"
#include "channelflow/tausolver.h"
#include "channelflow/utilfuncs.h"
#include "nsolver/nsolver.h"

using namespace std;

namespace channelflow {

enum class continuationParameter {
    T,
    Re,
    P,
    Ub,
    Uw,
    ReP,
    Theta,
    ThArc,
    ThLx,
    ThLz,
    Lx,
    Lz,
    Aspect,
    Diag,
    Lt,
    Vs,
    ReVs,
    H,
    HVs,
    Rot
};

Real GMRESHookstep_vector(FlowField& u, FieldSymmetry& sigma, PoincareCondition* h,
                          const nsolver::NewtonSearchFlags& searchflags, DNSFlags& dnsflags, TimeStep& dt, Real& CFL,
                          Real Unormalize);

// Real meandissipation (const FlowField& uarg, Real T, DNSFlags dnsflags,
//                       const TimeStep& dtarg, SolutionType solntype);

// converts the string from "fieldstats" in diffops to a vector of Reals
vector<Real> fieldstats_vector(const FlowField& u);

class cfDSI : public nsolver::DSI {
   public:
    /** \brief default constructor */
    cfDSI();
    virtual ~cfDSI() {}

    /** \brief Initialize cfDSI */
    cfDSI(DNSFlags& dnsflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch, bool xrelative,
          bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u, std::ostream* os = &std::cout);

    VectorXd eval(const VectorXd& x) override;
    VectorXd eval(const VectorXd& x0, const VectorXd& x1, bool symopt) override;
    void save(const VectorXd& x, const string filebase, const string outdir = "./",
              const bool fieldsonly = false) override;
    void saveEigenvec(const VectorXd& x, const string label, const string outdir) override;
    void saveEigenvec(const VectorXd& x1, const VectorXd& x2, const string label1, const string label2,
                      const string outdir) override;

    Real DSIL2Norm(const VectorXd& x) override;
    string stats(const VectorXd& x) override;
    pair<string, string> stats_minmax(const VectorXd& x) override;
    string statsHeader() override;
    void makeVector(const channelflow::FlowField& u, const FieldSymmetry& sigma, const Real T, VectorXd& x);
    void extractVector(const VectorXd& x, FlowField& u, FieldSymmetry& sigma, Real& T);
    void toVector(const vector<FlowField>& u, const FieldSymmetry& sigma, const Real T, VectorXd& x){};

    /// \name Compute derivatives of FlowField corresponding to this vector
    VectorXd xdiff(const VectorXd& a) override;
    VectorXd zdiff(const VectorXd& a) override;
    VectorXd tdiff(const VectorXd& a, Real epsDt) override;

    /// \name Handle continuation parameter
    void updateMu(Real mu) override;
    void chooseMu(std::string muName);
    void chooseMu(continuationParameter mu);
    string printMu() override;  // document
    void saveParameters(string searchdir) override;
    continuationParameter s2cPar(std::string muName);
    string cPar2s(continuationParameter cPar);
    void phaseShift(VectorXd& x) override;
    void phaseShift(MatrixXd& y) override;
    inline void setPhaseShifts(bool xphasehack, bool zphasehack, bool uUbasehack);
    Real observable(VectorXd& x) override;

    Real tph_observable(VectorXd& x) override;
    Real extractT(const VectorXd& x) override;
    Real extractXshift(const VectorXd& x) override;
    Real extractZshift(const VectorXd& x) override;

    Real getCFL() const { return CFL_; };
    bool XrelSearch() const override { return xrelative_; };
    bool ZrelSearch() const override { return zrelative_; };
    bool Tsearch() const override { return Tsearch_; };

   protected:
    DNSFlags dnsflags_;
    CfMPI* cfmpi_;
    FieldSymmetry sigma_;
    PoincareCondition* h_;
    TimeStep dt_;
    bool Tsearch_;
    bool xrelative_;
    bool zrelative_;
    Real Tinit_;
    Real axinit_;
    Real azinit_;
    bool Tnormalize_;
    Real Unormalize_;
    int fcount_;
    int Nx_;
    int Ny_;
    int Nz_;
    int Nd_;
    Real Lx_;
    Real Lz_;
    Real ya_;
    Real yb_;
    Real CFL_;
    int uunk_;
    bool xphasehack_;
    bool zphasehack_;
    bool uUbasehack_;

    continuationParameter cPar_ = continuationParameter::T;
};

inline void cfDSI::setPhaseShifts(bool xphasehack, bool zphasehack, bool uUbasehack) {
    xphasehack_ = xphasehack;
    zphasehack_ = zphasehack;
    uUbasehack_ = uUbasehack;
}

void f(const FlowField& u, int N, Real dt, FlowField& f_u, const DNSFlags& flags_, ostream& os);

// Versions of f, G, DG that handle Poincare section calculations, additionally.
void f(const FlowField& u, Real& T, PoincareCondition* h, FlowField& fu, const DNSFlags& flags, const TimeStep& dt,
       int& fcount, Real& CFL, std::ostream& os = std::cout);

void G(const FlowField& u, Real& T, PoincareCondition* h, const FieldSymmetry& sigma, FlowField& Gu,
       const DNSFlags& flags, const TimeStep& dt, bool Tnormalize, Real Unormalize, int& fcount, Real& CFL,
       std::ostream& os = std::cout);

}  // namespace channelflow

#endif
