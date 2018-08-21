/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */

#ifndef NSOLVER_NEWTONALGORITHM_H
#define NSOLVER_NEWTONALGORITHM_H

#include <memory>
#include "cfbasics/arglist.h"
#include "cfbasics/cfbasics.h"

#include "cfbasics/brent.h"
#include "nsolver/bicgstabl.h"
#include "nsolver/fgmres.h"
#include "nsolver/gmres.h"
#include "nsolver/newton.h"

using namespace Eigen;

namespace nsolver {
//   /*==================================================================================*/
//   /*            function  GMRESHookstep and friends                                   */
//   /*==================================================================================*/
//
//   // Find solution of equation sigma f^T(u) - u == 0. Inputs are initial
//   // guess and return as solution values. Real return value is residual
//   // L2Norm(sigma f^T(u) - u)
//   Real GMRESHookstep(FlowField& u, Real& T, FieldSymmetry& sigma,
//                PoincareCondition* h, const GMRESHookstepFlags& searchflags,
//                const DNSFlags& dnsflags, const TimeStep& dt, Real& CFL);
//
//

//   // Versions of f, G, DG that handle Poincare section calculations, additionally.
//   void fp(const FlowField& u, Real& T, PoincareCondition* h,
//        FlowField& fu, const DNSFlags& flags, const TimeStep& dt, int& fcount, Real& CFL,
//        std::ostream& os=std::cout);
//
//   // sf(u,sigma,T) = sigma f^T(u)
//   void sfp(const FlowField& u, Real T, PoincareCondition* h,
//         const FieldSymmetry& sigma, FlowField& sfu, const DNSFlags& flags,
//         const TimeStep& dt, int& fcount, Real& CFL, std::ostream& os=std::cout);
//
//   // Dsf(u,sigma,T,du) = sigma f^T(u+du) - sigma f^T(u)
//   void Dsfp(const FlowField& u, Real T, PoincareCondition* h,
//          const FieldSymmetry& sigma, const FlowField& sfu, const FlowField& du,
//          FlowField& Dsf_du, const DNSFlags& flags, const TimeStep& dt, Real eps,
//          bool centerdiff, int& fcount, Real& CFL, std::ostream& os=std::cout);
//
//   void Gp(const FlowField& u, Real& T, PoincareCondition* h,
//        const FieldSymmetry& sigma, FlowField& Gu, const DNSFlags& flags, const TimeStep& dt,
//        bool Tnormalize, bool Unormalize, int& fcount, Real& CFL, std::ostream& os=std::cout) ;
//
//   void DGp(const FlowField& u, const FlowField& du, Real& T, Real dT,
//         PoincareCondition* h, const FieldSymmetry& sigma, const FieldSymmetry& dsigma,
//         const FlowField& Gu, FlowField& DG_dx, const DNSFlags& flags, const TimeStep& dt,
//         bool Tnormalize, bool Unormalize, Real epsDu, bool centdiff, int& fcount, Real& CFL, std::ostream&
//         os=std::cout);

class NewtonSearchFlags {
   public:
    SolverMethod solver = SolverGMRES;
    OptimizationMethod optimization = Hookstep;
    SolutionType solntype;
    bool xrelative;  // new
    bool zrelative;  // new
    Real epsSearch;
    Real epsKrylov;
    Real epsDx;
    Real epsDt;
    Real epsSolver;
    Real epsSolverF;
    bool centdiff;
    int Nnewton;
    int Nsolver;
    int Nhook;
    Real delta;
    Real deltaMin;
    Real deltaMax;
    Real deltaFuzz;
    Real lambdaMin;
    Real lambdaMax;
    Real lambdaRequiredReduction;
    Real improvReq;
    Real improvOk;
    Real improvGood;
    Real improvAcc;
    int lBiCGStab;
    int nShot;
    bool fixtphase;
    Real TRef;
    Real axRef;
    Real azRef;
    Real gRatio;
    std::string outdir;
    std::ostream* logstream;
    bool verbose;
    bool orbit;
    bool laurette;

    NewtonSearchFlags(SolutionType solntype = Equilibrium, bool xrelative = false, bool zrelative = false,
                      Real epsSearch = 1e-13, Real epsKrylov = 1e-14, Real epsDx = 1e-7, Real epsDt = 1e-5,
                      Real epsSolver = 1e-3, Real epsSolverF = 0.05, bool centdiff = false, int Nnewton = 20,
                      int Ngmres = 500, int Nhook = 20, Real delta = 1e-2, Real deltaMin = 1e-12, Real deltaMax = 1e-1,
                      Real deltaFuzz = 1e-6, Real lambdaMin = 0.2, Real lambdaMax = 1.5,
                      Real lambdaRequiredReduction = 0.5, Real improvReq = 1e-3, Real improvOk = 1e-1,
                      Real improvGood = 0.75, Real improvAcc = 1e-1, int lBiCGStab = 2, int nShot = 1,
                      bool fixtphase = false, Real TRef = 1.0, Real axRef = 1.0, Real azRef = 1.0, Real gRatio = 10.0,
                      std::string outdir = "./", std::ostream* logstream = &std::cout, bool laurette = false);

    NewtonSearchFlags(ArgList& args);
    SolverMethod string2solver(std::string);
    string solver2string() const;
    OptimizationMethod string2optimization(std::string);
    string optimization2string() const;
    SolutionType string2solntype(string s) const;
    string solntype2string() const;
    void save(const string& outdir = "") const;
    void load(int taskid, const string indir);
    const vector<string> getFlagList();
};

/** Combined Newton Hookstep algorithm
 *
 * \param[in] dsiG DynamicalSystemsInterface specifying the equation to integrate
 * \param[in] x0 initial guess
 * \param[in] searchflags parameters for the search algorithms
 *
 * \return fixed point/one point on periodic orbit
 */
VectorXd hookstepSearch(DSI& dsiG, const VectorXd& x0, const NewtonSearchFlags& searchflags, Real& gx);

class NewtonAlgorithm : public Newton {
   public:
    NewtonAlgorithm(NewtonSearchFlags searchflags);
    virtual VectorXd solve(DSI& dsi, const VectorXd& x, Real& residual);
    NewtonSearchFlags searchflags;

    virtual void setLogstream(std::ostream* os);
    virtual void setOutdir(std::string od);

    Eigen::MatrixXd jacobi(const Eigen::VectorXd& x, const Real epsilon, const bool centerdiff, int& fcount);

   private:
    int linear(VectorXd& dxOpt, VectorXd& GxOpt, const VectorXd& x);
    int hookstep(VectorXd& dxH, VectorXd& GxH, const VectorXd& x, const VectorXd& b);
    int convergenceCheckAC(const VectorXd& dx, VectorXd& Gx, const VectorXd& x);

    ostream* os;
    int fcount_newton_;
    int fcount_opt_;

    // required for Hookstep algorithm
    unique_ptr<GMRES> gmres_;
    unique_ptr<FGMRES> fgmres_;
    Real delta_;
    Real rx_;  // Dennis & Schnabel residual r(x) = 1/2 ||f(x)||^2
};

}  // namespace nsolver

#endif
