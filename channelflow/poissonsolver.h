/**
 * Solver for the Poisson equation
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#ifndef CHANNELFLOW_POISSONSOLVER_H
#define CHANNELFLOW_POISSONSOLVER_H

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"
#include "channelflow/diffops.h"
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/helmholtz.h"

namespace channelflow {

// Solve lapl u = f with dirichlet BCs. u is unknown, f is given.
class PoissonSolver {
   public:
    PoissonSolver();
    PoissonSolver(const FlowField& u);
    PoissonSolver(int Nx, int Ny, int Nz, int Nd, cfbasics::Real Lx, cfbasics::Real Lz, cfbasics::Real a, cfbasics::Real b, CfMPI* cfmpi);

    PoissonSolver(const PoissonSolver& ps);
    PoissonSolver& operator=(const PoissonSolver& ps);
    ~PoissonSolver();

    // solve lapl u = f with boundary conditions u=0 or u=bc at y=a and y=b
    void solve(FlowField& u, const FlowField& f) const;  // u=0 BC
    void solve(FlowField& u, const FlowField& f, const FlowField& bc) const;

    cfbasics::Real verify(const FlowField& u, const FlowField& f) const;
    cfbasics::Real verify(const FlowField& u, const FlowField& f, const FlowField& bc) const;

    bool geomCongruent(const FlowField& u) const;
    bool congruent(const FlowField& u) const;

   protected:
    int Mx_;  // number of X modes
    int My_;  // number of Chebyshev T(y) modes
    int Mz_;  // number of Z modes
    cfbasics::lint Mxloc_;
    cfbasics::lint Mzloc_;
    cfbasics::lint mxlocmin_;
    cfbasics::lint mzlocmin_;
    int Nd_;  // vector dimension
    cfbasics::Real Lx_;
    cfbasics::Real Lz_;
    cfbasics::Real a_;
    cfbasics::Real b_;

    HelmholtzSolver** helmholtz_;  // 2d cfarray of HelmHoltz solvers
};

// Solve lapl p = -div(nonl(u+U)) with BCs dp/dy = -nu d^2 (u+U)/dy^2 at y=a,b.
// The pressure p is unknown and u+U is a given velocity field.
// IMPORTANT NOTE: If nonl_method = Rotational, p is the modified pressure
// p + 1/2 u dot u. If you want the true pressure, set nonl_method to
// Convection, SkewSymmetric, Divergence, or Alternating.

class PressureSolver : public PoissonSolver {
   public:
    PressureSolver();
    PressureSolver(int Nx, int Ny, int Nz, cfbasics::Real Lx, cfbasics::Real Lz, cfbasics::Real a, cfbasics::Real b, const ChebyCoeff& U, const ChebyCoeff& W,
                   cfbasics::Real nu, cfbasics::Real Vsuck, NonlinearMethod nonl_method, CfMPI* cfmpi);

    PressureSolver(const FlowField& u, cfbasics::Real nu, cfbasics::Real Vsuck, NonlinearMethod nonl_method);

    PressureSolver(const FlowField& u, const ChebyCoeff& U, const ChebyCoeff& W, cfbasics::Real nu, cfbasics::Real Vsuck,
                   NonlinearMethod nonl_method);

    ~PressureSolver();

    // PressureSolver(const PoissonSolver& ps);
    // PressureSolver& operator=(const PoissonSolver& ps);

    // solve lapl p = -div nonl(u+U) with BCs dp/dy = d^2 (u+U)/dy^2 at y=a,b
    FlowField solve(const FlowField& u);
    void solve(FlowField& p, FlowField u);

    cfbasics::Real verify(const FlowField& p, const FlowField& u);

   private:
    ChebyCoeff U_;
    ChebyCoeff W_;
    ChebyTransform trans_;
    FlowField nonl_;
    FlowField tmp_;
    FlowField div_nonl_;
    cfbasics::Real nu_;
    cfbasics::Real Vsuck_;
    NonlinearMethod nonl_method_;
};

// Normally I avoid derived classes, because it's hard to apprehend the
// data structure and functionality of the derived class, particluarly when
// the two classes reside in different files. In this case the classes are
// simple enough and can go in the same file, since they naturally fall in
// the same place in the software dependency graph.
// jfg Wed May 18 18:12:13 EDT 2005

}  // namespace channelflow
#endif
