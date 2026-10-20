/**
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <iostream>
#include "cfbasics/cfbasics.h"
#include "channelflow/diffops.h"
#include "channelflow/flowfield.h"
#include "modules/settling/settlingdsi.h"

using namespace std;
using namespace Eigen;
using namespace chflow;

// Largest violation, over the locally owned Fourier modes, of the constraints that
// SettlingDSI::vectorToField must impose on a 4-component state (velocity + settling scalar):
//   no-slip velocity u = v = w = 0 at the walls,
//   no-flux scalar rho' + alpha rho = g at the walls (g = 1 -/+ alpha for the (0,0) mode, else 0),
//   mean(rho) = 1 (unit mass loading).
// Divergence is checked separately with divNorm.
Real constraintError(const FlowField& u, Real alpha) {
    Real err = 0;
    for (int mx = u.mxlocmin(); mx < u.mxlocmin() + u.Mxloc(); ++mx)
        for (int mz = u.mzlocmin(); mz < u.mzlocmin() + u.Mzloc(); ++mz) {
            const bool mean = (u.kx(mx) == 0 && u.kz(mz) == 0);
            for (int i = 0; i < 3; ++i) {
                ComplexChebyCoeff f = u.profile(mx, mz, i);
                err = Greater(err, Greater(abs(f.eval_a()), abs(f.eval_b())));
            }
            ComplexChebyCoeff c = u.profile(mx, mz, 3);
            const Complex ga = mean ? Complex(1.0 - alpha, 0.0) : Complex(0.0, 0.0);
            const Complex gb = mean ? Complex(1.0 + alpha, 0.0) : Complex(0.0, 0.0);
            err = Greater(err, abs(c.slope_a() + alpha * c.eval_a() - ga));
            err = Greater(err, abs(c.slope_b() + alpha * c.eval_b() - gb));
            if (mean)
                err = Greater(err, abs(c.re.mean() - 1.0));
        }
#ifdef HAVE_MPI
    Real globalerr = 0;
    MPI_Allreduce(&err, &globalerr, 1, MPI_DOUBLE, MPI_MAX, *u.comm_world());
    err = globalerr;
#endif
    return err;
}

// Test that SettlingDSI's fieldToVector and vectorToField are inverse to each other on states
// satisfying the constraints, and that vectorToField imposes the constraints on an arbitrary
// vector. On 2 MPI processes this also covers the conjugate (-kx, 0) modes owned by another process.
int main(int argc, char* argv[]) {
    int failure = 0;
    cfMPI_Init(&argc, &argv);
    {
        CfMPI* cfmpi = &CfMPI::getInstance();
        // eq3 regression reference: travelling wave at Ri = 0.0025, vs/kappa = 1
        FlowField u("data/settle_eq3.nc", cfmpi);
        const Real alpha = 1.0;  // vs/kappa
        DNSFlags flags;
        flags.nu = 1.0 / 400;
        flags.Ri = 0.0025;
        flags.vs = 2.5e-3;
        flags.kappa = 2.5e-3;
        TimeStep dt(0.02, 0.02, 0.02, 1.0, 0.4, 0.6, false);
        SettlingDSI dsi(flags, FieldSymmetry(), 0, dt, false, false, false, false, 0.0, u, &cout);
        const Real unorm = L2Norm(u);
        // The scalar BC involves the wall slope sum(n^2 a_n), which amplifies rounding by up to
        // (Ny-1)^2 ~ 1e3: the stored solution itself meets its constraints only to ~3e-13. Real
        // packing errors show up at O(1e-2).
        const Real tol = 1e-11;

        // 1. Round trip from a field that satisfies the constraints
        cout << "constraint error of input field   = " << constraintError(u, alpha) << endl;
        VectorXd v, v2;
        dsi.fieldToVector(u, v);
        if (v.size() != dsi.fieldVectorSize(u)) {
            cout << "fieldToVector size " << v.size() << " != fieldVectorSize " << dsi.fieldVectorSize(u) << endl;
            failure = 1;
        }
        FlowField u2(u);
        dsi.vectorToField(v, u2);
        dsi.fieldToVector(u2, v2);
        Real dvector2 = (v2 - v).squaredNorm(), vnorm2 = v.squaredNorm();
#ifdef HAVE_MPI
        Real buf[2] = {dvector2, vnorm2}, sum[2];
        MPI_Allreduce(buf, sum, 2, MPI_DOUBLE, MPI_SUM, *u.comm_world());
        dvector2 = sum[0];
        vnorm2 = sum[1];
#endif
        const Real dfield = L2Dist(u, u2) / unorm;
        const Real dvector = sqrt(dvector2 / vnorm2);
        cout << "round trip: rel L2Dist(u, u2)     = " << dfield << endl;
        cout << "round trip: rel |v2 - v|          = " << dvector << endl;
        if (dfield > tol || dvector > tol)
            failure = 1;

        // 2. An arbitrary vector must give a field satisfying every constraint
        srand(1 + u.taskid());
        VectorXd r = VectorXd::Random(v.size()) * 0.1;
        FlowField ur(u);
        dsi.vectorToField(r, ur);
        FlowField uvel, c;
        splitVelocityScalar(ur, uvel, c);
        const Real cerr = constraintError(ur, alpha);
        const Real derr = divNorm(uvel) / L2Norm(uvel);
        cout << "random vector: constraint error   = " << cerr << endl;
        cout << "random vector: rel divNorm        = " << derr << endl;
        if (cerr > tol || derr > tol)
            failure = 1;

        // 3. fixdivnoslipSettling: an interpolated or perturbed state is projected onto the constraints
        FlowField up(u);
        up += ur;
        up *= 0.5;
        if (up.taskid() == up.task_coeff(0, 0))
            up.cmplx(0, 0, 0, 3) += Complex(0.01, 0.0);  // breaks the BC and the mass
        fixdivnoslipSettling(up, settlingScalars(flags)[0]);
        const Real perr = constraintError(up, alpha);
        cout << "fixdivnoslipSettling: constraint error = " << perr << endl;
        if (perr > tol)
            failure = 1;

        cout << (failure ? "\t** FAIL **" : "\t   pass   ") << endl;
    }
    cfMPI_Finalize();
    return failure;
}
