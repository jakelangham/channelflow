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

// As constraintError, for a stratified scalar: no slip, and c = 0 (Dirichlet) or c' = 0 (Neumann) at
// the walls in every mode, and with pinMean mean(c) = 0
Real stratConstraintError(const FlowField& u, const SettlingScalar& sc) {
    Real err = 0;
    for (int mx = u.mxlocmin(); mx < u.mxlocmin() + u.Mxloc(); ++mx)
        for (int mz = u.mzlocmin(); mz < u.mzlocmin() + u.Mzloc(); ++mz) {
            for (int i = 0; i < 3; ++i) {
                ComplexChebyCoeff f = u.profile(mx, mz, i);
                err = Greater(err, Greater(abs(f.eval_a()), abs(f.eval_b())));
            }
            ComplexChebyCoeff c = u.profile(mx, mz, 3);
            if (sc.wallbc == NeumannWalls)
                err = Greater(err, Greater(abs(c.slope_a()), abs(c.slope_b())));
            else
                err = Greater(err, Greater(abs(c.eval_a()), abs(c.eval_b())));
            if (sc.pinMean() && u.kx(mx) == 0 && u.kz(mz) == 0)
                err = Greater(err, abs(c.re.mean()));
        }
#ifdef HAVE_MPI
    Real globalerr = 0;
    MPI_Allreduce(&err, &globalerr, 1, MPI_DOUBLE, MPI_MAX, *u.comm_world());
    err = globalerr;
#endif
    return err;
}

// Packing of a stratified scalar (gravity along y or z, Dirichlet or Neumann walls): vectorToField
// imposes the constraints on an arbitrary vector, and fieldToVector inverts it; and a symmetry acts on
// the 4-component state with the scalar's parity (as u_g). Returns 1 on failure.
int stratifiedCase(const FlowField& u0, ScalarModel model, ScalarWallBC wallbc, Real tol) {
    DNSFlags flags;
    flags.nu = 1.0 / 400;
    flags.Ri = 0.01;
    flags.Sc = 7.0;
    flags.scalarmodel = model;
    flags.scalarbc = wallbc;
    const SettlingScalar sc = settlingScalars(flags)[0];
    TimeStep dt(0.02, 0.02, 0.02, 1.0, 0.4, 0.6, false);
    SettlingDSI dsi(flags, FieldSymmetry(), 0, dt, false, false, false, false, 0.0, u0, &cout);
    int failure = 0;

    VectorXd v, v2;
    dsi.fieldToVector(u0, v);
    srand(3 + u0.taskid());
    VectorXd r = VectorXd::Random(v.size()) * 0.1;
    FlowField u(u0), u2(u0);
    dsi.vectorToField(r, u);
    const Real cerr = stratConstraintError(u, sc);
    dsi.fieldToVector(u, v);
    dsi.vectorToField(v, u2);
    dsi.fieldToVector(u2, v2);
    Real d2 = (v2 - v).squaredNorm(), n2 = v.squaredNorm();
#ifdef HAVE_MPI
    Real buf[2] = {d2, n2}, sum[2];
    MPI_Allreduce(buf, sum, 2, MPI_DOUBLE, MPI_SUM, *u.comm_world());
    d2 = sum[0];
    n2 = sum[1];
#endif
    const Real dfield = L2Dist(u, u2) / L2Norm(u);
    const Real dvector = sqrt(d2 / n2);

    // the scalar changes sign with the reflection of its gravity axis: s1 = [u, v, -w](x + Lx/2, y, -z)
    // flips it along z, the shift-rotation [-u, -v, w](-x + Lx/2, -y, z + Lz/2) along y
    const int g = sc.gravity();
    const FieldSymmetry sigma = (g == 2) ? FieldSymmetry(1, 1, -1, 0.5, 0.0) : FieldSymmetry(-1, -1, 1, 0.5, 0.5);
    FlowField su4 = sigma(u), uv, cv, suv, scv;
    splitVelocityScalar(u, uv, cv);
    splitVelocityScalar(su4, suv, scv);
    FlowField mscv = sigma(cv);
    mscv *= -1.0;
    const Real svel = L2Dist(suv, sigma(uv)) / L2Norm(uv);
    const Real sscal = L2Dist(scv, mscv) / L2Norm(cv);
    // and the DNS's symmetry projection list for the scalar matches
    cfarray<FieldSymmetry> one(1);
    one[0] = sigma;
    FlowField pcv = scalarSymmetries(one, sc.parityAxis())[0](cv);
    const Real sproj = L2Dist(pcv, mscv) / L2Norm(cv);

    cout << scalarmodel2string(model) << " " << scalarwallbc2string(wallbc) << ": constraint error " << cerr
         << ", round trip " << dfield << " " << dvector << ", symmetry velocity " << svel << " scalar " << sscal
         << " projection " << sproj << endl;
    if (cerr > tol || dfield > tol || dvector > tol || svel > tol || sscal > tol || sproj > tol)
        failure = 1;
    return failure;
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

        // 4. A symmetry acts on the velocity as a vector and on the scalar as a scalar (overall sign
        // only), e.g. x-reflection with a shift: (u, v, w, c)(x, y, z) -> (-u, v, w, c)(-x + ax, y, z + az)
        const FieldSymmetry sigma(-1, 1, 1, 0.25, 0.5);
        FlowField su4 = sigma(u), uv, cv, suv, scv;
        splitVelocityScalar(u, uv, cv);
        splitVelocityScalar(su4, suv, scv);
        const Real svel = L2Dist(suv, sigma(uv)) / L2Norm(uv);
        const Real sscal = L2Dist(scv, sigma(cv)) / L2Norm(cv);
        cout << "symmetry: velocity " << svel << ", scalar " << sscal << endl;
        if (svel > tol || sscal > tol)
            failure = 1;

        // 5. Stratified scalars, on the same grid (their SettlingDSI sets the scalar's parity axis)
        for (ScalarModel model : {StratifiedY, StratifiedZ})
            for (ScalarWallBC wallbc : {DirichletWalls, NeumannWalls})
                failure |= stratifiedCase(u, model, wallbc, tol);
        setScalarParityAxis(-1);

        // 6. dnsflags.txt: the stratified scalar's flags survive a save/load, and a file from before
        // them (%Pr in place of %Sc, no scalar model or BC) loads as a settling scalar
        {
            DNSFlags fs;
            fs.Sc = 7.0;
            fs.Ri = 1.5;
            fs.scalarmodel = StratifiedZ;
            fs.scalarbc = NeumannWalls;
            mkdir("flagtest");
            fs.save("flagtest/");
            DNSFlags fl, fo;
            fl.load(u.taskid(), "flagtest/");
            const bool okNew = fl.Sc == 7.0 && fl.Ri == 1.5 && fl.scalarmodel == StratifiedZ &&
                               fl.scalarbc == NeumannWalls && scalarParityAxis() == 2;
            fo.scalarmodel = StratifiedY;  // must be reset by the legacy load
            fo.load(u.taskid(), "data/legacy/");
            const bool okOld = fo.scalarmodel == SettlingModel && fo.kappa == 0.0025 && fo.vs == 2.5e-5 &&
                               fo.Sc == 1.0 && scalarParityAxis() == -1;
            cout << "dnsflags.txt: new format " << (okNew ? "ok" : "WRONG") << ", legacy format "
                 << (okOld ? "ok" : "WRONG") << endl;
            if (!okNew || !okOld)
                failure = 1;
        }

        cout << (failure ? "\t** FAIL **" : "\t   pass   ") << endl;
    }
    cfMPI_Finalize();
    return failure;
}
