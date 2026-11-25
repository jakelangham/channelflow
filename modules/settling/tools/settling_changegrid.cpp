/**
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <fstream>
#include <iomanip>
#include <iostream>
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/utilfuncs.h"
#include "modules/settling/settlingdsi.h"

using namespace std;
using namespace Eigen;

using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose(
            "interpolate a given flowfield onto a different grid. A 4-component field (velocity and "
            "settling scalar in component 3) keeps its scalar; with -vsok, -dv also imposes the "
            "scalar's no-flux BC and unit mass, with -strat y|z a stratified scalar's BC (-sbc)");
        ArgList args(argc, argv, purpose);  // Assign parameters

        const bool padded = args.getbool("-p", "--padded", true, "set padding modes to zero");
        const bool fixdiv = args.getbool("-dv", "--divergence", true, "fix divergence and Dirichlet BCs");
        const int Nxarg = args.getint("-Nx", "--Nx", 0, "new # x gridpoints [default==no change]");
        const int Nyarg = args.getint("-Ny", "--Ny", 0, "new # y gridpoints [default==no change]");
        const int Nzarg = args.getint("-Nz", "--Nz", 0, "new # z gridpoints [default==no change]");
        const Real aarg = args.getreal("-a", "--a", NAN, "new lower wall position [default==no change]");
        const Real barg = args.getreal("-b", "--b", NAN, "new upper wall position [default==no change]");
        const Real alpha = args.getreal("-al", "--alpha", 0.0, "new alpha == 2pi/Lx");
        const Real gamma = args.getreal("-ga", "--gamma", 0.0, "new gamma == 2pi/Lz");
        const Real vs_o_kappa =
            args.getreal("-vsok", "--vsoverkappa", 0.0, "vs/kappa of the settling scalar (for -dv; 0 == leave it)");
        const ScalarModel model = s2scalarmodel(args.getstr(
            "-strat", "--stratification", "none", "stratified scalar along y or z (for -dv; none == see -vsok)"));
        const ScalarWallBC wallbc =
            s2scalarwallbc(args.getstr("-sbc", "--scalarbc", "diri", "a stratified scalar's wall BC, [diri|neum]"));
        if (model != SettlingModel && vs_o_kappa != 0)
            cferror("settling_changegrid : give -vsok (settling scalar) or -strat (stratified scalar), not both");
        const Real lx = (alpha == 0) ? args.getreal("-lx", "--lx", 0.0, "new Lx = 2 pi lx") : 1.0 / alpha;
        const Real lz = (gamma == 0) ? args.getreal("-lz", "--lz", 0.0, "new Lz = 2 pi lz") : 1.0 / gamma;
        /***/ Real Lx = (lx == 0.0) ? args.getreal("-Lx", "--Lx", 0.0, "streamwise (x) box length") : 2 * pi * lx;
        /***/ Real Lz = (lz == 0.0) ? args.getreal("-Lz", "--Lz", 0.0, "spanwise   (z) box length") : 2 * pi * lz;
        const int nproc0 =
            args.getint("-np0", "--nproc0", 0, "number of MPI-processes for transpose/number of parallel ffts");
        const int nproc1 = args.getint("-np1", "--nproc1", 0, "number of MPI-processes for one fft");

        const string u0name = args.getstr(2, "<infield>", "input flow field");
        const string u1name = args.getstr(1, "<outfield>", "output flow field");

        args.check();
        args.save("./");

        CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);
        FlowField u0(u0name, cfmpi);

        const int Nx = (Nxarg == 0 ? u0.Nx() : Nxarg);
        const int Ny = (Nyarg == 0 ? u0.Ny() : Nyarg);
        const int Nz = (Nzarg == 0 ? u0.Nz() : Nzarg);
        const Real a = (std::isnan(aarg) ? u0.a() : aarg);
        const Real b = (std::isnan(barg) ? u0.b() : barg);

        FlowField u1(Nx, Ny, Nz, u0.Nd(), u0.Lx(), u0.Lz(), u0.a(), u0.b(), cfmpi);
        u1.interpolate(u0);

        if (Lx == 0)
            Lx = u0.Lx();
        if (Lz == 0)
            Lz = u0.Lz();

        if (Lx != u0.Lx() || Lz != u0.Lz())
            u1.rescale(Lx, Lz);

        cout << setprecision(16);
        cout << "L2Norm(u0)  == " << L2Norm(u0) << endl;
        cout << "L2Norm(u1)  == " << L2Norm(u1) << endl;
        cout << "bcNorm(u0)  == " << bcNorm(u0) << endl;
        cout << "bcNorm(u1)  == " << bcNorm(u1) << endl;

        if (u0.Nd() == 3 || u0.Nd() == 4) {
            // velocity (and, for 4 components, scalar) of u0 and u1
            FlowField v0, v1, c0, c1;
            if (u0.Nd() == 4) {
                splitVelocityScalar(u0, v0, c0);
                splitVelocityScalar(u1, v1, c1);
            } else {
                v0 = u0;
                v1 = u1;
            }
            cout << "divNorm(u0) == " << divNorm(v0) << endl;
            cout << "divNorm(u1) == " << divNorm(v1) << endl;

            // if ((Ny < u0.Ny() || a != u0.a() || b != u0.b()) & fixdiv) {
            if (fixdiv) {
                VectorXd v;
                FlowField foo(Nx, Ny, Nz, 3, Lx, Lz, a, b, cfmpi);
                field2vector(v1, v);
                vector2field(v, foo);
                // cout << "L2Dist(u1,2)== " << L2Dist(u1,foo) << endl;   //Commented by Sajjad: becuase this program
                // hit an assertion failure when it wants to rescale field in y direction!!!
                v1 = foo;
                cout << "L2Norm(u2)  == " << L2Norm(v1) << endl;
                cout << "divNorm(u2) == " << divNorm(v1) << endl;
                cout << "bcNorm(u2)  == " << bcNorm(v1) << endl;
                if (u0.Nd() == 4) {
                    // the scalar's coefficients on the walls a, b, as for the velocity above
                    FlowField c2(Nx, Ny, Nz, 1, Lx, Lz, a, b, cfmpi);
                    c1.makeSpectral();
                    for (int mx = c1.mxlocmin(); mx < c1.mxlocmin() + c1.Mxloc(); ++mx)
                        for (int mz = c1.mzlocmin(); mz < c1.mzlocmin() + c1.Mzloc(); ++mz)
                            for (int ny = 0; ny < Ny; ++ny)
                                c2.cmplx(mx, ny, mz, 0) = c1.cmplx(mx, ny, mz, 0);
                    c1 = c2;
                    if (vs_o_kappa != 0) {
                        projectScalar(c1, SettlingScalar{0.0, vs_o_kappa, 1.0});
                        cout << "scalar: imposed the no-flux BC (vs/kappa == " << vs_o_kappa << ") and unit mass"
                             << endl;
                    } else if (model != SettlingModel) {
                        const SettlingScalar sc(0.0, 0.0, 1.0, model, wallbc);
                        projectScalar(c1, sc);
                        cout << "scalar: imposed the stratified scalar's " << scalarwallbc2string(wallbc) << " BC"
                             << (sc.pinMean() ? " and zero mean" : "") << endl;
                    } else {
                        cout << "scalar: BC and mass unchanged (give -vsok or -strat to impose them)" << endl;
                    }
                }
            } else {
                cout << "Not fixing divergence..." << endl;
            }
            if (u0.Nd() == 4)
                combineVelocityScalar(v1, c1, u1);
            else
                u1 = v1;
        }

        if (padded || (Nx >= (3 * u0.Nx()) / 2 && Nx >= (3 * u0.Nz()) / 2))
            u1.zeroPaddedModes();

        u1.save(u1name);
    }
    cfMPI_Finalize();
}
