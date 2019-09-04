/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: Tobias Kreilos
 */

#include <channelflow/laurettedsi.h>
#include "channelflow/cfdsi.h"
#include "channelflow/flowfield.h"
#include "nsolver/nsolver.h"

using namespace std;
using namespace Eigen;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        ArgList args(argc, argv, "find an invariant solution using Newton-Krylov-hookstep algorithm");

        /** Choose the Newton algorithm to be used. Currently, two options are available: simple Newton without any
         * trust region optimization, and Newton with Hookstep (default). For the simple Newton, you can choose either a
         * full-space algorithm to solve the Newton equations (-solver "eigen") or between the two iterative algorithms
         * GMRES and BiCGStab. Newton-Hookstep requires GMRES. Note that the available parameters depend on your choice
         * of the algorithm.
         */

        unique_ptr<Newton> N;
        NewtonSearchFlags searchflags(args);
        searchflags.save(searchflags.outdir);
        N = unique_ptr<Newton>(new NewtonAlgorithm(searchflags));

        DNSFlags dnsflags(args, searchflags.laurette);
        TimeStep dt(dnsflags);

        bool Rxsearch, Rzsearch, Tsearch;
        Rxsearch = searchflags.xrelative;
        Rzsearch = searchflags.zrelative;
        Tsearch = searchflags.solntype == PeriodicOrbit ? true : false;

        const bool Tnormalize = (Tsearch || searchflags.laurette) ? false : true;

        /** Read in remaining arguments */

        args.section("Program options");
        const string sigmastr =
            args.getstr("-sigma", "--sigma", "", "file containing sigma of sigma f^T(u) - u = 0 (default == identity)");
        const Real unormalize = args.getreal("-un", "--unormalize", 0.0, "lower bound in energy for search");
        const int nproc0 =
            args.getint("-np0", "--nproc0", 0, "number of MPI-processes for transpose/number of parallel ffts");
        const int nproc1 = args.getint("-np1", "--nproc1", 0, "number of MPI-processes for one fft");
        const bool msinit =
            args.getflag("-MSinit", "--MSinitials", "read different files as the initial guesses for different shoots");
        const string rhofile = args.getstr(1, "<flowfield>", "initial guess for the solution");
        const string velfile = args.getstr(2, "<flowfield>", "precomputed velocity field");

        args.check();
        args.save();
        WriteProcessInfo(argc, argv);
        dnsflags.save();
        cout << dnsflags << endl;

        CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);

        FlowField rhoin(rhofile, cfmpi);
        Real vs_o_kappa = dnsflags.vs / dnsflags.kappa;
        BoundaryCond bc(Mixed, rhoin.Mx(), rhoin.Mz(), 1.0 + vs_o_kappa, vs_o_kappa);
        FlowField rho(rhoin.Nx(), rhoin.Ny(), rhoin.Nz(), 1, 
                      rhoin.Lx(), rhoin.Lz(), rhoin.a(), rhoin.b(), bc, cfmpi);
        // JL if rhofile has 4 dimensions assume 4th is the initial density,
        // else assume the 1st
        if (rhoin.Nd() == 4) {
            rho.copySubfields(rhoin, {3}, {0});
        } else if (rhoin.Nd() == 1) {
            rho = rhoin;
        } else {
            cout << "rho file must have 1 or 4 dimensions\n";
            cfMPI_Finalize();
            return 1;
        }

        // upscale u if necessary for the later dot product in physical space
        FlowField uin(velfile, cfmpi);
        FlowField u(rho.Nx(), rho.Ny(), rho.Nz(), uin.Nd(), uin.Lx(), uin.Lz(), uin.a(), uin.b(), uin.BC(), rho.cfmpi());
        u.interpolate(uin);

        rho.set_nonconstant_ga(u);

        FieldSymmetry sigma;
        if (sigmastr.length() != 0)
            sigma = FieldSymmetry(sigmastr);

        /** Construct the dynamical-systems interface object depending on the given parameters. Current options are
         * either standard (f(u) via forward time integration) or Laurette (f(u) via Laurettes method)
         */
        unique_ptr<cfDSI> dsi;
        dsi = unique_ptr<cfDSI>(new cfDSI(dnsflags, sigma, 0, dt, Tsearch, 
            Rxsearch, Rzsearch, Tnormalize, unormalize,
            u, rho, N->getLogstream()));

        VectorXd x_singleShot;
        VectorXd x;
        VectorXd yvec;
        MatrixXd y;
        MultishootingDSI* msDSI = N->getMultishootingDSI();
        dsi->makeVector(rho, sigma, dnsflags.T, x_singleShot);
        msDSI->setDSI(*dsi, x_singleShot.size());
        if (msinit) {
            int nSh = msDSI->nShot();
            y.resize(x_singleShot.size(), nSh);
            Real Tms = dnsflags.T / nSh;
            vector<FlowField> u_ms(nSh);
            u_ms[0] = rho;
            for (int i = 1; i < nSh; i++) {
                string uname_ms = "./Multishooting/" + rhofile + i2s(i);
                FlowField ui(uname_ms, cfmpi);
                u_ms[i] = ui;
            }
            for (int i = 0; i < nSh; i++) {
                dsi->makeVector(u_ms[i], sigma, Tms, yvec);
                y.col(i) = yvec;
            }
            x = msDSI->toVector(y);
        } else {
            x = msDSI->makeMSVector(x_singleShot);
        }

        int Nunk = x.size();
        int Nunk_total = Nunk;
#ifdef HAVE_MPI
        MPI_Allreduce(&Nunk, &Nunk_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
        cout << Nunk_total << " unknowns" << endl;

        Real residual = 0;
        N->solve(*dsi, x, residual);
    }

    cfMPI_Finalize();

    return 0;
}
