/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/symmetry.h"
#include "channelflow/tausolver.h"
#include "channelflow/utilfuncs.h"

using namespace std;
using namespace chflow;

string printdiagnostics(FlowField& rho, FlowField& u, const DNS& dns, Real t, const TimeStep& dt, Real nu, Real umin, bool vardt,
                        bool pl2norm, bool pchnorm, bool pdissip, bool pshear, bool pcfl);

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        WriteProcessInfo(argc, argv);
        string purpose(
            "integrate plane Couette or channel flow from a given "
            "initial condition and save velocity fields to disk.");

        ArgList args(argc, argv, purpose);

        DNSFlags flags(args);
        TimeStep dt(flags);

        args.section("Program options");
        const string outdir = args.getpath("-o", "--outdir", "data/", "output directory");
        const string label = args.getstr("-l", "--label", "u", "output field prefix");

        const bool pcfl = args.getflag("-cfl", "--cfl", "print CFL number each dT");
        const bool pl2norm = args.getflag("-l2", "--l2norm", "print L2Norm(u) each dT");
        const bool pchnorm = args.getbool("-ch", "--chebyNorm", true, "print chebyNorm(u) each dT");
        const bool pdissip = args.getflag("-D", "--dissipation", "print dissipation each dT");
        const bool pshear = args.getflag("-I", "--input", "print wall shear power input each dT");
        const Real rhomin = args.getreal("-rm", "--rhomin", 0.0, "stop if chebyNorm(rho) < rhomin");

        const int saveint = args.getint("-s", "--saveinterval", 1, "save fields every s dT");

        const int nproc0 =
            args.getint("-np0", "--nproc0", 0, "number of MPI-processes for transpose/number of parallel ffts");
        const int nproc1 = args.getint("-np1", "--nproc1", 0, "number of MPI-processes for one fft");
        const string rhofile = args.getstr(1, "<flowfield>", "initial concentration/density field");
        const string velfile = args.getstr(2, "<flowfield>", "precomputed velocity field");

        args.check();
        args.save("./");
        mkdir(outdir);
        args.save(outdir);
        flags.save(outdir);

        CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);

        printout("Constructing u,q, and optimizing FFTW...");

        // load concentration/density field
        FlowField rhoin(rhofile, cfmpi);
        Real vs_o_kappa = flags.vs / flags.kappa;
        BoundaryCond bc(Mixed, 0.0, 1.0 + vs_o_kappa, vs_o_kappa);
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

        const int Nx = rho.Nx();
        const int Ny = rho.Ny();
        const int Nz = rho.Nz();
        const Real Lx = rho.Lx();
        const Real Lz = rho.Lz();
        const Real a = rho.a();
        const Real b = rho.b();
 
        // load velocity field, upscaling 
        // if necessary for the later dot product in physical space
        FlowField uin(velfile, cfmpi);
        if (Lx != uin.Lx() || Lz != uin.Lz() || a != uin.a() || b != uin.b()) {
            cout << "velocity field geometry must match rho geometry\n";
            cfMPI_Finalize();
            return 1;
        }
        FlowField u(Nx, Ny, Nz, uin.Nd(), Lx, Lz, a, b, uin.BC(), rho.cfmpi());
        u.interpolate(uin);

        const bool inttime =
            (abs(saveint * dt.dT() - int(saveint * dt.dT())) < 1e-12) && (abs(flags.t0 - int(flags.t0)) < 1e-12)
                ? true
                : false;

        cout << "Uwall == " << flags.uupperwall << endl;
        cout << "Wwall == " << flags.wupperwall << endl;
        cout << "dnsflags == " << flags << endl;
        cout << "constructing DNS..." << endl;
        vector<FlowField> fields = {rho, u};
        DNS dns(fields, flags);

        dns.Ubase().save(outdir + "Ubase");
        dns.Wbase().save(outdir + "Wbase");

        ChebyCoeff Ubase = laminarProfile(flags, u.a(), u.b(), u.Ny());

        ios::openmode openflag = (flags.t0 > 0) ? ios::app : ios::out;

        ofstream eout, x0out;
        openfile(eout, outdir + "energy.asc", openflag);
        eout << fieldstatsheader_t() << endl;

        int i = 0;
        for (Real t = flags.t0; t <= flags.T; t += dt.dT()) {
            string s;
            s = printdiagnostics(fields[0], fields[1], dns, t, dt, flags.nu, rhomin, 
                dt.variable(), pl2norm, pchnorm, pdissip,
                pshear, pcfl);

            cout << s;

            if (saveint != 0 && i % saveint == 0) {
                fields[0].save(outdir + label + t2s(t, inttime));
                //if (savep)
                //    fields[1].save(outdir + "p" + t2s(t, inttime));
            }
            i++;

            dns.advance(fields, dt.n());

            if (dt.variable() &&
                dt.adjust(dns.CFL(fields[1])))  // TODO: dt.variable()==true is checked twice here, remove it.
                dns.reset_dt(dt);
        }
        cout << "done!" << endl;
    }
    cfMPI_Finalize();
}

string printdiagnostics(FlowField& rho, FlowField& u, const DNS& dns, Real t, const TimeStep& dt, Real nu, Real rhomin, bool vardt,
                        bool pl2norm, bool pchnorm, bool pdissip, bool pshear, //bool pUbulk, bool pubulk,
                        //bool pdPdx, 
                        bool pcfl) {
    // Printing diagnostics
    stringstream sout;
    sout << "           t == " << t << endl;
    if (vardt)
        sout << "          dt == " << Real(dt) << endl;
    if (pl2norm)
        sout << "   L2Norm(rho) == " << L2Norm(rho) << endl;

    if (pchnorm || rhomin != 0.0) {
        Real chnorm = chebyNorm(rho);
        sout << "chebyNorm(rho) == " << chnorm << endl;
        if (chnorm < rhomin) {
            cout << "Exiting: chebyNorm(rho) < rhomin." << endl;
            exit(0);
        }
    }

    Real h = 0.5 * (u.b() - u.a());
    rho -= dns.Ubase();
    if (pl2norm)
        sout << "   energy(rho+rBase) == " << 0.5 * L2Norm(rho) << endl;
    if (pdissip)
        sout << "   dissip(rho+rBase) == " << dissipation(rho) << endl;
    if (pshear)
        sout << "wallshear(rho+rBase) == " << abs(wallshearLower(rho)) + abs(wallshearUpper(rho)) << endl;
    rho += dns.Ubase();
    if (pl2norm)
        sout << "     L2Norm(rho) == " << L2Norm(rho) << endl;
    if (pl2norm)
        sout << "   L2Norm3d(rho) == " << L2Norm3d(rho) << endl;

    Real cfl = dns.CFL(u);
    if (u.taskid() == u.task_coeff(0, 0)) {
        ChebyCoeff U = dns.Ubase();
        ChebyCoeff W = dns.Wbase();

        U.makeSpectral();
        U += Re(u.profile(0, 0, 0));
        Real Ucenter = U.eval(0.5 * (u.a() + u.b()));
        Real Uwall = pythag(0.5 * (U.eval_b() - U.eval_a()), 0.5 * (W.eval_b() - W.eval_a()));
        Real Umean = U.mean();
        sout << "        1/nu == " << 1 / nu << endl;
        sout << "  Uwall h/nu == " << Uwall * h / nu << endl;
        sout << "  Umean h/nu == " << Umean << " * " << h << " / " << nu << endl;
        sout << "  Umean h/nu == " << Umean * h / nu << endl;
        sout << "Ucenter h/nu == " << Ucenter * h / nu << endl;
    }
    sout << "         CFL == " << cfl << endl;
    return sout.str();
}
