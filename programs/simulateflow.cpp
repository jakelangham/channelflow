/**
 * This file is a part of channelflow version 2.0 https://channelflow.ch.
 * License is GNU GPL version 2 or later: ./LICENCE
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
#include "channelflow/poissonsolver.h"
#include "channelflow/symmetry.h"
#include "channelflow/tausolver.h"
#include "channelflow/utilfuncs.h"

using namespace std;
using namespace chflow;

string printdiagnostics(FlowField& u, const DNS& dns, Real t, const TimeStep& dt, Real nu, Real umin, bool vardt,
                        bool pl2norm, bool pchnorm, bool pdissip, bool pshear, bool pdiverge, bool pUbulk, bool pubulk,
                        bool pdPdx, bool pcfl);
void saveprofiles(const FlowField& u, const FlowField& p, const ChebyCoeff& Ubase, const string outdir, Real t, const bool inttime);

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
        const bool savep = args.getflag("-sp", "--savepressure", "save pressure fields");

        const bool pcfl = args.getflag("-cfl", "--cfl", "print CFL number each dT");
        const bool pl2norm = args.getflag("-l2", "--l2norm", "print L2Norm(u) each dT");
        const bool pchnorm = args.getbool("-ch", "--chebyNorm", true, "print chebyNorm(u) each dT");
        const bool pdissip = args.getflag("-D", "--dissipation", "print dissipation each dT");
        const bool pshear = args.getflag("-I", "--input", "print wall shear power input each dT");
        const bool pdiverge = args.getflag("-dv", "--divergence", "print divergence each dT");
        const bool pubulk = args.getflag("-u", "--ubulk", "print ubulk each dT");
        const bool pUbulk = args.getflag("-Up", "--Ubulk-print", "print Ubulk each dT");
        const bool pdPdx = args.getflag("-p", "--pressure", "print pressure gradient each dT");
        const Real umin = args.getreal("-u", "--umin", 0.0, "stop if chebyNorm(u) < umin");

        const Real ecfmin = args.getreal("-e", "--ecfmin", 0.0, "stop if Ecf(u) < ecfmin");
        const int saveint = args.getint("-s", "--saveinterval", 1, "save fields every s dT");
        const int saveprofint = args.getint("-sprof", "--saveprofinterval", 1, "save profile data every s dT");

        const int nproc0 =
            args.getint("-np0", "--nproc0", 0, "number of MPI-processes for transpose/number of parallel ffts");
        const int nproc1 = args.getint("-np1", "--nproc1", 0, "number of MPI-processes for one fft");
        const string uname = args.getstr(1, "<flowfield>", "initial condition");

        args.check();
        args.save("./");
        mkdir(outdir);
        args.save(outdir);
        flags.save(outdir);

        CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);

        printout("Constructing u,q, and optimizing FFTW...");
        FlowField u(uname, cfmpi);
        // JL check if there's a density field. If not, add a zero density field
        Real vs_o_kappa = flags.vs / flags.kappa;
        //BoundaryCond bc(DiriRobin, u.Mx(), u.Mz(), 0.0, 1.0 + vs_o_kappa, vs_o_kappa);
        //BoundaryCond bc(NeumRobin, u.Mx(), u.Mz(), 1.0, 1.0 + vs_o_kappa, vs_o_kappa);
        BoundaryCond bc(NoFlux, u.Mx(), u.Mz(), 1.0 - vs_o_kappa, 1.0 + vs_o_kappa, vs_o_kappa);
        //BoundaryCond bc(Diri, 0.0, 0.0);
        FlowField u_with_density(u.Nx(), u.Ny(), u.Nz(), 4, u.Lx(), u.Lz(), u.a(), u.b(), bc, cfmpi);
        if (u.padded())
            u_with_density.setPadded(true);
        vector<int> vel_indices = {0, 1, 2};
        vector<int> all_indices = {0, 1, 2, 3};
        if (u.Nd() == 3) {
            u_with_density.copySubfields(u, vel_indices, vel_indices);
        } else {
            u_with_density.copySubfields(u, all_indices, all_indices);
        }    

        const int Nx = u.Nx();
        const int Ny = u.Ny();
        const int Nz = u.Nz();
        const Real Lx = u.Lx();
        const Real Lz = u.Lz();
        const Real a = u.a();
        const Real b = u.b();

        FlowField q(Nx, Ny, Nz, 1, Lx, Lz, a, b, BoundaryCond(), cfmpi);
        FlowField p(Nx, Ny, Nz, 1, Lx, Lz, a, b, BoundaryCond(), cfmpi);
        const bool inttime =
            (abs(saveint * dt.dT() - int(saveint * dt.dT())) < 1e-12) && (abs(flags.t0 - int(flags.t0)) < 1e-12)
                ? true
                : false;

        cout << "Uwall == " << flags.uupperwall << endl;
        cout << "Wwall == " << flags.wupperwall << endl;
        cout << "dnsflags == " << flags << endl;
        cout << "constructing DNS..." << endl;
        DNS dns({u_with_density, q}, flags);
        //     u.setnu (flags.nu);
        //
        cout << "uwd " << u_with_density.taskid() << endl;

        dns.Ubase().save(outdir + "Ubase");
        dns.Wbase().save(outdir + "Wbase");

        ChebyCoeff Ubase = laminarProfile(flags, u.a(), u.b(), u.Ny());
        // Ubase.save("Ubase2");

        PressureSolver psolver(u_with_density, dns.Ubase(), dns.Wbase(), flags.nu, flags.Vsuck, flags.nonlinearity);
        psolver.solve(q, u_with_density);
        vector<FlowField> fields = {u_with_density, q};

        ios::openmode openflag = (flags.t0 > 0) ? ios::app : ios::out;

        ofstream eout, x0out, concout;
        openfile(eout, outdir + "energy.asc", openflag);
        eout << fieldstatsheader_t() << endl;
        openfile(concout, outdir + "cstats.asc", openflag);
        concout << concstatsheader_t() << endl;

        FlowField u0, du, tmp;
        FlowField vel(u.Nx(), u.Ny(), u.Nz(), 4, u.Lx(), u.Lz(), u.a(), u.b(), bc, cfmpi);
        FlowField uc(u.Nx(), u.Ny(), u.Nz(), 4, u.Lx(), u.Lz(), u.a(), u.b(), bc, cfmpi);

        int i = 0;
        for (Real t = flags.t0; t <= flags.T; t += dt.dT()) {
            string s;
            s = printdiagnostics(fields[0], dns, t, dt, flags.nu, umin, dt.variable(), pl2norm, pchnorm, pdissip,
                                 pshear, pdiverge, pUbulk, pubulk, pdPdx, pcfl);

            if (ecfmin > 0 && Ecf(fields[0]) < ecfmin) {
                cferror("Ecf < ecfmin == " + r2s(ecfmin) + ", exiting");
            }

            cout << s;
            // fields[0] += dns.Ubase(); //////////////////// ONLY
            s = fieldstats_t(fields[0], t);
            // fields[0] -= dns.Ubase(); //////////////////// ONLY
            eout << s << endl;

            s = concstats_t(fields[0], vel, uc, t);
            concout << s << endl;

            if (saveprofint != 0 && i % saveprofint == 0) {
                // save <u>
                dns.uq2p(fields[0], fields[1], p);
                saveprofiles(fields[0], p, dns.Ubase(), outdir, t, inttime);
            }
            if (saveint != 0 && i % saveint == 0) {
                fields[0].save(outdir + label + t2s(t, inttime));
                if (savep) {
                    fields[1].save(outdir + "p" + t2s(t, inttime));
                }
            }
            i++;

            dns.advance(fields, dt.n());

            if (dt.variable() &&
                dt.adjust(dns.CFL(fields[0])))  // TODO: dt.variable()==true is checked twice here, remove it.
                dns.reset_dt(dt);
        }
        cout << "done!" << endl;
    }
    cfMPI_Finalize();
}

string printdiagnostics(FlowField& u, const DNS& dns, Real t, const TimeStep& dt, Real nu, Real umin, bool vardt,
                        bool pl2norm, bool pchnorm, bool pdissip, bool pshear, bool pdiverge, bool pUbulk, bool pubulk,
                        bool pdPdx, bool pcfl) {
    // Printing diagnostics
    stringstream sout;
    sout << "           t == " << t << endl;
    if (vardt)
        sout << "          dt == " << Real(dt) << endl;
    if (pl2norm)
        sout << "   L2Norm(u) == " << L2Norm(u) << endl;

    if (pchnorm || umin != 0.0) {
        Real chnorm = chebyNorm(u);
        sout << "chebyNorm(u) == " << chnorm << endl;
        if (chnorm < umin) {
            cout << "Exiting: chebyNorm(u) < umin." << endl;
            exit(0);
        }
    }
    Real h = 0.5 * (u.b() - u.a());
    u += dns.Ubase();
    if (pl2norm)
        sout << "   energy(u+U) == " << 0.5 * L2Norm(u) << endl;
    if (pdissip)
        sout << "   dissip(u+U) == " << dissipation(u) << endl;
    if (pshear)
        sout << "wallshear(u+U) == " << abs(wallshearLower(u)) + abs(wallshearUpper(u)) << endl;
    if (pdiverge)
        sout << "  divNorm(u+U) == " << divNorm(u) << endl;
    if (pUbulk)
        sout << "mean u+U Ubulk == " << dns.Ubulk() << endl;
    u -= dns.Ubase();
    if (u.taskid() == u.task_coeff(0, 0)) {
        if (pubulk)
            sout << "         ubulk == " << Re(u.profile(0, 0, 0)).mean() << endl;
    }
    if (pdPdx)
        sout << "          dPdx == " << dns.dPdx() << endl;
    if (pl2norm)
        sout << "     L2Norm(u) == " << L2Norm(u) << endl;
    if (pl2norm)
        sout << "   L2Norm3d(u) == " << L2Norm3d(u) << endl;

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
        sout << "  Ubulk h/nu == " << dns.Ubulk() * h / nu << endl;
        sout << "  Umean h/nu == " << Umean << " * " << h << " / " << nu << endl;
        sout << "  Umean h/nu == " << Umean * h / nu << endl;
        sout << " Uparab h/nu == " << 1.5 * dns.Ubulk() * h / nu << endl;
        sout << "Ucenter h/nu == " << Ucenter * h / nu << endl;
    }
    sout << "         CFL == " << cfl << endl;
    return sout.str();
}
void saveprofiles(const FlowField& u, const FlowField& p, const ChebyCoeff& Ubase, const string outdir, Real t, const bool inttime) {
    ChebyCoeff uprof = Re(u.profile(0, 0, 0));
    uprof += Ubase;
    ChebyCoeff wprof = Re(u.profile(0, 0, 2));
    ChebyCoeff cprof = Re(u.profile(0, 0, 3));
    ChebyCoeff dcxzdy = diff(cprof);
    ChebyCoeff d2uxzdy2 = diff(uprof, 2);

    //uprof.makePhysical(); wprof.makePhysical(); cprof.makePhysical();
    if (u.taskid() == u.task_coeff(0, 0)) {
        uprof.save(outdir + "uxz" + t2s(t, inttime), Physical);
        wprof.save(outdir + "wxz" + t2s(t, inttime), Physical);
        cprof.save(outdir + "cxz" + t2s(t, inttime), Physical);
        dcxzdy.save(outdir + "dcxzdy" + t2s(t, inttime), Physical);
        d2uxzdy2.save(outdir + "d2uxzdy2" + t2s(t, inttime), Physical);
    }

    FlowField e = TKE(u);
    ChebyCoeff TKEprof = Re(e.profile(0, 0, 0));
    ChebyCoeff visc_diffusion_prof = diff2(TKEprof);
    if (u.taskid() == u.task_coeff(0, 0)) {
//        TKEprof.save(outdir + "tkexz" + t2s(t, inttime), Physical);
        visc_diffusion_prof.save(outdir + "viscdiffxz" + t2s(t, inttime), Physical);
    }
    FlowField upup = Reynolds_stress_uu(u);
    FlowField upvp = Reynolds_stress_uv(u);
    FlowField upwp = Reynolds_stress_uw(u);
    FlowField vpvp = Reynolds_stress_vv(u);
    FlowField vpwp = Reynolds_stress_vw(u);
    FlowField wpwp = Reynolds_stress_ww(u);
    FlowField upcp = turbulent_sediment_flux_uc(u);
    FlowField vpcp = turbulent_sediment_flux_vc(u);
    FlowField cpcp = turbulent_sediment_flux_cc(u);
//    FlowField uskew = u_skewness(u);
//    FlowField ukurt = u_kurtosis(u);
//    FlowField cskew = c_skewness(u);
//    FlowField ckurt = c_kurtosis(u);
    ChebyCoeff upup_prof = Re(upup.profile(0, 0, 0));
    ChebyCoeff upvp_prof = Re(upvp.profile(0, 0, 0));
    ChebyCoeff upwp_prof = Re(upwp.profile(0, 0, 0));
    ChebyCoeff vpvp_prof = Re(vpvp.profile(0, 0, 0));
    ChebyCoeff vpwp_prof = Re(vpwp.profile(0, 0, 0));
    ChebyCoeff wpwp_prof = Re(wpwp.profile(0, 0, 0));
    ChebyCoeff upcp_prof = Re(upcp.profile(0, 0, 0));
    ChebyCoeff vpcp_prof = Re(vpcp.profile(0, 0, 0));
    ChebyCoeff cpcp_prof = Re(cpcp.profile(0, 0, 0));
//    ChebyCoeff uskew_prof = Re(uskew.profile(0, 0, 0));
//    ChebyCoeff ukurt_prof = Re(ukurt.profile(0, 0, 0));
//    ChebyCoeff cskew_prof = Re(cskew.profile(0, 0, 0));
//    ChebyCoeff ckurt_prof = Re(ckurt.profile(0, 0, 0));
    ChebyCoeff dupvpxzdy_prof = diff(upvp_prof);
    ChebyCoeff diff_cpcp_prof = diff(cpcp_prof);
    ChebyCoeff diff2_cpcp_prof = diff(cpcp_prof, 2);
    if (u.taskid() == u.task_coeff(0, 0)) {
        upup_prof.save(outdir + "upupxz" + t2s(t, inttime), Physical);
        upvp_prof.save(outdir + "upvpxz" + t2s(t, inttime), Physical);
        dupvpxzdy_prof.save(outdir + "dupvpxzdy" + t2s(t, inttime), Physical);
        upwp_prof.save(outdir + "upwpxz" + t2s(t, inttime), Physical);
        vpvp_prof.save(outdir + "vpvpxz" + t2s(t, inttime), Physical);
        vpwp_prof.save(outdir + "vpwpxz" + t2s(t, inttime), Physical);
        wpwp_prof.save(outdir + "wpwpxz" + t2s(t, inttime), Physical);
        upcp_prof.save(outdir + "upcpxz" + t2s(t, inttime), Physical);
        vpcp_prof.save(outdir + "vpcpxz" + t2s(t, inttime), Physical);
//        cpcp_prof.save(outdir + "cpcpxz" + t2s(t, inttime), Physical);
//        uskew_prof.save(outdir + "uskewxz" + t2s(t, inttime), Physical);
//        ukurt_prof.save(outdir + "ukurtxz" + t2s(t, inttime), Physical);
//        cskew_prof.save(outdir + "cskewxz" + t2s(t, inttime), Physical);
//        ckurt_prof.save(outdir + "ckurtxz" + t2s(t, inttime), Physical);
        diff_cpcp_prof.save(outdir + "dcpcpxzdy" + t2s(t, inttime), Physical);
        diff2_cpcp_prof.save(outdir + "d2cpcpxzdy" + t2s(t, inttime), Physical);
    }
    ChebyCoeff duxzdy = diff(uprof);
    upvp_prof.makePhysical(); duxzdy.makePhysical();
    ChebyCoeff prod = upvp_prof;
    prod *= duxzdy;
    if (u.taskid() == u.task_coeff(0, 0)) {
        prod.save(outdir + "prod" + t2s(t, inttime), Physical);
    }

    FlowField dissipxz = dissipation_field(u);
    FlowField sed_dissipxz = sed_dissipation_field(u);
    ChebyCoeff dissipxz_prof = Re(dissipxz.profile(0, 0, 0));
    ChebyCoeff sed_dissipxz_prof = Re(sed_dissipxz.profile(0, 0, 0));
    if (u.taskid() == u.task_coeff(0, 0)) {
        dissipxz_prof.save(outdir + "dissipxz" + t2s(t, inttime), Physical);
        sed_dissipxz_prof.save(outdir + "seddissipxz" + t2s(t, inttime), Physical);
    }
    FlowField turb_conv = turbulent_convection(u);
    FlowField sed_conv = sediment_fluctuation_convection(u);
    ChebyCoeff turb_conv_prof = Re(turb_conv.profile(0, 0, 0));
    ChebyCoeff turb_conv_prof_diff = diff(turb_conv_prof);
    ChebyCoeff sed_conv_prof = Re(sed_conv.profile(0, 0, 0));
    ChebyCoeff sed_conv_prof_diff = diff(sed_conv_prof);
    if (u.taskid() == u.task_coeff(0, 0)) {
        turb_conv_prof_diff.save(outdir + "tconvxz" + t2s(t, inttime), Physical);
        sed_conv_prof_diff.save(outdir + "sedconvxz" + t2s(t, inttime), Physical);
    }
    FlowField p_transp = pressure_transport(u, p);
    ChebyCoeff p_transp_prof = Re(p_transp.profile(0, 0, 0));
    ChebyCoeff diffvp = diff(p_transp_prof);
    if (u.taskid() == u.task_coeff(0, 0)) {
        diffvp.save(outdir + "ptranspxz" + t2s(t, inttime), Physical);
    }
    FlowField intermittency = intermittency_factor(u, 1e-4);
    ChebyCoeff intermittency_prof = Re(intermittency.profile(0, 0, 0));
    if (u.taskid() == u.task_coeff(0, 0)) {
        intermittency_prof.save(outdir + "intermittencyxz" + t2s(t, inttime), Physical);
    }

    // xz slices
//    string preface = outdir + "slice_t" + t2s(t, inttime);
//    int ny_mid = (u.Ny() - 1) / 2;
//    // For Ny = 141, these give y ~= -0.3, -0.6, -0.9
//    //int ny_L1 = 84, ny_L2 = 99, ny_L3 = 120;
//    //int ny_L2 = 99, ny_L3 = 120;
//    int ny_L2 = 35, ny_L3 = 42;
//    u.saveSlice(1, 0, ny_mid, preface + "_u_mid_xz", 1, 1, 1);
//    //u.saveSlice(1, 0, ny_L1, preface + "_u_L1_xz", 1, 1, 1);
//    u.saveSlice(1, 0, ny_L2, preface + "_u_L2_xz", 1, 1, 1);
//    u.saveSlice(1, 0, ny_L3, preface + "_u_L3_xz", 1, 1, 1);
//
//    vpcp.saveSlice(1, 0, ny_mid, preface + "_vpcp_mid_xz", 1, 1, 1);
//    //vpcp.saveSlice(1, 0, ny_L1, preface + "_vpcp_L1_xz", 1, 1, 1);
//    vpcp.saveSlice(1, 0, ny_L2, preface + "_vpcp_L2_xz", 1, 1, 1);
//    vpcp.saveSlice(1, 0, ny_L3, preface + "_vpcp_L3_xz", 1, 1, 1);

    //upvp.saveSlice(1, 0, ny_mid, preface + "_upvp_mid_xz", 1, 1, 1);
    //upvp.saveSlice(1, 0, ny_L1, preface + "_upvp_L1_xz", 1, 1, 1);
    //upvp.saveSlice(1, 0, ny_L2, preface + "_upvp_L2_xz", 1, 1, 1);
    //upvp.saveSlice(1, 0, ny_L3, preface + "_upvp_L3_xz", 1, 1, 1);
}
