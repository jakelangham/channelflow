/**
 * The periodic response of a passive settling scalar (Ri = 0) to a periodic orbit of the velocity.
 *
 * The velocity is integrated from the same u(0) every period, so it stays on the (unstable) orbit; the
 * scalar is integrated continuously, with the velocity prescribed, until c(T) = c(0) to a tolerance.
 * The scalar can have its own grid (e.g. a larger Ny for the settling boundary layer at high
 * vs/kappa): each step the velocity is interpolated onto it spectrally.
 *
 * This file is a part of channelflow version 2.0 https://channelflow.ch.
 * License is GNU GPL version 2 or later: ./LICENCE
 */

#include <cmath>
#include <fstream>
#include <iomanip>
#include <memory>
#include <iostream>
#include <string>

#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/poissonsolver.h"
#include "channelflow/utilfuncs.h"
#include "modules/settling/settlingdns.h"
#include "modules/settling/settlingdsi.h"
#include "modules/settling/settlingnse.h"

using namespace std;
using namespace chflow;

// x-z mean profiles (spectral) of C, (u + U) C and v C on c's grid, where C = c - y is the
// concentration (the stored scalar is C + y) and u the velocity perturbation interpolated onto c's grid.
// The products are formed at c's collocation points, where u is exact, so the x-z means are exact.
void fluxProfiles(const FlowField& u, const FlowField& c, const ChebyCoeff& Ubase, FlowField& tmp, ChebyCoeff& C,
                  ChebyCoeff& uC, ChebyCoeff& vC) {
    FlowField ut(u), Ct(c);
    ut.makeSpectral();
    Ct.makeSpectral();
    if (ut.taskid() == ut.task_coeff(0, 0)) {
        ChebyCoeff U(Ubase);
        U.makeSpectral();
        for (int ny = 0; ny < ut.Ny(); ++ny)
            ut.cmplx(0, ny, 0, 0) += Complex(U[ny], 0.0);
        // y = (b + a)/2 T0 + (b - a)/2 T1
        Ct.cmplx(0, 0, 0, 0) -= Complex(0.5 * (c.b() + c.a()), 0.0);
        Ct.cmplx(0, 1, 0, 0) -= Complex(0.5 * (c.b() - c.a()), 0.0);
    }
    C = Re(Ct.profile(0, 0, 0));

    ut.makePhysical();
    Ct.makePhysical();
    if (!tmp.geomCongruent(c) || tmp.Nd() != 2)
        tmp.resize(c.Nx(), c.Ny(), c.Nz(), 2, c.Lx(), c.Lz(), c.a(), c.b(), c.cfmpi());
    tmp.setState(Physical, Physical);
    for (lint ny = ut.nylocmin(); ny < ut.nylocmax(); ++ny)
        for (lint nx = ut.nxlocmin(); nx < ut.nxlocmin() + ut.Nxloc(); ++nx)
            for (lint nz = 0; nz < ut.Nz(); ++nz) {
                tmp(nx, ny, nz, 0) = ut(nx, ny, nz, 0) * Ct(nx, ny, nz, 0);
                tmp(nx, ny, nz, 1) = ut(nx, ny, nz, 1) * Ct(nx, ny, nz, 0);
            }
    tmp.makeSpectral();
    uC = Re(tmp.profile(0, 0, 0));
    vC = Re(tmp.profile(0, 0, 1));
}

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        WriteProcessInfo(argc, argv);
        string purpose(
            "periodic response of a passive settling scalar (Ri = 0) to a periodic orbit u(t) of period T: "
            "integrate u from u(0) every period and c continuously until c(T) = c(0)");
        ArgList args(argc, argv, purpose);

        DNSFlags flags(args);

        args.section("Program options");
        const string outdir = args.getpath("-o", "--outdir", "scalarorbit/", "output directory");
        const string cname = args.getstr("-c", "--scalar", "",
                                         "initial scalar (1 or 4 components; default: component 3 of <velocity>)");
        const int cNx = args.getint("-cNx", "--cNx", 0, "x gridpoints for the scalar (0 == as the velocity)");
        const int cNy = args.getint("-cNy", "--cNy", 0, "y gridpoints for the scalar (0 == as the velocity)");
        const int cNz = args.getint("-cNz", "--cNz", 0, "z gridpoints for the scalar (0 == as the velocity)");
        const int Nperiods = args.getint("-Np", "--Nperiods", 20, "maximum number of periods");
        const Real tol = args.getreal("-tol", "--tolerance", 1e-10, "stop when L2Dist(c(T), c(0))/L2Norm(c) < tol");
        const int Nsave = args.getint("-ns", "--Nsave", 0,
                                      "save this many states per period (overwritten each period; 0 == none)");
        const int Navg = args.getint("-na", "--Naverage", 400,
                                     "samples per period for the period-averaged flux profiles (trapezoid rule)");
        const int Nex = args.getint("-ex", "--extrapolate", 0,
                                    "every this many periods, extrapolate c past its dominant decaying mode "
                                    "(Aitken; >= 2, 0 == never)");
        const Real exmin = args.getreal("-exmin", "--extrapolateMin", 1e-4,
                                        "extrapolate only while |c(T) - c(0)|/|c(0)| > exmin: each extrapolation "
                                        "restarts the multistep history, which costs ~2e-6 at vs/kappa = 32");
        const int nproc0 = args.getint("-np0", "--nproc0", 0, "number of MPI-processes for transpose");
        const int nproc1 = args.getint("-np1", "--nproc1", 0, "number of MPI-processes for one fft");
        const string uname = args.getstr(1, "<velocity>", "velocity at t = 0 on the orbit (3 or 4 components)");

        args.check();
        mkdir(outdir);
        args.save(outdir);
        flags.save(outdir);

        const Real T = flags.T;
        const vector<SettlingScalar> scalars = settlingScalars(flags);
        if (scalars.size() != 1 || scalars[0].Ri != 0.0)
            cferror("scalarorbit : needs one passive scalar (-ri 0)");

        CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);

        // velocity at t = 0, and the scalar's initial condition
        FlowField uin(uname, cfmpi);
        FlowField u0, cin;
        if (uin.Nd() == 4)
            splitVelocityScalar(uin, u0, cin);
        else
            u0 = uin;
        if (cname != "") {
            FlowField cfile(cname, cfmpi);
            if (cfile.Nd() == 4) {
                FlowField tmp;
                splitVelocityScalar(cfile, tmp, cin);
            } else
                cin = cfile;
        }
        if (cin.Nd() != 1)
            cferror("scalarorbit : no initial scalar (give -c, or a 4-component <velocity>)");
        u0.makeSpectral();

        // the scalar's grid, the velocity interpolated onto it, and the scalar's initial condition there
        const int Nxc = cNx > 0 ? cNx : u0.Nx();
        const int Nyc = cNy > 0 ? cNy : u0.Ny();
        const int Nzc = cNz > 0 ? cNz : u0.Nz();
        FlowField uc(Nxc, Nyc, Nzc, 3, u0.Lx(), u0.Lz(), u0.a(), u0.b(), cfmpi);
        FlowField c(Nxc, Nyc, Nzc, 1, u0.Lx(), u0.Lz(), u0.a(), u0.b(), cfmpi);
        FlowField qc(Nxc, Nyc, Nzc, 1, u0.Lx(), u0.Lz(), u0.a(), u0.b(), cfmpi);
        c.interpolate(cin);
        projectScalar(c, scalars[0]);
        cout << "velocity grid " << u0.Nx() << " x " << u0.Ny() << " x " << u0.Nz() << ", scalar grid " << Nxc
             << " x " << Nyc << " x " << Nzc << endl;

        // the velocity's own start-up: the scalar's (multistep only) need not suit it
        DNSFlags uflags = flags;
        uflags.initstepping = SMRK2;

        // a whole number of steps per period, and of steps between saves and between flux samples
        auto gcd = [](int p, int q) {
            while (q != 0) {
                const int r = p % q;
                p = q;
                q = r;
            }
            return p;
        };
        const int Nmult = Nsave > 0 ? Nsave / gcd(Nsave, Navg) * Navg : Navg;
        auto stepsPerPeriod = [&](Real dtwant) { return Nmult * int(ceil(T / dtwant / Nmult - 1e-9)); };
        int N = stepsPerPeriod(flags.dt);

        // Check the CFL numbers at t = 0 before the first step (the scalar's need not be linear in dt,
        // so iterate); only reduce dt here, since the velocity's varies along the orbit
        if (flags.variabledt) {
            uc.interpolate(u0);
            for (int iter = 0; iter < 10; ++iter) {
                DNSFlags f0 = flags;
                f0.dt = T / N;
                FlowField q0(u0.Nx(), u0.Ny(), u0.Nz(), 1, u0.Lx(), u0.Lz(), u0.a(), u0.b(), cfmpi);
                const Real cfl = Greater(DNS({u0, q0}, f0).CFL(u0), SettlingDNS({uc, c, qc}, f0, true).CFL(uc));
                cout << "t = 0: dt == " << f0.dt << ", CFL == " << cfl << endl;
                if (cfl <= flags.CFLmax)
                    break;
                N = stepsPerPeriod(f0.dt * 0.5 * (flags.CFLmin + flags.CFLmax) / cfl);
            }
        }

        FlowField u4(Nxc, Nyc, Nzc, 4, u0.Lx(), u0.Lz(), u0.a(), u0.b(), cfmpi);
        FlowField cstart;
        unique_ptr<SettlingDNS> cdns;
        vector<FlowField> cfields;
        Real dtc = 0.0;

        ofstream conv((outdir + "convergence.asc").c_str());
        conv << "% period  dt  steps  maxCFL(velocity)  maxCFL(scalar)  L2Dist(u(T),u(0))/L2Norm(u(0))  "
                "L2Dist(c(T),c(0))/L2Norm(c(0))  <(u+U)C>  <vC>  max|F|/max|vs<C>|"
             << endl;
        const Real vs = scalars[0].vs;
        const Real kappa = scalars[0].kappa;
        FlowField ftmp;
        ChebyCoeff C, uC, vC;

        // Extrapolation: one period maps c to P c + b (linear, the velocity fixed), so after a few
        // periods the error is ~ lambda^n v, lambda the dominant eigenvalue of P. From the last two
        // period-to-period changes d1, d2, lambda ~ <d2, d1>/<d1, d1> and c* ~ c + lambda/(1 - lambda) d2.
        // Needs two free periods since the last change of c (start or extrapolation).
        FlowField cstartPrev;
        int sinceJump = 0;

        for (int period = 0; period < Nperiods; ++period) {
            const Real dt = T / N;
            uflags.dt = dt;
            FlowField u = u0;
            FlowField q(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), cfmpi);
            DNS udns({u, q}, uflags);
            PressureSolver psolver(u, udns.Ubase(), udns.Wbase(), uflags.nu, uflags.Vsuck, uflags.nonlinearity);
            psolver.solve(q, u);
            vector<FlowField> ufields = {u, q};

            // the scalar's DNS is kept across periods while dt is unchanged (so c's history is kept)
            if (dt != dtc) {
                if (cdns)
                    c = cfields[1];
                DNSFlags cflags = flags;
                cflags.dt = dt;
                uc.interpolate(u0);
                cfields = {uc, c, qc};
                cdns.reset(new SettlingDNS(cfields, cflags, true));
                dtc = dt;
            }
            cstartPrev = cstart;
            cstart = cfields[1];

            Real cflu = 0.0, cflc = 0.0;
            const int checkint = Greater(1, int(round(1.0 / dt)));  // CFL check about once per time unit
            ChebyCoeff Cavg(Nyc, c.a(), c.b(), Spectral), uCavg(Cavg), vCavg(Cavg);
            for (int n = 0; n < N; ++n) {
                cfields[0].interpolate(ufields[0]);
                if (n % (N / Navg) == 0) {
                    fluxProfiles(cfields[0], cfields[1], cdns->Ubase(), ftmp, C, uC, vC);
                    Cavg += C;
                    uCavg += uC;
                    vCavg += vC;
                }
                if (n % checkint == 0) {
                    cflu = Greater(cflu, udns.CFL(ufields[0]));
                    cflc = Greater(cflc, cdns->CFL(cfields[0]));
                }
                if (Nsave > 0 && n % (N / Nsave) == 0) {
                    combineVelocityScalar(cfields[0], cfields[1], u4);
                    u4.setPadded(flags.dealias_xz());
                    u4.save(outdir + "u" + i2s(n / (N / Nsave)));
                }
                cdns->advance(cfields, 1);
                udns.advance(ufields, 1);
            }
            cfields[0].interpolate(ufields[0]);  // u(T), for the saved state

            const Real du = L2Dist(ufields[0], u0) / L2Norm(u0);
            const Real dc = L2Dist(cfields[1], cstart) / L2Norm(cstart);
            cout << "  L2Norm(c(0)) == " << L2Norm(cstart) << " (state " << cstart.xzstate() << cstart.ystate()
                 << "), L2Norm(c(T)) == " << L2Norm(cfields[1]) << " (state " << cfields[1].xzstate()
                 << cfields[1].ystate() << ")" << endl;

            // period averages (periodic trapezoid rule), and the total vertical flux
            // F = <vC> - vs <C> - kappa d<C>/dy, which vanishes at every y for a converged periodic state
            Cavg *= 1.0 / Navg;
            uCavg *= 1.0 / Navg;
            vCavg *= 1.0 / Navg;
            ChebyCoeff dCdy = diff(Cavg);
            const Real uCbulk = uCavg.mean();
            const Real vCbulk = vCavg.mean();
            Cavg.makePhysical();
            uCavg.makePhysical();
            vCavg.makePhysical();
            dCdy.makePhysical();
            Vector y = chebypoints(Nyc, c.a(), c.b());
            ofstream prof((outdir + "profiles.asc").c_str());
            prof << "% period " << period + 1 << " averages: y  <C>  <(u+U)C>  <vC>  -vs<C>  -kappa d<C>/dy  F" << endl;
            Real Fmax = 0.0, settmax = 0.0;
            for (int ny = 0; ny < Nyc; ++ny) {
                const Real F = vCavg[ny] - vs * Cavg[ny] - kappa * dCdy[ny];
                Fmax = Greater(Fmax, abs(F));
                settmax = Greater(settmax, abs(vs * Cavg[ny]));
                prof << setprecision(15) << y[ny] << ' ' << Cavg[ny] << ' ' << uCavg[ny] << ' ' << vCavg[ny] << ' '
                     << -vs * Cavg[ny] << ' ' << -kappa * dCdy[ny] << ' ' << F << endl;
            }
            const Real Frel = settmax > 0 ? Fmax / settmax : Fmax;

            cout << "period " << period + 1 << ": dt == " << dt << " (" << N << " steps), max CFL velocity "
                 << cflu << ", scalar " << cflc << ", |u(T) - u(0)|/|u(0)| == " << du
                 << ", |c(T) - c(0)|/|c(0)| == " << dc << ", <(u+U)C> == " << setprecision(10) << uCbulk
                 << ", <vC> == " << vCbulk << ", max|F|/max|vs<C>| == " << Frel << setprecision(6) << endl;
            conv << setw(4) << period + 1 << ' ' << setprecision(10) << dt << ' ' << N << ' ' << cflu << ' ' << cflc
                 << ' ' << du << ' ' << dc << ' ' << uCbulk << ' ' << vCbulk << ' ' << Frel << endl;

            combineVelocityScalar(cfields[0], cfields[1], u4);
            u4.setPadded(flags.dealias_xz());
            u4.save(outdir + "uT");
            if (dc < tol) {
                cout << "converged: |c(T) - c(0)|/|c(0)| < " << tol << endl;
                break;
            }

            ++sinceJump;
            if (Nex >= 2 && sinceJump >= Nex && dc > exmin) {
                FlowField d1(cstart), d2(cfields[1]);
                d1 -= cstartPrev;
                d2 -= cstart;
                const Real lambda = L2IP(d2, d1) / L2IP(d1, d1);
                if (abs(lambda) < 0.9 && abs(lambda) > 1e-3) {  // either sign: c* = c + lambda/(1 - lambda) d2
                    // c + a (c - c_prev) keeps the affine BC and the mass of c (coefficients sum to 1)
                    d2 *= lambda / (1.0 - lambda);
                    cfields[1] += d2;
                    dtc = 0.0;  // restart the scalar's multistep history from the new c
                    sinceJump = 0;
                    cout << "  extrapolated: lambda == " << lambda << ", |dc|/|c| == " << L2Norm(d2) / L2Norm(cfields[1])
                         << endl;
                } else
                    cout << "  no extrapolation: lambda == " << lambda << endl;
            }

            // keep the larger CFL number in [CFLmin, CFLmax] for the next period
            const Real cfl = Greater(cflu, cflc);
            if (flags.variabledt && (cfl > flags.CFLmax || cfl < flags.CFLmin)) {
                const Real dtwant = Greater(flags.dtmin, lesser(flags.dtmax, dt * 0.5 * (flags.CFLmin + flags.CFLmax) / cfl));
                N = stepsPerPeriod(dtwant);
            }
        }
    }
    cfMPI_Finalize();
}
