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
#include "modules/settling/settlingstats.h"
#include "nsolver/gmres.h"

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

// g = f on g's grid, for spectral fields that differ at most in Ny (f.Ny() <= g.Ny()) with the same x-z
// distribution: each process copies its own modes' Chebyshev coefficients and zeroes the rest. Same result
// as FlowField::interpolate (including its zeroed x-Nyquist mode), which loops over every global mode on
// every process (and copies f) and so does not scale with the number of processes; called every step.
// Falls back to it otherwise.
static void interpolateY(const FlowField& f, FlowField& g) {
    const bool local = f.xzstate() == Spectral && f.ystate() == Spectral && f.Nx() == g.Nx() &&
                       f.Nz() == g.Nz() && f.Ny() <= g.Ny() && f.Nd() == g.Nd() && f.Lx() == g.Lx() &&
                       f.Lz() == g.Lz() && f.a() == g.a() && f.b() == g.b() && f.cfmpi() == g.cfmpi() &&
                       f.mxlocmin() == g.mxlocmin() && f.Mxloc() == g.Mxloc() && f.mzlocmin() == g.mzlocmin() &&
                       f.Mzloc() == g.Mzloc();
    if (!local) {
        g.interpolate(f);
        return;
    }
    g.setState(Spectral, Spectral);
    const int fNy = f.Ny(), gNy = g.Ny();
    const int mxNyq = f.Nx() / 2;  // interpolate's kx loop stops before kxmax() = Nx/2, odd Nx included
    for (int i = 0; i < f.Nd(); ++i)
        for (int mx = f.mxlocmin(); mx < f.mxlocmin() + f.Mxloc(); ++mx)
            for (int mz = f.mzlocmin(); mz < f.mzlocmin() + f.Mzloc(); ++mz) {
                const bool nyq = (mx == mxNyq);
                for (int ny = 0; ny < fNy; ++ny)
                    g.cmplx(mx, ny, mz, i) = nyq ? Complex(0.0, 0.0) : f.cmplx(mx, ny, mz, i);
                for (int ny = fNy; ny < gNy; ++ny)
                    g.cmplx(mx, ny, mz, i) = Complex(0.0, 0.0);
            }
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
                                      "save this many states per period, equally spaced over T, with energy.asc and "
                                      "cstats.asc as simulateflow writes them (overwritten each period; 0 == none)");
        const int Navg = args.getint("-na", "--Naverage", 400,
                                     "samples per period for the period-averaged flux profiles (trapezoid rule)");
        const int Nex = args.getint("-ex", "--extrapolate", 0,
                                    "every this many periods, extrapolate c past its dominant decaying mode "
                                    "(Aitken; >= 2, 0 == never)");
        const Real exmin = args.getreal("-exmin", "--extrapolateMin", 0.0,
                                        "extrapolate only while |c(T) - c(0)|/|c(0)| > exmin");
        const int Ng = args.getint("-gm", "--gmres", 0,
                                   "solve for the periodic c with GMRES (Newton steps of at most this many Krylov "
                                   "periods each) instead of iterating; 0 == iterate");
        string chdir = args.getstr("-ch", "--history", "",
                                   "directory with uTm1.nc, uTm2.nc of a previous run (same grid and dt): "
                                   "continue without restarting the multistep history");
        if (chdir != "" && chdir.back() != '/')
            chdir += '/';
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
        if (scalars[0].stratified())
            cferror("scalarorbit : its profiles and fluxes assume a settling scalar (-strat none)");

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
        ofstream conv((outdir + "convergence.asc").c_str());
        conv << "% period  dt  steps  maxCFL(velocity)  maxCFL(scalar)  L2Dist(u(T),u(0))/L2Norm(u(0))  "
                "L2Dist(c(T),c(0))/L2Norm(c(0))  <(u+U)C>  <vC>  max|F|/max|vs<C>|"
             << endl;
        const Real vs = scalars[0].vs;
        const Real kappa = scalars[0].kappa;
        FlowField ftmp;
        ChebyCoeff C, uC, vC;

        // The state of the scalar's SBDF scheme at the start of a period: c at t = 0, -dt and -2dt (and the
        // velocity at -dt, -2dt, i.e. at T - dt, T - 2dt on the orbit). Pushing the earlier levels into the
        // multistep history (DNS::push) continues the scheme exactly, without a lower-order restart; without
        // a history (valid == false) the scheme restarts with -is, which costs ~1e-6.
        struct State {
            FlowField c0, c1, c2;  // c at t = 0, -dt, -2dt
            bool valid = false;
        };
        FlowField uh1, uh2;  // velocity on c's grid at T - dt, T - 2dt (recorded each period)
        bool uhValid = false;

        State X;
        X.c0 = c;
        // Placeholders until the first period records the history (unused while !valid). States are copied,
        // and under MPI copying an empty FlowField fails (FFTW-MPI's layout with zero block sizes).
        X.c1 = c;
        X.c2 = c;
        if (chdir != "") {
            FlowField h1(chdir + "uTm1"), h2(chdir + "uTm2");
            if (!h1.geomCongruent(c) || !h2.geomCongruent(c))
                cferror("scalarorbit : the history in " + chdir + " is not on the scalar's grid");
            splitVelocityScalar(h1, uh1, X.c1);
            splitVelocityScalar(h2, uh2, X.c2);
            X.valid = uhValid = true;
            cout << "continuing the multistep history from " << chdir << " (it must have the same dt)" << endl;
        }

        int periodCount = 0;
        Real lastDu = 0.0, lastCfl = 0.0;

        // One period from X; returns the state at T. With stats, also writes the period averages, the
        // snapshots (-ns) and a line of convergence.asc.
        auto period = [&](const State& Xin, bool stats) {
            const Real dt = T / N;
            uflags.dt = dt;
            FlowField u = u0;
            FlowField q(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), cfmpi);
            DNS udns({u, q}, uflags);
            PressureSolver psolver(u, udns.Ubase(), udns.Wbase(), uflags.nu, uflags.Vsuck, uflags.nonlinearity);
            psolver.solve(q, u);
            vector<FlowField> ufields = {u, q};

            DNSFlags cflags = flags;
            cflags.dt = dt;
            uc.interpolate(u0);
            vector<FlowField> cfields = {uc, Xin.c0, qc};
            SettlingDNS cdns(cfields, cflags, true);
            if (Xin.valid && uhValid) {
                cdns.push({uh2, Xin.c2, qc});
                cdns.push({uh1, Xin.c1, qc});
            }

            State Y;
            Real cflu = 0.0, cflc = 0.0;
            const int checkint = Greater(1, int(round(1.0 / dt)));  // CFL check about once per time unit
            ChebyCoeff Cavg(Nyc, c.a(), c.b(), Spectral), uCavg(Cavg), vCavg(Cavg);
            ofstream eout, cout_;
            FlowField svel, suc;
            if (stats && Nsave > 0) {  // as simulateflow's energy.asc and cstats.asc, at the saved states
                eout.open((outdir + "energy.asc").c_str());
                eout << fieldstatsheader_t() << endl;
                cout_.open((outdir + "cstats.asc").c_str());
                cout_ << concstatsheader_t() << endl;
                cdns.Ubase().save(outdir + "Ubase");
                cdns.Wbase().save(outdir + "Wbase");
            }
            // The velocity of each step comes from udns through the provider, so the scalar is advanced
            // from one step that needs the fields (flux sample, snapshot, CFL check, the last two steps)
            // to the next in one call: DNS::advance copies all fields in and out per call.
            cdns.setVelocityProvider([&](FlowField& ustep) {
                interpolateY(ufields[0], ustep);  // u(t_n) for the step from t_n
                udns.advance(ufields, 1);         // then the velocity to t_{n+1}
            });
            auto event = [&](int n) {
                return n >= N - 2 || (stats && n % (N / Navg) == 0) || n % checkint == 0 ||
                       (stats && Nsave > 0 && n % (N / Nsave) == 0);
            };
            for (int n = 0; n < N;) {
                interpolateY(ufields[0], cfields[0]);
                if (n == N - 2) {
                    uh2 = cfields[0];
                    Y.c2 = cfields[1];
                } else if (n == N - 1) {
                    uh1 = cfields[0];
                    Y.c1 = cfields[1];
                }
                if (stats && n % (N / Navg) == 0) {
                    fluxProfiles(cfields[0], cfields[1], cdns.Ubase(), ftmp, C, uC, vC);
                    Cavg += C;
                    uCavg += uC;
                    vCavg += vC;
                }
                if (n % checkint == 0) {
                    cflu = Greater(cflu, udns.CFL(ufields[0]));
                    cflc = Greater(cflc, cdns.CFL(cfields[0]));
                }
                if (stats && Nsave > 0 && n % (N / Nsave) == 0) {
                    combineVelocityScalar(cfields[0], cfields[1], u4);
                    u4.setPadded(flags.dealias_xz());
                    u4.save(outdir + "u" + i2s(n / (N / Nsave)));
                    eout << fieldstats_t(u4, n * dt) << endl;
                    cout_ << concstats_t(u4, svel, suc, n * dt) << endl;
                }
                // While the multistep history fills, DNS::advance pushes the given fields (evaluating their
                // nonlinear term) before each start-up step, so cfields[0] must hold u(t_n): one step at a time.
                int next = n + 1;
                if (cdns.full())
                    while (next < N && !event(next))
                        ++next;
                cdns.advance(cfields, next - n);  // advances udns too, step by step
                n = next;
            }
            interpolateY(ufields[0], cfields[0]);  // u(T), for the saved state
            Y.c0 = cfields[1];
            Y.valid = uhValid = true;
            lastDu = L2Dist(ufields[0], u0) / L2Norm(u0);
            lastCfl = Greater(cflu, cflc);
            ++periodCount;

            const Real dc = L2Dist(Y.c0, Xin.c0) / L2Norm(Xin.c0);
            if (!stats) {
                cout << "period " << periodCount << " (Krylov): |c(T) - c(0)|/|c(0)| == " << dc << endl;
                return Y;
            }

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
            prof << "% period " << periodCount << " averages: y  <C>  <(u+U)C>  <vC>  -vs<C>  -kappa d<C>/dy  F"
                 << endl;
            Real Fmax = 0.0, settmax = 0.0;
            for (int ny = 0; ny < Nyc; ++ny) {
                const Real F = vCavg[ny] - vs * Cavg[ny] - kappa * dCdy[ny];
                Fmax = Greater(Fmax, abs(F));
                settmax = Greater(settmax, abs(vs * Cavg[ny]));
                prof << setprecision(15) << y[ny] << ' ' << Cavg[ny] << ' ' << uCavg[ny] << ' ' << vCavg[ny] << ' '
                     << -vs * Cavg[ny] << ' ' << -kappa * dCdy[ny] << ' ' << F << endl;
            }
            const Real Frel = settmax > 0 ? Fmax / settmax : Fmax;
            cout << "period " << periodCount << ": dt == " << dt << " (" << N << " steps), max CFL velocity "
                 << cflu << ", scalar " << cflc << ", |u(T) - u(0)|/|u(0)| == " << lastDu
                 << ", |c(T) - c(0)|/|c(0)| == " << dc << ", <(u+U)C> == " << setprecision(10) << uCbulk
                 << ", <vC> == " << vCbulk << ", max|F|/max|vs<C>| == " << Frel << setprecision(6)
                 << (Xin.valid ? "" : " (restarted)") << endl;
            conv << setw(4) << periodCount << ' ' << setprecision(10) << dt << ' ' << N << ' ' << cflu << ' '
                 << cflc << ' ' << lastDu << ' ' << dc << ' ' << uCbulk << ' ' << vCbulk << ' ' << Frel << endl;

            // the state at T, and the history levels for continuing (-ch)
            combineVelocityScalar(cfields[0], Y.c0, u4);
            u4.setPadded(flags.dealias_xz());
            u4.save(outdir + "uT");
            combineVelocityScalar(uh1, Y.c1, u4);
            u4.setPadded(flags.dealias_xz());
            u4.save(outdir + "uTm1");
            combineVelocityScalar(uh2, Y.c2, u4);
            u4.setPadded(flags.dealias_xz());
            u4.save(outdir + "uTm2");
            return Y;
        };

        auto relChange = [](const State& A, const State& B) { return L2Dist(A.c0, B.c0) / L2Norm(B.c0); };

        if (Ng > 0) {
            // GMRES: Phi(X) = P X + b is affine, so the periodic state solves (P - I) d = -(Phi(X) - X), with
            // P v = (Phi(X + s v) - Phi(X))/s (exact up to round-off). The unknown is the whole SBDF state
            // (c at 0, -dt, -2dt), so the fixed point is the periodic solution of the scheme itself.
            auto pack = [](const State& S, Eigen::VectorXd& v) {
                const FlowField* f[3] = {&S.c0, &S.c1, &S.c2};
                const lint n1 = 2 * f[0]->Mxloc() * f[0]->My() * f[0]->Mzloc();
                v.resize(3 * n1);
                lint k = 0;
                for (int l = 0; l < 3; ++l)
                    for (int mx = f[l]->mxlocmin(); mx < f[l]->mxlocmin() + f[l]->Mxloc(); ++mx)
                        for (int my = 0; my < f[l]->My(); ++my)
                            for (int mz = f[l]->mzlocmin(); mz < f[l]->mzlocmin() + f[l]->Mzloc(); ++mz) {
                                const Complex z = f[l]->cmplx(mx, my, mz, 0);
                                v(k++) = Re(z);
                                v(k++) = Im(z);
                            }
            };
            auto unpack = [](const Eigen::VectorXd& v, State& S) {
                FlowField* f[3] = {&S.c0, &S.c1, &S.c2};
                lint k = 0;
                for (int l = 0; l < 3; ++l)
                    for (int mx = f[l]->mxlocmin(); mx < f[l]->mxlocmin() + f[l]->Mxloc(); ++mx)
                        for (int my = 0; my < f[l]->My(); ++my)
                            for (int mz = f[l]->mzlocmin(); mz < f[l]->mzlocmin() + f[l]->Mzloc(); ++mz) {
                                f[l]->cmplx(mx, my, mz, 0) = Complex(v(k), v(k + 1));
                                k += 2;
                            }
            };
            if (!X.valid) {
                cout << "GMRES: one period to set up the multistep history" << endl;
                X = period(X, true);
            }
            for (int newton = 0; newton < Nperiods; ++newton) {
                State Y = period(X, true);
                const Real dc = relChange(Y, X);
                if (dc < tol) {
                    cout << "converged: |c(T) - c(0)|/|c(0)| == " << dc << " < " << tol << endl;
                    break;
                }
                Eigen::VectorXd x, yv;
                pack(X, x);
                pack(Y, yv);
                const Real s = L2Norm(x);  // over all processes (x holds this process's modes only)
                Eigen::VectorXd b = x - yv;  // -(Phi(X) - X)
                GMRES gmres(b, Ng, 1e-14);
                for (int it = 0; it < Ng; ++it) {
                    State Xq = X;
                    unpack(x + s * gmres.testVector(), Xq);
                    Xq.valid = true;
                    State Yq = period(Xq, false);
                    Eigen::VectorXd yq;
                    pack(Yq, yq);
                    gmres.iterate((yq - yv) / s - gmres.testVector());
                    cout << "  GMRES " << it + 1 << ": residual " << gmres.residual() << endl;
                    if (gmres.residual() * dc < 0.1 * tol)
                        break;
                }
                unpack(x + gmres.solution(), X);
                X.valid = true;
            }
        } else {
            // Iterate c -> Phi(c), carrying the multistep history from period to period; optionally
            // extrapolate past the dominant decaying mode (Aitken), on all history levels at once: one
            // period maps the state to P X + b, so after a few periods the error is ~ lambda^n v; from the
            // last two changes d1, d2 of c(0), lambda ~ <d2, d1>/<d1, d1>, X* ~ Y + lambda/(1 - lambda) (Y - X).
            State Xprev;
            int sinceJump = 0;
            for (int p = 0; p < Nperiods; ++p) {
                const int Nold = N;
                State Y = period(X, true);
                const Real dc = relChange(Y, X);
                if (dc < tol) {
                    cout << "converged: |c(T) - c(0)|/|c(0)| < " << tol << endl;
                    break;
                }
                ++sinceJump;
                State Xnext = Y;
                if (Nex >= 2 && sinceJump >= Nex && dc > exmin && Xprev.valid) {
                    FlowField d1(X.c0), d2(Y.c0);
                    d1 -= Xprev.c0;
                    d2 -= X.c0;
                    const Real lambda = L2IP(d2, d1) / L2IP(d1, d1);
                    if (abs(lambda) < 0.9 && abs(lambda) > 1e-3) {
                        // Y + a (Y - X) keeps the affine BC and the mass (coefficients sum to 1)
                        const Real a = lambda / (1.0 - lambda);
                        FlowField* yl[3] = {&Xnext.c0, &Xnext.c1, &Xnext.c2};
                        const FlowField* xl[3] = {&X.c0, &X.c1, &X.c2};
                        for (int l = 0; l < 3; ++l) {
                            FlowField d(*yl[l]);
                            d -= *xl[l];
                            d *= a;
                            *yl[l] += d;
                        }
                        sinceJump = 0;
                        cout << "  extrapolated: lambda == " << lambda << endl;
                    } else
                        cout << "  no extrapolation: lambda == " << lambda << endl;
                }
                Xprev = sinceJump == 0 ? State() : X;
                X = Xnext;

                // keep the larger CFL number in [CFLmin, CFLmax] for the next period (a new dt restarts)
                if (flags.variabledt && (lastCfl > flags.CFLmax || lastCfl < flags.CFLmin)) {
                    const Real dtwant = Greater(
                        flags.dtmin, lesser(flags.dtmax, (T / N) * 0.5 * (flags.CFLmin + flags.CFLmax) / lastCfl));
                    N = stepsPerPeriod(dtwant);
                }
                if (N != Nold) {  // the history belongs to the old dt
                    X.valid = false;
                    Xprev = State();
                }
            }
        }
    }
    cfMPI_Finalize();
}
