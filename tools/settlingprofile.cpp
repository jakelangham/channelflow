/**
 * Write the exact steady settling profile as a 4-component field, for testing.
 *
 * With no velocity perturbation and a concentration C(y) that depends on y only, the steady
 * settling-diffusion equation kappa C'' + vs C' = 0 with no-flux walls (kappa C' + vs C = 0 at
 * y = a, b) and mean(C) = 1 has the exact solution
 *     C(y) = A exp(-r y),  r = vs/kappa,  A chosen so that mean(C) = 1,
 * which on [-1, 1] gives A = r / sinh(r). The scalar component is stored relative to the
 * background -y used by the current code, i.e. u[3] = C(y) + y.
 *
 * Optionally adds an exact decaying perturbation eps phi_m(y) cos(2 pi n z / Lz) exp(sigma t). With
 * no velocity perturbation and Ri = 0 it solves sigma phi = vs phi' + kappa (phi'' - k^2 phi) with
 * kappa phi' + vs phi = 0 at y = -1, 1 (k = 2 pi n / Lz; the base flow does not advect z-dependence).
 * Writing phi = exp(-r y / 2) psi gives psi'' = -mu^2 psi, psi' + (r/2) psi = 0 at the walls, whose
 * eigenvalues satisfy sin(mu) cos(mu) = 0:
 *   m = 0:  phi = exp(-r y),                         sigma = -kappa k^2
 *   m >= 1: mu = m pi / 2, psi = cos(mu y) + (2 mu / r) sin(mu y) (m odd)
 *                          psi = sin(mu y) - (2 mu / r) cos(mu y) (m even),
 *           sigma = -kappa (k^2 + r^2/4 + mu^2).
 *
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <iostream>
#include "channelflow/flowfield.h"
#include "channelflow/utilfuncs.h"

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose("write the exact steady settling profile C(y) = A exp(-(vs/kappa) y), mean(C) = 1,\n"
                       "\twith zero velocity perturbation (scalar stored as C + y)");
        ArgList args(argc, argv, purpose);
        const int Nx = args.getint("-Nx", "--Nx", 16, "# x gridpoints");
        const int Ny = args.getint("-Ny", "--Ny", 35, "# y gridpoints");
        const int Nz = args.getint("-Nz", "--Nz", 32, "# z gridpoints");
        const Real Lx = args.getreal("-Lx", "--Lx", 2 * pi, "streamwise (x) box length");
        const Real Lz = args.getreal("-Lz", "--Lz", pi, "spanwise (z) box length");
        const Real r = args.getreal("-vsok", "--vsoverkappa", 1.0, "vs/kappa == settling velocity / diffusivity");
        const Real eps = args.getreal("-eps", "--epsilon", 0.0, "amplitude of the decaying perturbation (0: none)");
        const int m = args.getint("-m", "--ymode", 1, "wall-normal mode number m of the perturbation");
        const int n = args.getint("-n", "--zmode", 1, "spanwise wavenumber index n of the perturbation, k = 2 pi n/Lz");
        const Real kappa = args.getreal("-kap", "--kappa", 2.5e-3, "diffusivity, sets the decay rate");
        const Real t = args.getreal("-t", "--time", 0.0, "time at which to evaluate the perturbation");
        const string outname = args.getstr(1, "<flowfield>", "output field");
        args.check();

        const Real a = -1, b = 1;
        // A = r (b - a) / (exp(-r a) - exp(-r b)), i.e. r / sinh(r) on [-1, 1]
        const Real A = (r == 0) ? 1.0 : r * (b - a) / (exp(-r * a) - exp(-r * b));

        // decaying perturbation: shape phi(y) and rate sigma
        const Real k = 2 * pi * n / Lz;
        const Real mu = m * pi / 2;
        const Real sigma = (m == 0) ? -kappa * k * k : -kappa * (k * k + r * r / 4 + mu * mu);
        auto phi = [&](Real y) {
            if (m == 0)
                return exp(-r * y);
            const Real psi = (m % 2 == 1) ? cos(mu * y) + (2 * mu / r) * sin(mu * y)
                                          : sin(mu * y) - (2 * mu / r) * cos(mu * y);
            return exp(-r * y / 2) * psi;
        };
        if (eps != 0)
            cout << "perturbation: m = " << m << ", k = " << k << ", sigma = " << sigma << ", factor exp(sigma t) = "
                 << exp(sigma * t) << endl;

        FlowField u(Nx, Ny, Nz, 4, Lx, Lz, a, b);
        u.setToZero();
        u.makePhysical();
        for (lint ny = u.nylocmin(); ny < u.nylocmax(); ++ny) {
            const Real y = u.y(ny);
            const Real c = A * exp(-r * y);
            const Real p = eps * phi(y) * exp(sigma * t);
            for (lint nx = u.nxlocmin(); nx < u.nxlocmin() + u.Nxloc(); ++nx)
                for (lint nz = 0; nz < u.Nz(); ++nz)
                    u(nx, ny, nz, 3) = c + y + p * cos(k * u.z(nz));
        }
        u.makeSpectral();
        u.save(outname);
    }
    cfMPI_Finalize();
}
