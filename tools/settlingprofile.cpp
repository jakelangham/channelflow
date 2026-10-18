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
        const string outname = args.getstr(1, "<flowfield>", "output field");
        args.check();

        const Real a = -1, b = 1;
        // A = r (b - a) / (exp(-r a) - exp(-r b)), i.e. r / sinh(r) on [-1, 1]
        const Real A = (r == 0) ? 1.0 : r * (b - a) / (exp(-r * a) - exp(-r * b));

        FlowField u(Nx, Ny, Nz, 4, Lx, Lz, a, b);
        u.setToZero();
        u.makePhysical();
        for (lint ny = u.nylocmin(); ny < u.nylocmax(); ++ny) {
            const Real y = u.y(ny);
            const Real c = A * exp(-r * y);
            for (lint nx = u.nxlocmin(); nx < u.nxlocmin() + u.Nxloc(); ++nx)
                for (lint nz = 0; nz < u.Nz(); ++nz)
                    u(nx, ny, nz, 3) = c + y;
        }
        u.makeSpectral();
        u.save(outname);
    }
    cfMPI_Finalize();
}
