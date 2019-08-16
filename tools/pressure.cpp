/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "channelflow/flowfield.h"
#include "channelflow/poissonsolver.h"
#include "channelflow/utilfuncs.h"

using namespace std;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose("compute pressure field of a given velocity field");

        ArgList args(argc, argv, purpose);
        const string nonlstr = args.getstr("-nl", "--nonlinearity", "rot",
                                           "method of calculating "
                                           "nonlinearity, one of [rot conv div skew alt]");
        const string uname = args.getstr(2, "<flowfield>", "input velocity field (deviation from laminar)");
        const string pname =
            args.getstr(1, "<flowfield>", "output pressure field (not including const pressure gradient)");

        string Uname, Wname;
        DNSFlags baseflags = setBaseFlowFlags(args, Uname, Wname);
        baseflags.nonlinearity = s2nonlmethod(nonlstr);
        args.check();

        // define all input for pressure
        FlowField u(uname);
        vector<ChebyCoeff> base_Flow = baseFlow(u.Ny(), u.a(), u.b(), baseflags, Uname, Wname);

        FlowField uvel(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
        vector<int> vel_indices = {0, 1, 2};
        uvel.copySubfields(u, vel_indices, vel_indices);

        // compute pressure
        PressureSolver poisson(uvel, base_Flow[0], base_Flow[1], baseflags.nu, baseflags.Vsuck, baseflags.nonlinearity);

        FlowField q = poisson.solve(uvel);

        q.setPadded(true);
        q.save(pname);
    }
    cfMPI_Finalize();
}
