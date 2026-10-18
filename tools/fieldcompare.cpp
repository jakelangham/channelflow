/**
 * Compare two congruent fields component group by component group, for regression testing.
 *
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include <iomanip>
#include <iostream>
#include "channelflow/diffops.h"
#include "channelflow/flowfield.h"
#include "channelflow/utilfuncs.h"

using namespace std;
using namespace chflow;

// Print "<label> <L2Dist(u,ref)> <L2Dist(u,ref)/L2Norm(ref)> <L2Norm(u)> <L2Norm(ref)>"
void printdist(const string& label, const FlowField& u, const FlowField& ref) {
    const Real d = L2Dist(u, ref);
    const Real nref = L2Norm(ref);
    cout << setw(9) << left << label << right << setw(24) << d << setw(24) << (nref > 0 ? d / nref : d) << setw(24)
         << L2Norm(u) << setw(24) << nref << endl;
}

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        string purpose(
            "print L2 distances between a field and a reference: one line each for the velocity\n"
            "\t(components 0-2), the scalar (component 3, if present) and all components.\n"
            "\tColumns: label, L2Dist, L2Dist/L2Norm(ref), L2Norm(u), L2Norm(ref)");
        ArgList args(argc, argv, purpose);
        const string uname = args.getstr(2, "<flowfield>", "field to test");
        const string refname = args.getstr(1, "<flowfield>", "reference field");
        args.check();

        FlowField u(uname);
        FlowField ref(refname);
        u.makeSpectral();
        ref.makeSpectral();
        if (!u.congruent(ref))
            cferror("fieldcompare: the two fields are not congruent");

        cout << setprecision(8) << scientific;
        if (u.Nd() >= 3) {
            const vector<int> vel = {0, 1, 2};
            FlowField uvel(u.Nx(), u.Ny(), u.Nz(), 3, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
            FlowField refvel(uvel);
            uvel.copySubfields(u, vel, vel);
            refvel.copySubfields(ref, vel, vel);
            printdist("velocity", uvel, refvel);
        }
        if (u.Nd() == 4)
            printdist("scalar", u[3], ref[3]);
        printdist("all", u, ref);
    }
    cfMPI_Finalize();
}
