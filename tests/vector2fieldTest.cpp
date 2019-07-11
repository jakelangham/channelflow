/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#include <iostream>
#include "cfbasics/cfbasics.h"
#include "channelflow/diffops.h"
#include "channelflow/flowfield.h"

using namespace std;
using namespace Eigen;
using namespace chflow;

// Test if vector2field and field2vector doesn't alter either vector or flowfield
int main(int argc, char* argv[]) {
    int failure = 0;
    cfMPI_Init(&argc, &argv);
    {
        CfMPI* cfmpi = &CfMPI::getInstance();
        FlowField u("data/eq1strat", cfmpi);
        FlowField rho(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), cfmpi);
        VectorXd v, v2, v3;
        FlowField u2(rho);
        FlowField u3(rho);

        vector<int> rho_index = {3};
        vector<int> rho_pslim_index = {0};
        rho.copySubfields(u, rho_index, rho_pslim_index);

        // Check if a conversion of the original field to a vector and back doesn't diverge too much.
        field2vector(rho, v);
        cout << "L2Norm(v)        = " << L2Norm(v) << endl;
        vector2field(v, u2);
        Real err = 0;
        cout << "L2Dist(rho, u2)    = " << L2Dist(rho, u2) << endl;

        err += L2Dist(rho, u2);

        // Check if a conversion to another vector gives the same
        field2vector(u2, v2);
        v3 = v2;
        v2 -= v;
        cout << "L2Dist(v, v2)    = " << L2Norm(v2) << endl;
        cout << "L2Norm(v2)       = " << L2Norm(v3) << endl;

        err += L2Norm(v2);

        // Check if a reconversion to a field gives the same
        vector2field(v3, u3);
        cout << "L2Dist(u2, u3)    = " << L2Dist(u3, u2) << endl;
        err += L2Dist(u3, u2);
        cout << "err = " << err << endl;

        if (err > 2e-16) {
            cerr << "\t** FAIL **" << endl;
            cout << "\t** FAIL **" << endl;
            failure = 1;
        } else {
            cerr << "\t   pass   " << endl;
            cout << "\t   pass   " << endl;
        }
    }
    cfMPI_Finalize();
    return failure;
}
