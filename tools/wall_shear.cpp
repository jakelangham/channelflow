#include "cfbasics/arglist.h"
#include "channelflow/diffops.h"

using namespace std;
using namespace Eigen;
using namespace chflow;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        ArgList args(argc, argv, "output du/dy");

        const int nproc0 =
            args.getint("-np0", "--nproc0", 0, "number of MPI-processes for transpose/number of parallel ffts");
        const int nproc1 = args.getint("-np1", "--nproc1", 0, "number of MPI-processes for one fft");
        const string uname = args.getstr(1, "<flowfield>", "data file");

        args.check();

        CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);

        FlowField u(uname, cfmpi);
        FlowField dudy(u.Nx(), u.Ny(), u.Nz(), u.Nd(), u.Lx(), u.Lz(), 
                        u.a(), u.b(), u.BC(), cfmpi);
        dudy.setPadded(u.padded());

        // add on linear base state
        if (u.taskid() == u.task_coeff(0, 0))
            u.cmplx(0, 1, 0, 0) += Complex(1.0, 0.0);
        if (u.Nd() == 4 && u.taskid() == u.task_coeff(0, 0))
            u.cmplx(0, 1, 0, 3) += Complex(-1.0, 0.0);

        ydiff(u, dudy, 1);

        string outfile = uname;
        outfile.erase(outfile.end() - 3, outfile.end());
        outfile += "_diffy.nc";
        dudy.save(outfile);
    }

    cfMPI_Finalize();

    return 0;
}
