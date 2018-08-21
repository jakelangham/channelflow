/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#ifndef NSOLVER_LANCZOS_H
#define NSOLVER_LANCZOS_H
#include "arnoldi.h"
#include "cfbasics/cfbasics.h"

using namespace Eigen;
using namespace cfbasics;
namespace nsolver {

/*==================================================================================*/
/*            Class Lanczos                                                         */
/*==================================================================================*/

// Lanczos iteration to estimate eigenvalues of Hermitian matrix A. Usage:
// Lanczos lanczos(N, b);
// for (int n=0; n<N; ++n) {
//   VectorXd q  = lanczos.testVector();
//   VectorXd Aq = A*q; // or however else you calculate A*q
//   lanczos.iterate(Aq);
//   VectorXd ew = lanczos.ew(); // current estimate of eigenvalues
// }

class Lanczos : public Arnoldi {
   public:
    Lanczos();
    Lanczos(const VectorXd& b, int Niterations, Real minCondition = 1e-13);

    void iterate(const VectorXd& Aq) override;  // tell Lanczos the value of Aq

   private:
    MatrixXd T_;  // Tridiagonal matrix

    void eigencalc() override;
};

}  // namespace nsolver
#endif
