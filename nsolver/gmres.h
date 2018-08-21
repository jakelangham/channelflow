/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */

#ifndef NSOLVER_GMRES_H
#define NSOLVER_GMRES_H

#include "cfbasics/cfbasics.h"
using namespace Eigen;
using namespace cfbasics;

namespace nsolver {

/*==================================================================================*/
/*            Class GMRES                                                           */
/*==================================================================================*/

// Iterative GMRES solution to Ax=b.
// Usage for classic iterative GMRES solution of Ax=b:
// GMRES gmres(b, N);
// for (int n=0; n<N; ++n) {
//   VectorXd q  = gmres.testVector();
//   VectorXd Aq = A*q; // or however else you calculate A*q
//   gmres.iterate(Aq);
//   VectorXd x = gmres.solution(); // current estimate of soln
//   cout << "krylov residual == " << gmres.residual() << endl;
// }
//
// Additional functionality:
// Find approx solution x' of Ax'=b' projected into current Krylov subspace.
// Real residual;
// Vector xprime = grmes.solve(bprime, residual);

class GMRES {
   public:
    GMRES();
    GMRES(const VectorXd& b, int Niterations, Real minCondition = 1e-13);

    const VectorXd& testVector() const;
    void iterate(const VectorXd& A_testvec);

    const VectorXd& solution() const;  // current best approx to soln x of Ax=b
    const VectorXd& guess() const;
    Real residual() const;  // |Ax-b|/|b|

    int n() const;
    int Niter() const;

    MatrixXd Hn() const;        // Hn  = (n+1) x n submatrix of H
    MatrixXd Qn() const;        // Qn  = M x n     submatrix of Q
    MatrixXd Qn1() const;       // Qn1 = M x (n+1) submatrix of Q
    void resetQ();              // reset Qn to size(b)
    const MatrixXd& Q() const;  // Q is large, so avoid copying it

    VectorXd solve(const VectorXd& bprime, Real& residual);

   private:
    int M_;      // dimension of vector space
    int Niter_;  // max number of iterations
    int n_;      // current iteration number
    Real condition_;

    MatrixXd H_;
    MatrixXd Q_;
    MatrixXd Vn_;
    VectorXd x0_;  // guess;
    VectorXd v_;
    VectorXd qn_;
    VectorXd xn_;
    Real bnorm_;
    Real residual_;
};
}  // namespace nsolver
#endif
