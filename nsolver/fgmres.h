/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */

#ifndef NSOLVER_FGMRES_H
#define NSOLVER_FGMRES_H

#include "cfbasics/cfbasics.h"
using namespace Eigen;
using namespace cfbasics;

namespace nsolver {

/*==================================================================================*/
/*            Class FGMRES                                                           */
/*==================================================================================*/

class FGMRES {
    // Flexible GMRES

   public:
    FGMRES();
    FGMRES(const VectorXd& b, int Niterations, Real minCondition = 1e-13);

    const VectorXd& testVector() const;
    void iterate(const VectorXd& testvec, const VectorXd& A_testvec);

    const VectorXd& solution() const;  // current best approx to soln x of Ax=b
    const VectorXd& guess() const;
    Real residual() const;  // |Ax-b|/|b|

    int n() const;
    int Niter() const;

    MatrixXd Hn() const;  // Hn  = (n+1) x n submatrix of H
    MatrixXd Zn() const;  // Qn  = M x n     submatrix of Q
    MatrixXd AZn() const;
    MatrixXd Vn() const;        // Qn1 = M x (n+1) submatrix of Q
    void resetV();              // reset Qn to size(b)
    const MatrixXd& V() const;  // Q is large, so avoid copying it

    VectorXd solve(const VectorXd& bprime, Real& residual);
    VectorXd b();

   private:
    int M_;      // dimension of vector space
    int Niter_;  // max number of iterations
    int n_;      // current iteration number
    Real condition_;

    MatrixXd H_;   // Hessenberg Matrix
    MatrixXd V_;   // Orthogonalaized basis
    MatrixXd Z_;   // Prefined vectors
    MatrixXd AZ_;  // AZ_(:,j) = A*Z_(:,j)
    VectorXd x0_;  // guess;
    VectorXd v_;
    VectorXd qn_;
    VectorXd xn_;
    Real bnorm_;
    Real residual_;
};
}  // namespace nsolver
#endif
