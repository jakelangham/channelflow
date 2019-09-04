/**
 * Class to store boundary condition data
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 */

#ifndef CHANNELFLOW_BOUNDARYCOND_H
#define CHANNELFLOW_BOUNDARYCOND_H

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"

using namespace Eigen;

namespace chflow {

enum BC { Free, Diri, Zero, Neum, Mixed };

class BoundaryCond {
   public:
    BoundaryCond();
    BoundaryCond(int type, int Mx, int Mz, Real gb);
    BoundaryCond(int type, int Mx, int Mz, Real gb, Real alpha);
    void setgaxz(Vector& g_re, Vector& g_im, int Mx);
    Complex gaxz(int mx, int mz) const;

    int type_;
    int Mx_;
    int Mz_;
    Real gb_;
    Real alpha_;

    // these are required in the case of nonconstant ga
    Vector gaxz_re_;
    Vector gaxz_im_;
};

} // namespace chflow

#endif
