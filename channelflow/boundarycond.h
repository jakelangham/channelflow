/**
 * Class to store boundary condition data
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 */

#ifndef CHANNELFLOW_BOUNDARYCOND_H
#define CHANNELFLOW_BOUNDARYCOND_H

#include "cfbasics/mathdefs.h"

namespace chflow {

enum BC { Free, Diri, Zero, Neum, Mixed };

class BoundaryCond {
   public:
    BoundaryCond();
    BoundaryCond(int type, Real ga, Real gb);
    BoundaryCond(int type, Real ga, Real gb, Real alpha);

    int type_;
    Real ga_;
    Real gb_;
    Real alpha_;
};

} // namespace chflow

#endif
