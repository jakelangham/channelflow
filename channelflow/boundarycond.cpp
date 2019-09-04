#include "channelflow/boundarycond.h"

using namespace std;
using namespace Eigen;

namespace chflow {

BoundaryCond::BoundaryCond()
    : type_(Diri), gb_(0) { 
}

BoundaryCond::BoundaryCond(int type, int Mx, int Mz, Real gb)
    : type_(type), Mx_(Mx), Mz_(Mz), gb_(gb), alpha_(0),
    gaxz_re_(Mx_ * Mz_), gaxz_im_(Mx_ * Mz_) { 
}

BoundaryCond::BoundaryCond(int type, int Mx, int Mz, Real gb, Real alpha)
    : type_(type), Mx_(Mx), Mz_(Mz), gb_(gb), alpha_(alpha),
    gaxz_re_(Mx_ * Mz_), gaxz_im_(Mx_ * Mz_) { 
}

// Set Fourier components of ga - the bottom wall boundary condition
void BoundaryCond::setgaxz(Vector& g_re, Vector& g_im, int Mx)
{
    gaxz_re_ = g_re;
    gaxz_im_ = g_im;
    Mx_ = Mx;
}

Complex BoundaryCond::gaxz(int mx, int mz) const
{
    return Complex(gaxz_re_[mx + mz * Mx_], gaxz_im_[mx + mz * Mx_]);
}

} // chflow
