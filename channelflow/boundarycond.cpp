#include "channelflow/boundarycond.h"

using namespace std;

namespace chflow {

BoundaryCond::BoundaryCond()
    : type_(Diri), ga_(0), gb_(0) { ; }

BoundaryCond::BoundaryCond(int type, Real ga, Real gb)
    : type_(type), ga_(ga), gb_(gb), alpha_(0) { ; }

BoundaryCond::BoundaryCond(int type, Real ga, Real gb, Real alpha)
    : type_(type), ga_(ga), gb_(gb), alpha_(alpha) { ; }

} // chflow
