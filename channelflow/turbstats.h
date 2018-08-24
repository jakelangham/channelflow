/**
 * Class for computing basic statistics of turbulent flows
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author:
 */

#ifndef CHANNELFLOW_TURBSTATS_H
#define CHANNELFLOW_TURBSTATS_H

#include "cfbasics/mathdefs.h"
#include "channelflow/basisfunc.h"
#include "channelflow/chebyshev.h"
#include "channelflow/flowfield.h"

namespace channelflow {

class TurbStats {
   public:
    TurbStats();
    TurbStats(const std::string& filebase);
    TurbStats(const ChebyCoeff& Ubase, cfbasics::Real nu);

    void reset();
    void addData(FlowField& u, FlowField& tmp);
    void msave(const std::string& filebase, bool wallunits = false) const;

    // Terminology: utot = U + u = Ubase + ubase
    ChebyCoeff U() const;  // average of total flow
    ChebyCoeff dUdy() const;
    ChebyCoeff Ubase() const;  // base flow (e.g. parabola)
    ChebyCoeff ubase() const;  // mean fluctuation on base flow
    ChebyCoeff uu() const;
    ChebyCoeff uv() const;
    ChebyCoeff uw() const;
    ChebyCoeff vv() const;
    ChebyCoeff vw() const;
    ChebyCoeff ww() const;

    // wall unit stuff
    cfbasics::Real ustar() const;               // sqrt(nu <d/dy utot>)
    cfbasics::Real parabolicReynolds() const;   // h Uparab/nu, center vel of parab w = flux
    cfbasics::Real bulkReynolds() const;        // h Ubulk/nu
    cfbasics::Real centerlineReynolds() const;  // h Ucenterline/nu
    cfbasics::Real hplus() const;               // ustar/nu (b-a)/2
    cfbasics::Vector yplus() const;             // ustar/nu (y-a)

   private:
    int count_;
    cfbasics::Real nu_;  // All ChebyCoeff quantities are sums for means, in utot
    ChebyCoeff Ubase_;   // base flow (parabolas, etc).
    ChebyCoeff ubase_;   // mean fluc above base flow: utot = Ubase + ubase
    ChebyCoeff U_;       // mean flow = avg(utot)
    ChebyCoeff uu_;      // sum_1^count utot utot
    ChebyCoeff uv_;
    ChebyCoeff uw_;
    ChebyCoeff vv_;
    ChebyCoeff vw_;
    ChebyCoeff ww_;
};

}  // namespace channelflow

#endif  // TURBSTATS
