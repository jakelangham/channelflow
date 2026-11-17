/**
 * DNS of the settling-suspension equations (see settlingnse.h) on fields {u, c_0, ..., c_n-1, p}.
 *
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#ifndef SETTLINGDNS_H
#define SETTLINGDNS_H

#include "channelflow/dns.h"
#include "modules/settling/settlingnse.h"

namespace chflow {

class SettlingDNS : public DNS {
   public:
    // prescribedVelocity: step only the scalars (see SettlingNSE::setPrescribedVelocity)
    SettlingDNS(const std::vector<FlowField>& fields, const DNSFlags& flags, bool prescribedVelocity = false);
    virtual ~SettlingDNS() = default;

    // The velocity CFL number, or the settling one if larger: settling is explicit, and unlike v it does
    // not vanish at the walls, where the grid is finest. The settling number allows for the damping by
    // the implicit diffusion (see settlingdns.cpp).
    Real CFL(FlowField& u) const override;

    // With a prescribed velocity: provider(u) is called at the start of every step (main and start-up
    // schemes) to store the velocity of that step in u; see SettlingNSE::setVelocityProvider.
    void setVelocityProvider(const std::function<void(FlowField&)>& provider);

   private:
    std::vector<Real> vs_;
    std::vector<Real> kappa_;
};

}  // namespace chflow
#endif  // SETTLINGDNS_H
