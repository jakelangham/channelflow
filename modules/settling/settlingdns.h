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
    SettlingDNS(const std::vector<FlowField>& fields, const DNSFlags& flags);
    virtual ~SettlingDNS() = default;
};

}  // namespace chflow
#endif  // SETTLINGDNS_H
