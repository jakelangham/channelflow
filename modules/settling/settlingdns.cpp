/**
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include "modules/settling/settlingdns.h"

using namespace std;

namespace chflow {

SettlingDNS::SettlingDNS(const vector<FlowField>& fields, const DNSFlags& flags) : DNS() {
    main_nse_ = make_shared<SettlingNSE>(fields, flags);
    main_algorithm_ = newAlgorithm(fields, main_nse_, flags);
    if (!main_algorithm_->full() && flags.initstepping != flags.timestepping) {
        DNSFlags initflags = flags;
        initflags.timestepping = flags.initstepping;
        initflags.dt = flags.dt;
        init_nse_ = make_shared<SettlingNSE>(fields, initflags);
        init_algorithm_ = newAlgorithm(fields, init_nse_, initflags);
        if (init_algorithm_->Ninitsteps() != 0)
            cerr << "SettlingDNS : " << flags.initstepping << " can't initialize " << flags.timestepping
                 << " since it needs initialization itself.\n";
    }
}

}  // namespace chflow
