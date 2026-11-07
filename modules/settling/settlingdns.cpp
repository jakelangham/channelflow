/**
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include "modules/settling/settlingdns.h"

using namespace std;

namespace chflow {

SettlingDNS::SettlingDNS(const vector<FlowField>& fields, const DNSFlags& flags, bool prescribedVelocity) : DNS() {
    for (const SettlingScalar& s : settlingScalars(flags)) {
        vs_.push_back(s.vs);
        kappa_.push_back(s.kappa);
    }
    if (prescribedVelocity) {
        // MultistepDNS only: Runge-Kutta and CNAB-style substeps would need the velocity between steps
        auto multistep = [](TimeStepMethod m) {
            return m == SBDF1 || m == SBDF2 || m == SBDF3 || m == SBDF4 || m == CNFE1;
        };
        if (!multistep(flags.timestepping) || !multistep(flags.initstepping))
            cferror("SettlingDNS : a prescribed velocity needs multistep schemes (e.g. -ts sbdf3 -is sbdf1)");
    }
    auto nse = make_shared<SettlingNSE>(fields, flags);
    nse->setPrescribedVelocity(prescribedVelocity);
    main_nse_ = nse;
    main_algorithm_ = newAlgorithm(fields, main_nse_, flags);
    if (!main_algorithm_->full() && flags.initstepping != flags.timestepping) {
        DNSFlags initflags = flags;
        initflags.timestepping = flags.initstepping;
        initflags.dt = flags.dt;
        auto initnse = make_shared<SettlingNSE>(fields, initflags);
        initnse->setPrescribedVelocity(prescribedVelocity);
        init_nse_ = initnse;
        init_algorithm_ = newAlgorithm(fields, init_nse_, initflags);
        if (init_algorithm_->Ninitsteps() != 0)
            cerr << "SettlingDNS : " << flags.initstepping << " can't initialize " << flags.timestepping
                 << " since it needs initialization itself.\n";
    }
}

Real SettlingDNS::CFL(FlowField& u) const {
    Real cfl = DNS::CFL(u);

    // Explicit settling at speed vs in y (the advection by v is the velocity's, covered above), damped
    // by the implicit diffusion: at wavenumber k the explicit growth vs dt k is offset by kappa dt k^2,
    // so take vs dt k/(1 + kappa dt k^2) at k = pi/dy (dy as in FlowField::CFLfactor). It is vs dt k on
    // coarse grids and at most (vs/2) sqrt(dt/kappa) on fine ones, independent of Ny. SBDF3 is unstable
    // at about 0.37 (laminar Couette flow, vs/kappa = 160, Ny = 181 and 241: stable at dt = 0.008,
    // unstable at 0.009), so 1 marks the limit. The same for every process: no MPI reduction needed.
    const Real limit = 0.37;
    const Real h = dt();
    Vector y = chebypoints(u.Ny(), u.a(), u.b());
    Real cflc = 0.0;
    for (int ny = 0; ny < u.Ny(); ++ny) {
        Real dy = (ny == 0 || ny == u.Ny() - 1) ? y[0] - y[1] : (y[ny - 1] - y[ny + 1]) / 2.0;
        const Real k = pi / dy;
        for (uint i = 0; i < vs_.size(); ++i)
            cflc = Greater(cflc, abs(vs_[i]) * h * k / (1.0 + kappa_[i] * h * k * k));
    }
    return Greater(cfl, cflc / limit);
}

}  // namespace chflow
