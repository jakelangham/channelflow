/**
 * Dynamical-systems interface (Newton, continuation, eigenvalues) for the settling module.
 *
 * The state is the 4-component field of the existing files (velocity + scalar in component 3). The
 * Newton vector is the velocity's independent coefficients (upstream packing) followed by a block
 * for each scalar. f^T integrates SettlingNSE on the separate fields {u, c, p}.
 *
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#ifndef SETTLINGDSI_H
#define SETTLINGDSI_H

#include "channelflow/cfdsi.h"
#include "modules/settling/settlingdns.h"

namespace chflow {

/** Newton-vector packing of one settling scalar c (stored as c + y): the coefficients left free by
 * the no-flux BC and, for the mean mode, the unit-mass constraint. */
int scalarVectorSize(const FlowField& c);
void scalarToVector(const FlowField& c, Eigen::VectorXd& x, int& n);  // writes from x(n), advances n
void vectorToScalar(const Eigen::VectorXd& x, int& n, FlowField& c, const SettlingScalar& sc);

/** Impose the scalar's no-flux BC and unit mass on c (a round trip through its packing) */
void projectScalar(FlowField& c, const SettlingScalar& sc);
/** As fixdivnoslip for a 4-component state: divergence-free no-slip velocity, and projectScalar */
void fixdivnoslipSettling(FlowField& u4, const SettlingScalar& sc);

/** Continuation parameters of the settling equations; none defers to cfDSI's (Re, T, P, ...) */
enum class settling_continuationParameter { none, Ri, vsett, kappa };

class SettlingDSI : public cfDSI {
   public:
    SettlingDSI(DNSFlags& dnsflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch,
                bool xrelative, bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u,
                std::ostream* os = &std::cout);

    void integrate(const FlowField& u, Real& T, PoincareCondition* h, FlowField& fu, const DNSFlags& flags,
                   const TimeStep& dt, int& fcount, Real& CFL, std::ostream& os) override;
    int fieldVectorSize(const FlowField& u) override;
    void fieldToVector(const FlowField& u, Eigen::VectorXd& x) override;
    void vectorToField(const Eigen::VectorXd& x, FlowField& u) override;

    /// \name Continuation in Ri, vsett, kappa (names as for -cont), or any of cfDSI's parameters.
    /// The no-flux BC follows vsett and kappa, since it is rebuilt from dnsflags_ when unpacking.
    void chooseMuSettling(std::string muName);
    void updateMu(Real mu) override;
    std::string printMu() override;
    std::string statsHeader() override;
    static settling_continuationParameter s2settling_cPar(std::string muName);
    static std::string settling_cPar2s(settling_continuationParameter cPar);

   protected:
    settling_continuationParameter settling_cPar_ = settling_continuationParameter::none;
};

}  // namespace chflow
#endif  // SETTLINGDSI_H
