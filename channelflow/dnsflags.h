/**
 * Control parameters for time-integration by spectral Navier-Stokes simulator
 * DNSFlags specifies all relevant parameters for integrating Navier-Stokes equation in standard channel domains.
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original author: John F. Gibson
 */

#ifndef CHANNELFLOW_DNSFLAGS_H
#define CHANNELFLOW_DNSFLAGS_H

#include "cfbasics/arglist.h"
#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"
#include "channelflow/flowfield.h"

namespace channelflow {

class DNSFlags;

// Enum types for specifying the behavior of DNS, fields of DNSFlags.

enum VelocityScale { WallScale, ParabolicScale };  // BulkScale
enum BaseFlow { ZeroBase, LinearBase, ParabolicBase, LaminarBase, SuctionBase, ArbitraryBase };
enum MeanConstraint { PressureGradient, BulkVelocity };
enum TimeStepMethod { CNFE1, CNAB2, CNRK2, SMRK2, SBDF1, SBDF2, SBDF3, SBDF4 };
enum NonlinearMethod {
    Rotational,
    Convection,
    Divergence,
    SkewSymmetric,
    Alternating,
    Alternating_,
    LinearAboutProfile
};
enum Dealiasing { NoDealiasing, DealiasXZ, DealiasY, DealiasXYZ };
enum Verbosity { Silent, PrintTime, PrintTicks, VerifyTauSolve, PrintAll };

// Urgh, boilerplate code
VelocityScale s2velocityscale(const std::string& s);
BaseFlow s2baseflow(const std::string& s);
MeanConstraint s2constraint(const std::string& s);
TimeStepMethod s2stepmethod(const std::string& s);
NonlinearMethod s2nonlmethod(const std::string& s);
Dealiasing s2dealiasing(const std::string& s);
Verbosity s2verbosity(const std::string& s);

int getBodyforcefromLine(int taskid, std::ifstream& is);

std::string velocityscale2string(VelocityScale vs);
std::string baseflow2string(BaseFlow bf);
std::string constraint2string(MeanConstraint mc);
std::string stepmethod2string(TimeStepMethod ts);
std::string nonlmethod2string(NonlinearMethod nm);
std::string dealiasing2string(Dealiasing d);
std::string verbosity2string(Verbosity v);

std::ostream& operator<<(std::ostream& os, VelocityScale v);
std::ostream& operator<<(std::ostream& os, BaseFlow b);
std::ostream& operator<<(std::ostream& os, MeanConstraint m);
std::ostream& operator<<(std::ostream& os, TimeStepMethod t);
std::ostream& operator<<(std::ostream& os, NonlinearMethod n);
std::ostream& operator<<(std::ostream& os, Dealiasing d);
std::ostream& operator<<(std::ostream& os, Verbosity v);

class BodyForce {
   public:
    BodyForce();

    /** \brief The infamous virtual destructor */
    virtual ~BodyForce() = default;

    cfbasics::Vector operator()(cfbasics::Real x, cfbasics::Real y, cfbasics::Real z, cfbasics::Real t);
    void eval(cfbasics::Real t, FlowField& f);
    virtual void eval(cfbasics::Real x, cfbasics::Real y, cfbasics::Real z, cfbasics::Real t, cfbasics::Real& fx, cfbasics::Real& fy, cfbasics::Real& fz);
    virtual bool isOn(cfbasics::Real t);
};

// Specify the behavior of NSIntegrators by setting fields of DNSFlags.
class DNSFlags {
   public:
    //        type name       default
    DNSFlags(cfbasics::Real nu = 0.0025, cfbasics::Real dPdx = 0.0, cfbasics::Real dPdz = 0.0, cfbasics::Real Ubulk = 0.0, cfbasics::Real Wbulk = 0.0, cfbasics::Real Uwall = 1.0,
             cfbasics::Real ulowerwall = 0.0, cfbasics::Real uupperwall = 0.0, cfbasics::Real wlowerwall = 0.0, cfbasics::Real wupperwall = 0.0,
             cfbasics::Real theta = 0.0, cfbasics::Real Vsuck = 0.0, cfbasics::Real rotation = 0.0, cfbasics::Real t0 = 0.0, cfbasics::Real T = 20.0, cfbasics::Real dT = 1.0,
             cfbasics::Real dt = 0.03125, bool variabledt = true, cfbasics::Real dtmin = 0.001, cfbasics::Real dtmax = 0.2, cfbasics::Real CFLmin = 0.4,
             cfbasics::Real CFLmax = 0.6, cfbasics::Real symmetryprojectioninterval = 100.0, BaseFlow baseflow = LaminarBase,
             MeanConstraint constraint = PressureGradient, TimeStepMethod timestepping = SBDF3,
             TimeStepMethod initstepping = SMRK2, NonlinearMethod nonlinearity = Rotational,
             Dealiasing dealiasing = DealiasXZ, BodyForce* bodyforce = 0, bool taucorrection = true,
             Verbosity verbosity = PrintTicks, std::ostream* logstream = &std::cout);

    DNSFlags(cfbasics::ArgList& args, const bool laurette = false);

    /** \brief The infamous virtual destructor */
    virtual ~DNSFlags() = default;

    bool dealias_xz() const;
    bool dealias_y() const;

    virtual void save(const std::string& outdir = "") const;  // save into file filebase.txt
    virtual void load(int taskid, const std::string indir);
    virtual const std::vector<std::string> getFlagList();

    BaseFlow baseflow;             // utot = u + Ubase(y) ex
    MeanConstraint constraint;     // Enforce const press grad or const bulk vel
    TimeStepMethod timestepping;   // Time-stepping algorithm
    TimeStepMethod initstepping;   // Algorithm for initializing multistep methods
    NonlinearMethod nonlinearity;  // Method of calculating nonlinearity of NS eqn
    Dealiasing dealiasing;         // Use 3/2 rule to eliminate aliasing
    BodyForce* bodyforce;          // Body force, zero if pointer set to 0
    bool taucorrection;            // Remove divergence caused by discretization

    cfbasics::Real nu;                            // Kinematic viscosity nu
    cfbasics::Real Vsuck;                         // suction velocity
    cfbasics::Real rotation;                      // dimensionless rotation around the z-axis
    cfbasics::Real theta;                         // tilt of the domain relative to downstream
    cfbasics::Real dPdx;                          // Constraint value for mean flow: pressure gradient in x
    cfbasics::Real dPdz;                          // Constraint value for mean flow: pressure gradient in z
    cfbasics::Real Ubulk;                         // Constraint value for mean flow: bulk velocity in x
    cfbasics::Real Wbulk;                         // Constraint value for mean flow: bulk velocity in z
    cfbasics::Real Uwall;                         // wall speed downstream
    cfbasics::Real ulowerwall;                    // lower wall speed along x, e.g. -1 for plane couette
    cfbasics::Real uupperwall;                    // upper wall speed along x, e.g. +1 for plane couette
    cfbasics::Real wlowerwall;                    // lower wall speed along z
    cfbasics::Real wupperwall;                    // upper wall speed along z
    cfbasics::Real t0;                            // start time
    cfbasics::Real T;                             // final time
    cfbasics::Real dT;                            // print interval
    cfbasics::Real dt;                            // time step
    bool variabledt;                    // use variable time step
    cfbasics::Real dtmin;                         // lower bound for time step
    cfbasics::Real dtmax;                         // upper bound for time step
    cfbasics::Real CFLmin;                        // lower bound for CFL number
    cfbasics::Real CFLmax;                        // upper bound for CFL number
    int symmetryprojectioninterval;     // Only project onto symmetries at this interval
    Verbosity verbosity;                // Print diagnostics, times, ticks, or nothing
    std::ostream* logstream;            // stream for output
    cfbasics::cfarray<FieldSymmetry> symmetries;  // restrict u(t) to these symmetries

   protected:
    void args2BC(cfbasics::ArgList& args);
    void args2numerics(cfbasics::ArgList& args, const bool laurette = false);
};

std::ostream& operator<<(std::ostream& os, const DNSFlags& flags);

// TimeStep keeps dt between dtmin and dtmax, and CFL between CFLminand CFLmax,
// in such a way that dt*n = dT for some integer n. That's useful if you
// want to plot/save data at regular dT intervals, but use a variable timestep
// dt for efficiency. You can mandate a fixed timestep by setting dtmin==dtmax.
// For example of use, see example codes.

class TimeStep {
   public:
    TimeStep();
    TimeStep(cfbasics::Real dt, cfbasics::Real dtmin, cfbasics::Real dtmax, cfbasics::Real dT, cfbasics::Real CFLmin, cfbasics::Real CFLmax, bool variable = true);
    TimeStep(DNSFlags& flags);

    // If variable, adjust dt to keep CFLmin<=CFL<=CFLmax (soft),
    // and dtmin<=dt<=dtmax (hard). Returns true if dt changes, false otherwise
    bool adjust(cfbasics::Real CFL, bool verbose = true, std::ostream& os = std::cout);
    bool adjustToMiddle(cfbasics::Real CFL, bool verbose = true, std::ostream& os = std::cout);

    // to bring any variable * dt below a maximum
    bool adjust(cfbasics::Real a, cfbasics::Real a_max, bool verbose = true, std::ostream& os = std::cout);
    bool adjustToDesired(cfbasics::Real a, cfbasics::Real a_des, bool verbose = true, std::ostream& os = std::cout);

    // tweak dT and dt to fit T exactly
    bool adjust_for_T(cfbasics::Real T, bool verbose = true, std::ostream& os = std::cout);

    int n() const;    // n*dt == dT
    int N() const;    // N*dT == T
    cfbasics::Real dt() const;  // integration timestep
    cfbasics::Real dtmin() const;
    cfbasics::Real dtmax() const;
    cfbasics::Real dT() const;  // plot/CFL-check interval
    cfbasics::Real T() const;   // total integration time
    cfbasics::Real CFL() const;
    cfbasics::Real CFLmin() const;
    cfbasics::Real CFLmax() const;
    bool variable() const;
    operator cfbasics::Real() const;  // same as dt()

   private:
    int n_;
    int N_;
    cfbasics::Real dt_;
    cfbasics::Real dtmin_;  //
    cfbasics::Real dtmax_;  //
    cfbasics::Real dT_;     // dT_ == n_*dt_, plot interval
    cfbasics::Real T_;      // T_  == N_*dt_, total integration time
    cfbasics::Real CFLmin_;
    cfbasics::Real CFL_;
    cfbasics::Real CFLmax_;
    bool variable_;
};

std::ostream& operator<<(std::ostream& os, const TimeStep& ts);

}  // namespace channelflow
#endif
