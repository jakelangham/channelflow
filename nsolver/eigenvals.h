/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 *
 * Original authors: Ayse Yesil, Mirko Farano
 */
#ifndef NSOLVER_EIGENVALS_H
#define NSOLVER_EIGENVALS_H

#include "cfbasics/arglist.h"
#include "cfbasics/cfbasics.h"
#include "nsolver/arnoldi.h"
#include "nsolver/dsi.h"
#include "nsolver/lanczos.h"

namespace nsolver {

class EigenvalsFlags {
   public:
    EigenvalsFlags();
    EigenvalsFlags(cfbasics::ArgList& args);

    std::ostream* logstream = &std::cout;

    bool isnormal = false;

    int Narnoldi = 100;
    int Nstable = 5;
    bool fixedNs = false;

    cfbasics::Real EPS_kry = 1e-10;
    bool centdiff = false;
    bool orthochk = false;

    std::string duname = "";
    std::string outdir = "./";

    // new flags:
    cfbasics::Real EPS_stab = 1e-06;

    void save(const std::string& outdir = "") const;  // save into file filebase.txt
    void load(int taskid, const std::string indir);
};

class Eigenvals {
   public:
    Eigenvals(cfbasics::ArgList& args);
    Eigenvals(EigenvalsFlags eigenflags);

    void solve(nsolver::DSI& dsi, const Eigen::VectorXd& x, Eigen::VectorXd& dx, cfbasics::Real T, cfbasics::Real eps);
    void checkConjugacy(const Eigen::VectorXcd& u, const Eigen::VectorXcd& v);
    std::ostream* getLogstream() { return eigenflags.logstream; }

   private:
    EigenvalsFlags eigenflags;
};

std::ostream& operator<<(std::ostream& os, const EigenvalsFlags& flags);

}  // namespace nsolver

#endif
