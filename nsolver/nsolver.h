/**
 * Includes all classed needed to interface NSolver
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#ifndef CHANNELFLOW_NSOLVER_H
#define CHANNELFLOW_NSOLVER_H

#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "Eigen/Dense"

#include "cfbasics/arglist.h"
#include "cfbasics/cfbasics.h"

#include "nsolver/arnoldi.h"
#include "nsolver/bicgstab.h"
#include "nsolver/bicgstabl.h"
#include "nsolver/continuation.h"
#include "nsolver/eigenvals.h"
#include "nsolver/fgmres.h"
#include "nsolver/gmres.h"
#include "nsolver/lanczos.h"
#include "nsolver/multiShootingDSI.h"
#include "nsolver/newton.h"
#include "nsolver/newtonalgorithm.h"

using namespace Eigen;

namespace nsolver {}  // namespace nsolver
#endif                // NSOLVER_H
