/**
 * Statistics of the settling-suspension state (4-component field: velocity and the settling scalar,
 * stored as c + y, in component 3): turbulent kinetic energy, dissipation, Reynolds stresses,
 * sediment fluxes, higher moments and the concentration statistics written by simulateflow.
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#ifndef SETTLINGSTATS_H
#define SETTLINGSTATS_H

#include <string>
#include "channelflow/flowfield.h"

namespace chflow {

FlowField TKE(const FlowField& u);
FlowField dissipation_field(const FlowField& u);
FlowField sed_dissipation_field(const FlowField& u);
FlowField turbulent_convection(const FlowField& u);
FlowField sediment_fluctuation_convection(const FlowField& u);
FlowField pressure_transport(const FlowField& u, const FlowField& p);
FlowField Reynolds_stress_uu(const FlowField& u);
FlowField Reynolds_stress_uv(const FlowField& u);
FlowField Reynolds_stress_uw(const FlowField& u);
FlowField Reynolds_stress_vv(const FlowField& u);
FlowField Reynolds_stress_vw(const FlowField& u);
FlowField Reynolds_stress_ww(const FlowField& u);
FlowField u_skewness(const FlowField& u);
FlowField u_kurtosis(const FlowField& u);
FlowField turbulent_sediment_flux_uc(const FlowField& u);
FlowField turbulent_sediment_flux_vc(const FlowField& u);
FlowField turbulent_sediment_flux_cc(const FlowField& u);
FlowField c_skewness(const FlowField& u);
FlowField c_kurtosis(const FlowField& u);
FlowField intermittency_factor(const FlowField& u, double thresh);
void TKE(const FlowField& u, FlowField& e);
void dissipation_field(const FlowField& u, FlowField& e);
void turbulent_convection(const FlowField& u, FlowField& e);
void pressure_transport(const FlowField& u, const FlowField& p, FlowField& e);
void Reynolds_stress_uu(const FlowField& u, FlowField& e);
void Reynolds_stress_uv(const FlowField& u, FlowField& e);
void Reynolds_stress_uw(const FlowField& u, FlowField& e);
void Reynolds_stress_vv(const FlowField& u, FlowField& e);
void Reynolds_stress_vw(const FlowField& u, FlowField& e);
void Reynolds_stress_ww(const FlowField& u, FlowField& e);
void u_skewness(const FlowField& u, FlowField& e);
void u_kurtosis(const FlowField& u, FlowField& e);
void c_skewness(const FlowField& u, FlowField& e);
void c_kurtosis(const FlowField& u, FlowField& e);
void turbulent_sediment_flux_uc(const FlowField& u, FlowField& e);
void turbulent_sediment_flux_vc(const FlowField& u, FlowField& e);
void intermittency_factor(const FlowField& u, FlowField& e, double thresh);
Real getbulk(const FlowField& u, int i);
Real getCbulk(const FlowField& u);
// bulk mean of (u_i + U) (c - y), U = y the Couette base flow (i = 0 only)
Real getScalarCorrelation(const FlowField& u, FlowField& vel, FlowField& uc, int i);
std::string concstats_t(const FlowField& u, FlowField& vel, FlowField& uc, Real t);           // return some information about c
std::string concstatsheader_t(const std::string tname = "t");  // header for concstats
std::string concstats(const FlowField& u, FlowField& vel, FlowField& uc);                     // return some information about c
std::string concstatsheader();                                 // header for concstats

}  // namespace chflow
#endif  // SETTLINGSTATS_H
