/**
 * Test DenseHelmholtzSolver (the settling module's scalar solver) for every boundary-condition type:
 * solve nu u'' - lambda u = f for a smooth exact u, with the BC values (and, for NoFlux, the mean)
 * taken from u, and compare with u. Also check fixBC / fixNoFluxNormaliseMass, which must leave an
 * expansion already satisfying its BC unchanged and impose it on one that does not.
 *
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#include <cmath>
#include <iostream>
#include "channelflow/chebyshev.h"
#include "modules/settling/boundarycond.h"
#include "modules/settling/densehelmholtz.h"

using namespace std;
using namespace chflow;

// u = exp(p y) cos(q y) + s y on [-1, 1] and its first two derivatives
const Real p = 0.3, q = 1.7, s = 0.2;
Real uex(Real y) { return exp(p * y) * cos(q * y) + s * y; }
Real duex(Real y) { return exp(p * y) * (p * cos(q * y) - q * sin(q * y)) + s; }
Real d2uex(Real y) { return exp(p * y) * ((p * p - q * q) * cos(q * y) - 2 * p * q * sin(q * y)); }

ChebyCoeff spectral(int N, Real (*func)(Real)) {
    ChebyCoeff u(N, -1, 1, Physical);
    Vector y = chebypoints(N, -1, 1);
    for (int n = 0; n < N; ++n)
        u[n] = func(y[n]);
    u.makeSpectral();
    return u;
}

int main() {
    const int N = 33;
    const Real nu = 0.7, lambda = 3.1, alpha = 1.3;
    const Real tol = 1e-11;
    int failure = 0;

    ChebyCoeff u = spectral(N, uex);
    ChebyCoeff f(N, -1, 1, Physical);
    Vector y = chebypoints(N, -1, 1);
    for (int n = 0; n < N; ++n)
        f[n] = nu * d2uex(y[n]) - lambda * uex(y[n]);
    f.makeSpectral();

    const Real ua = uex(-1), ub = uex(1), dua = duex(-1), dub = duex(1);
    struct Case {
        string name;
        ScalarBCType type;
        Real g0, g1;
        bool mass;
    };
    const Case cases[] = {
        {"Diri", ScalarBCType::Diri, ua, ub, false},
        {"DiriRobin", ScalarBCType::DiriRobin, ua, dub + alpha * ub, false},
        {"NeumRobin", ScalarBCType::NeumRobin, dua, dub + alpha * ub, false},
        {"NoFlux", ScalarBCType::NoFlux, dua + alpha * ua, dub + alpha * ub, false},
        {"NoFlux+mass", ScalarBCType::NoFlux, dua + alpha * ua, dub + alpha * ub, true},
    };

    cout << "solver: max |u - u_exact| (coefficients)" << endl;
    for (const Case& c : cases) {
        BoundaryCond bc(c.type, 1, 1, c.g0, c.g1, alpha);
        DenseHelmholtzSolver solver(N, bc, -1, 1, lambda, nu, c.mass);
        ChebyCoeff v(N, -1, 1, Spectral);
        solver.solve(v, f, c.g0, c.g1, u.mean());
        Real err = 0;
        for (int n = 0; n < N; ++n)
            err = max(err, abs(v[n] - u[n]));
        cout << "  " << c.name << " " << err << endl;
        if (!(err < tol))
            failure = 1;
    }

    // fixBC imposes the inhomogeneous BC through coefficients 0 and 1 (Diri, DiriRobin, NoFlux)
    cout << "fixBC: BC residual after perturbing coefficients 0, 1" << endl;
    for (const Case& c : cases) {
        if (c.type == ScalarBCType::NeumRobin || c.mass)
            continue;  // fixBC(f, ga, gb, bc) has no NeumRobin branch; the mass is checked below
        BoundaryCond bc(c.type, 1, 1, c.g0, c.g1, alpha);
        ChebyCoeff v(u);
        v[0] += 0.1;
        v[1] -= 0.2;
        fixBC(v, c.g0, c.g1, bc);
        const Real r0 = (c.type == ScalarBCType::NoFlux) ? v.slope_a() + alpha * v.eval_a() - c.g0 : v.eval_a() - c.g0;
        const Real r1 = (c.type == ScalarBCType::Diri) ? v.eval_b() - c.g1 : v.slope_b() + alpha * v.eval_b() - c.g1;
        const Real err = max(abs(r0), abs(r1));
        cout << "  " << c.name << " " << err << endl;
        if (!(err < tol))
            failure = 1;
    }

    // fixNoFluxNormaliseMass imposes the no-flux BC and mean(u) = 1 through coefficients 0, 1, 2.
    // The stored settling scalar has g = 1 -/+ alpha for the mean mode.
    {
        ChebyCoeff v(u);
        v[0] += 0.1;
        v[1] -= 0.2;
        v[2] += 0.05;
        const Real ga = 1 - alpha, gb = 1 + alpha;
        fixNoFluxNormaliseMass(v, ga, gb, alpha);
        const Real err = max(max(abs(v.slope_a() + alpha * v.eval_a() - ga), abs(v.slope_b() + alpha * v.eval_b() - gb)),
                             abs(v.mean() - 1.0));
        cout << "fixNoFluxNormaliseMass: BC and mass residual " << err << endl;
        if (!(err < tol))
            failure = 1;
    }

    cout << (failure ? "\t** FAIL **" : "\t   pass   ") << endl;
    return failure;
}
