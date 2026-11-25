/**
 * This file is a part of channelflow version 2.0.
 * License is GNU GPL version 2 or later: https://channelflow.org/license
 */
#include "modules/settling/settlingstats.h"

#include <iomanip>
#include <sstream>
#include "channelflow/diffops.h"

using namespace std;

namespace chflow {

Real getScalarCorrelation(const FlowField& u, FlowField& vel, FlowField& uc, int i) {
    Real ucbulk;

    assert(i >= 0 && i <= 2);
    //FlowField vel(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    //FlowField uc(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    //FlowField vel, uc;
    ChebyCoeff Ubase = ChebyCoeff(u.My(), u.a(), u.b(), Spectral);
    ChebyCoeff Cbase = ChebyCoeff(u.My(), u.a(), u.b(), Spectral);
    Ubase[1] = 1;
    Cbase[1] = -1;

    vel = u[i];
    uc = u[3];

    if (i == 0)
        vel += Ubase;
    uc += Cbase;

    vel.makePhysical();
    uc.makePhysical();
    uc *= vel;
    uc.makeSpectral();
#ifdef HAVE_MPI
    if (u.taskid() == u.task_coeff(0, 0))
        ucbulk = Re(uc.profile(0, 0, 0)).mean();
    MPI_Bcast(&ucbulk, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
#else
    ucbulk = Re(uc.profile(0, 0, 0)).mean();
#endif

    if (abs(ucbulk) < 1e-15)
        ucbulk = 0.0;
    return ucbulk;
}

void TKE(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    FlowField up(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    FlowField wp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    up = u[0];
    wp = u[2];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            up.cmplx(0, my, 0, 0) = 0.0;
            wp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    //Real ubar = getUbulk(u), wbar = getWbulk(u);
    //ChebyCoeff ubar = Re(u.profile(0, 0, 0));
    //ChebyCoeff wbar = Re(u.profile(0, 0, 2));
    u.makePhysical();
    up.makePhysical();
    wp.makePhysical();
    //ubar.makePhysical(); wbar.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    //     int Nx = u.Nx();
    //     int Ny = u.Ny();
    lint Nz = u.Nz();
    //lint Nd = u.Nd();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

//    for (lint i = 0; i < 0; ++i)
    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += 0.5 * square(up(nx, ny, nz, 0));
    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += 0.5 * square(u(nx, ny, nz, 1));
    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += 0.5 * square(wp(nx, ny, nz, 0));

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;
}

void dissipation_field(const FlowField& u_, FlowField& e) {
    //FlowField& u = (FlowField&)u_;
    FlowField up(u_);
    FlowField ui, grad_ui;
    //fieldstate sxz = u.xzstate();
    //fieldstate sy = u.ystate();

    //ChebyCoeff uprof = Re(u.profile(0, 0, 0));
    //ChebyCoeff wprof = Re(u.profile(0, 0, 2));
    //ChebyCoeff dubardy = diff(uprof);
    //ChebyCoeff dwbardy = diff(wprof);
    for (int my = 0; my < up.My(); ++my) {
        if (0 >= up.mxlocmin() && 0 < up.Mxloc() + up.mxlocmin() && 
            0 >= up.mzlocmin() && 0 < up.Mzloc() + up.mzlocmin()) {
            up.cmplx(0, my, 0, 0) = 0.0;
            up.cmplx(0, my, 0, 2) = 0.0;
        }
    }

    //dubardy.makePhysical();
    //dwbardy.makePhysical();

    up.makePhysical();

    if (!up.geomCongruent(e) || e.Nd() != 1)
        e.resize(up.Nx(), up.Ny(), up.Nz(), 1, up.Lx(), up.Lz(), up.a(), up.b(), up.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    lint Nz = up.Nz();
    lint nxlocmin = up.nxlocmin();
    lint nxlocmax = up.nxlocmin() + up.Nxloc();
    lint nylocmin = up.nylocmin();
    lint nylocmax = up.nylocmax();

    // u 
    ui = up[0];
    grad(ui, grad_ui);
    grad_ui.makePhysical();
    for (lint ny = nylocmin; ny < nylocmax; ++ny) {
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx) {
            for (int nz = 0; nz < Nz; ++nz) {
                e(nx, ny, nz, 0) += grad_ui(nx, ny, nz, 0) * grad_ui(nx, ny, nz, 0) +
                                    grad_ui(nx, ny, nz, 1) * grad_ui(nx, ny, nz, 1) + 
                                    grad_ui(nx, ny, nz, 2) * grad_ui(nx, ny, nz, 2);
            }
        }
    }

    // v 
    ui = up[1];
    grad(ui, grad_ui);
    grad_ui.makePhysical();
    for (lint ny = nylocmin; ny < nylocmax; ++ny) {
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx) {
            for (int nz = 0; nz < Nz; ++nz) {
                e(nx, ny, nz, 0) += grad_ui(nx, ny, nz, 0) * grad_ui(nx, ny, nz, 0) +
                                    grad_ui(nx, ny, nz, 1) * grad_ui(nx, ny, nz, 1) + 
                                    grad_ui(nx, ny, nz, 2) * grad_ui(nx, ny, nz, 2);
            }
        }
    }

    // w 
    ui = up[2];
    grad(ui, grad_ui);
    grad_ui.makePhysical();
    for (lint ny = nylocmin; ny < nylocmax; ++ny) {
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx) {
            for (int nz = 0; nz < Nz; ++nz) {
                e(nx, ny, nz, 0) += grad_ui(nx, ny, nz, 0) * grad_ui(nx, ny, nz, 0) +
                                    grad_ui(nx, ny, nz, 1) * grad_ui(nx, ny, nz, 1) + 
                                    grad_ui(nx, ny, nz, 2) * grad_ui(nx, ny, nz, 2);
            }
        }
    }

    e.makeSpectral();
    //u.makeState(sxz, sy);
    return;
}

void sed_dissipation_field(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    FlowField cp = u[3];
    FlowField grad_cp;
    //fieldstate sxz = u.xzstate();
    //fieldstate sy = u.ystate();

    for (int my = 0; my < cp.My(); ++my) {
        if (0 >= cp.mxlocmin() && 0 < cp.Mxloc() + cp.mxlocmin() && 
            0 >= cp.mzlocmin() && 0 < cp.Mzloc() + cp.mzlocmin()) {
            cp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    grad(cp, grad_cp);
    grad_cp.makePhysical();

    if (!cp.geomCongruent(e) || e.Nd() != 1)
        e.resize(cp.Nx(), cp.Ny(), cp.Nz(), 1, cp.Lx(), cp.Lz(), cp.a(), cp.b(), cp.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    lint Nz = cp.Nz();
    lint nxlocmin = cp.nxlocmin();
    lint nxlocmax = cp.nxlocmin() + cp.Nxloc();
    lint nylocmin = cp.nylocmin();
    lint nylocmax = cp.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny) {
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx) {
            for (int nz = 0; nz < Nz; ++nz) {
                e(nx, ny, nz, 0) += grad_cp(nx, ny, nz, 0) * grad_cp(nx, ny, nz, 0) +
                                    grad_cp(nx, ny, nz, 1) * grad_cp(nx, ny, nz, 1) + 
                                    grad_cp(nx, ny, nz, 2) * grad_cp(nx, ny, nz, 2);
            }
        }
    }

    e.makeSpectral();
    //u.makeState(sxz, sy);
    return;
}

void turbulent_convection(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    FlowField up(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    FlowField wp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    up = u[0];
    wp = u[2];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            up.cmplx(0, my, 0, 0) = 0.0;
            wp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    //Real ubar = getUbulk(u), wbar = getWbulk(u);
    //ChebyCoeff ubar = Re(u.profile(0, 0, 0));
    //ChebyCoeff wbar = Re(u.profile(0, 0, 2));
    u.makePhysical();
    up.makePhysical();
    wp.makePhysical();
    //ubar.makePhysical(); wbar.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    //     int Nx = u.Nx();
    //     int Ny = u.Ny();
    lint Nz = u.Nz();
    //lint Nd = u.Nd();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny) {
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx) {
            for (int nz = 0; nz < Nz; ++nz) {
                e(nx, ny, nz, 0) += 0.5 * u(nx, ny, nz, 1) * 
                    (up(nx, ny, nz, 0) * up(nx, ny, nz, 0) + 
                     u(nx, ny, nz, 1) * u(nx, ny, nz, 1) +
                     wp(nx, ny, nz, 0) * wp(nx, ny, nz, 0));
            }
        }
    }

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;

}

void sediment_fluctuation_convection(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    FlowField cp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    cp = u[3];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            cp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    u.makePhysical();
    cp.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    //     int Nx = u.Nx();
    //     int Ny = u.Ny();
    lint Nz = u.Nz();
    //lint Nd = u.Nd();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny) {
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx) {
            for (int nz = 0; nz < Nz; ++nz) {
                e(nx, ny, nz, 0) += 0.5 * u(nx, ny, nz, 1) * 
                    (cp(nx, ny, nz, 0) * cp(nx, ny, nz, 0));
            }
        }
    }

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;
}

void pressure_transport(const FlowField& u_, const FlowField& p_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    FlowField& p = (FlowField&)p_;
    fieldstate sxz = u.xzstate(), sxz_p = p.xzstate();
    fieldstate sy = u.ystate(), sy_p = p.ystate();

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    u.makePhysical();
    p.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    // p is assumed to be true pressure (not q = p + 0.5 u dot u)
    for (lint ny = nylocmin; ny < nylocmax; ++ny) {
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx) {
            for (int nz = 0; nz < Nz; ++nz) {
                e(nx, ny, nz, 0) += u(nx, ny, nz, 1) * p(nx, ny, nz, 0);
            }
        }
    }

    e.makeSpectral();

    u.makeState(sxz, sy);
    p.makeState(sxz_p, sy_p);
    return;
}

void Reynolds_stress_uu(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    FlowField up(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    up = u[0];
    for (int my = 0; my < u.My(); ++my)
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin())
            up.cmplx(0, my, 0, 0) = 0.0;

    u.makePhysical();
    up.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += up(nx, ny, nz, 0) * up(nx, ny, nz, 0);

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;

}

void u_skewness(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    FlowField up(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    up = u[0];
    for (int my = 0; my < u.My(); ++my)
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin())
            up.cmplx(0, my, 0, 0) = 0.0;

    u.makePhysical();
    up.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += up(nx, ny, nz, 0) * up(nx, ny, nz, 0) * up(nx, ny, nz, 0);

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;
}

void u_kurtosis(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    FlowField up(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    up = u[0];
    for (int my = 0; my < u.My(); ++my)
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin())
            up.cmplx(0, my, 0, 0) = 0.0;

    u.makePhysical();
    up.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += up(nx, ny, nz, 0) * up(nx, ny, nz, 0) * up(nx, ny, nz, 0) * up(nx, ny, nz, 0);

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;
}

void Reynolds_stress_uv(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    FlowField up(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    up = u[0];
    for (int my = 0; my < u.My(); ++my)
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin())
            up.cmplx(0, my, 0, 0) = 0.0;

    u.makePhysical();
    up.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += up(nx, ny, nz, 0) * u(nx, ny, nz, 1);

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;

}

void Reynolds_stress_uw(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;

    FlowField up(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    FlowField wp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    up = u[0]; wp = u[2];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            up.cmplx(0, my, 0, 0) = 0.0;
            wp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    up.makePhysical();
    wp.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += up(nx, ny, nz, 0) * wp(nx, ny, nz, 0);

    e.makeSpectral();
    return;

}

void Reynolds_stress_vv(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    u.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += u(nx, ny, nz, 1) * u(nx, ny, nz, 1);

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;
}

void Reynolds_stress_vw(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    FlowField wp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    wp = u[2];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            wp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    u.makePhysical();
    wp.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += u(nx, ny, nz, 1) * wp(nx, ny, nz, 0);

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;

}

void Reynolds_stress_ww(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    FlowField wp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    wp = u[2];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            wp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    u.makePhysical();
    wp.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += wp(nx, ny, nz, 0) * wp(nx, ny, nz, 0);

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;

}

void turbulent_sediment_flux_uc(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    //fieldstate sxz = u.xzstate();
    //fieldstate sy = u.ystate();

    FlowField up(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    FlowField cp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    up = u[0];
    cp = u[3];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            up.cmplx(0, my, 0, 0) = 0.0;
            cp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    //Real ubar = getUbulk(u), cbar = getCbulk(u);
    //ChebyCoeff ubar = Re(u.profile(0, 0, 0));
    //ChebyCoeff cbar = Re(u.profile(0, 0, 3));
    //u.makePhysical();
    up.makePhysical();
    cp.makePhysical();
    //ubar.makePhysical(); cbar.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    //     int Nx = u.Nx();
    //     int Ny = u.Ny();
    lint Nz = u.Nz();
    //lint Nd = u.Nd();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += up(nx, ny, nz, 0) * cp(nx, ny, nz, 0);

    e.makeSpectral();
    //u.makeState(sxz, sy);
    return;

}

void turbulent_sediment_flux_vc(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    fieldstate sxz = u.xzstate();
    fieldstate sy = u.ystate();

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    FlowField cp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    cp = u[3];
    for (int my = 0; my < u.My(); ++my)
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin())
            cp.cmplx(0, my, 0, 0) = 0.0;

    cp.makePhysical();
    //ChebyCoeff cbar = Re(u.profile(0, 0, 3));
    //cbar.makePhysical();
    u.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += u(nx, ny, nz, 1) * cp(nx, ny, nz, 0);

    e.makeSpectral();
    u.makeState(sxz, sy);
    return;
}

void turbulent_sediment_flux_cc(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;
    //fieldstate sxz = u.xzstate();
    //fieldstate sy = u.ystate();

    FlowField cp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    cp = u[3];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            cp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    cp.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    //     int Nx = u.Nx();
    //     int Ny = u.Ny();
    lint Nz = u.Nz();
    //lint Nd = u.Nd();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += cp(nx, ny, nz, 0) * cp(nx, ny, nz, 0);

    e.makeSpectral();
    //u.makeState(sxz, sy);
    return;
}

void c_skewness(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;

    FlowField cp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    cp = u[3];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            cp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    cp.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += cp(nx, ny, nz, 0) * cp(nx, ny, nz, 0) * cp(nx, ny, nz, 0);

    e.makeSpectral();
    return;
}

void c_kurtosis(const FlowField& u_, FlowField& e) {
    FlowField& u = (FlowField&)u_;

    FlowField cp(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi(), Spectral, Spectral);
    cp = u[3];
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            cp.cmplx(0, my, 0, 0) = 0.0;
        }
    }

    cp.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    for (lint ny = nylocmin; ny < nylocmax; ++ny)
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx)
            for (int nz = 0; nz < Nz; ++nz)
                e(nx, ny, nz, 0) += cp(nx, ny, nz, 0) * cp(nx, ny, nz, 0) * cp(nx, ny, nz, 0) * cp(nx, ny, nz, 0);

    e.makeSpectral();
    return;
}

void intermittency_factor(const FlowField& u_, FlowField& e, double thresh) {
    FlowField& u = (FlowField&)u_;

    lint Nz = u.Nz();
    lint nxlocmin = u.nxlocmin();
    lint nxlocmax = u.nxlocmin() + u.Nxloc();
    lint nylocmin = u.nylocmin();
    lint nylocmax = u.nylocmax();

    FlowField up = u;
    for (int my = 0; my < u.My(); ++my) {
        if (0 >= u.mxlocmin() && 0 < u.Mxloc() + u.mxlocmin() && 
            0 >= u.mzlocmin() && 0 < u.Mzloc() + u.mzlocmin()) {
            up.cmplx(0, my, 0, 0) = 0.0;
            up.cmplx(0, my, 0, 2) = 0.0;
        }
    }

    //FlowField vort = curl(up);
    //vort.makePhysical();
    up.makePhysical();

    if (!u.geomCongruent(e) || e.Nd() != 1)
        e.resize(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi());
    else
        e.setToZero();
    e.setState(Physical, Physical);

    for (lint ny = nylocmin; ny < nylocmax; ++ny) {
        for (lint nx = nxlocmin; nx < nxlocmax; ++nx) {
            for (int nz = 0; nz < Nz; ++nz) {
                //double enstrophy = vort(nx, ny, nz, 0) * vort(nx, ny, nz, 0) +
                //                   vort(nx, ny, nz, 1) * vort(nx, ny, nz, 1) +
                //                   vort(nx, ny, nz, 2) * vort(nx, ny, nz, 2);
                //e(nx, ny, nz, 0) = (enstrophy > thresh) ? 1 : 0;
                double tke = 0.5 * (up(nx, ny, nz, 0) * up(nx, ny, nz, 0) +
                                    up(nx, ny, nz, 1) * up(nx, ny, nz, 1) + 
                                    up(nx, ny, nz, 2) * up(nx, ny, nz, 2));
                e(nx, ny, nz, 0) = (tke > thresh) ? 1 : 0;
            }
        }
    }

    e.makeSpectral();
}

FlowField TKE(const FlowField& u) {
    FlowField e;
    TKE(u, e);
    return e;
}

FlowField dissipation_field(const FlowField& u) {
    FlowField e;
    dissipation_field(u, e);
    return e;
}

FlowField sed_dissipation_field(const FlowField& u) {
    FlowField e;
    sed_dissipation_field(u, e);
    return e;
}

FlowField turbulent_convection(const FlowField& u) {
    FlowField e;
    turbulent_convection(u, e);
    return e;
}

FlowField sediment_fluctuation_convection(const FlowField& u) {
    FlowField e;
    sediment_fluctuation_convection(u, e);
    return e;
}

FlowField pressure_transport(const FlowField& u, const FlowField& p) {
    FlowField e;
    pressure_transport(u, p, e);
    return e;
}

FlowField Reynolds_stress_uu(const FlowField& u) {
    FlowField e;
    Reynolds_stress_uu(u, e);
    return e;
}

FlowField Reynolds_stress_uv(const FlowField& u) {
    FlowField e;
    Reynolds_stress_uv(u, e);
    return e;
}

FlowField Reynolds_stress_uw(const FlowField& u) {
    FlowField e;
    Reynolds_stress_uw(u, e);
    return e;
}

FlowField Reynolds_stress_vv(const FlowField& u) {
    FlowField e;
    Reynolds_stress_vv(u, e);
    return e;
}

FlowField Reynolds_stress_vw(const FlowField& u) {
    FlowField e;
    Reynolds_stress_vw(u, e);
    return e;
}

FlowField Reynolds_stress_ww(const FlowField& u) {
    FlowField e;
    Reynolds_stress_ww(u, e);
    return e;
}

FlowField u_skewness(const FlowField& u) {
    FlowField e;
    u_skewness(u, e);
    return e;
}

FlowField u_kurtosis(const FlowField& u) {
    FlowField e;
    u_kurtosis(u, e);
    return e;
}

FlowField turbulent_sediment_flux_uc(const FlowField& u) {
    FlowField e;
    turbulent_sediment_flux_uc(u, e);
    return e;
}

FlowField turbulent_sediment_flux_vc(const FlowField& u) {
    FlowField e;
    turbulent_sediment_flux_vc(u, e);
    return e;
}

FlowField turbulent_sediment_flux_cc(const FlowField& u) {
    FlowField e;
    turbulent_sediment_flux_cc(u, e);
    return e;
}

FlowField c_skewness(const FlowField& u) {
    FlowField e;
    c_skewness(u, e);
    return e;
}

FlowField c_kurtosis(const FlowField& u) {
    FlowField e;
    c_kurtosis(u, e);
    return e;
}

FlowField intermittency_factor(const FlowField& u, double thresh) {
    FlowField e;
    intermittency_factor(u, e, thresh);
    return e;
}

Real getbulk(const FlowField& u, int i) {
    Real ubulk;

#ifdef HAVE_MPI
    if (u.taskid() == u.task_coeff(0, 0))
        ubulk = Re(u.profile(0, 0, i)).mean();
    MPI_Bcast(&ubulk, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
#else
    ubulk = Re(u.profile(0, 0, i)).mean();
#endif

    if (abs(ubulk) < 1e-15)
        ubulk = 0.0;
    return ubulk;
}

Real getCbulk(const FlowField& u) {
    Real cbulk;

#ifdef HAVE_MPI
    if (u.taskid() == u.task_coeff(0, 0))
        cbulk = Re(u.profile(0, 0, 3)).mean();
    MPI_Bcast(&cbulk, 1, MPI_DOUBLE, u.task_coeff(0, 0), *u.comm_world());
#else
    cbulk = Re(u.profile(0, 0, 3)).mean();
#endif

    if (abs(cbulk) < 1e-15)
        cbulk = 0.0;
    return cbulk;
}

string concstatsheader() {
    stringstream header;
    header << setw(14) << "mass_loading" << setw(14) << "ucbulk" 
           << setw(14) << "vcbulk" << setw(14) << "wcbulk";
    return header.str();
}

string concstatsheader_t(const string tname) {
    stringstream header;
    header << setw(8) << "#(" << tname << ")" << concstatsheader();
    return header.str();
}

string concstats(const FlowField& u, FlowField& vel, FlowField& uc) {
    stringstream s;

    s << setw(14) << getCbulk(u) << setw(14) << getScalarCorrelation(u, vel, uc, 0)
      << setw(14) << getScalarCorrelation(u, vel, uc, 1) << setw(14) << getScalarCorrelation(u, vel, uc,  2);
    return s.str();
}

string stratstatsheader_t(const string tname) {
    stringstream header;
    header << setw(8) << "#(" << tname << ")" << setw(14) << "cbulk" << setw(14) << "c2/2" << setw(14) << "uc"
           << setw(14) << "vc" << setw(14) << "wc" << setw(14) << "chi" << setw(14) << "buoyflux";
    return header.str();
}

string stratstats_t(const FlowField& u4, const SettlingScalar& sc, Real t) {
    FlowField u, c;
    splitVelocityScalar(u4, u, c);
    u.makeSpectral();
    c.makeSpectral();
    Real flux[3];
    for (int i = 0; i < 3; ++i) {
        FlowField ui = u[i];
        flux[i] = L2InnerProduct(ui, c);
    }
    const Real chi = sc.kappa * L2Norm2(grad(c));
    stringstream s;
    s << setw(8) << t << setw(14) << getCbulk(u4) << setw(14) << 0.5 * L2Norm2(c) << setw(14) << flux[0]
      << setw(14) << flux[1] << setw(14) << flux[2] << setw(14) << chi << setw(14) << -sc.Ri * flux[sc.gravity()];
    return s.str();
}

// Return some statistics about the concentration field
string concstats_t(const FlowField& u, FlowField& vel, FlowField& uc, Real t) {
    stringstream s;
    s << setw(8);
    s << t;
    s << concstats(u, vel, uc);
    return s.str();
}

}  // namespace chflow
