/**
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */
#include "modules/settling/settlingdsi.h"
#include "channelflow/utilfuncs.h"

using namespace std;
using namespace Eigen;

namespace chflow {

// The modes follow the velocity's packing (Gibson, Halcrow & Cvitanovic 2009, table 1): (0,0), then
// (kx,0) for 0 < kx <= Kx, (0,kz) for 0 < kz <= Kz, (kx,kz) for kx != 0, 0 < kz <= Kz. The (-kx,0)
// modes are conjugates of (kx,0) and are not stored. Each process packs the modes it owns.

int scalarVectorSize(const FlowField& c) {
    const int Kx = c.kxmaxDealiased(), Kz = c.kzmaxDealiased(), Ny = c.Ny();
    int N = 0;
    if (c.taskid() == c.task_coeff(0, 0))
        N += Ny - 3;  // c'' + ... : 0, 1 from the no-flux BC, 2 from the mass constraint
    for (int kx = 1; kx <= Kx; ++kx)
        if (c.taskid() == c.task_coeff(c.mx(kx), 0))
            N += 2 * (Ny - 2);  // 0, 1 from the no-flux BC
    for (int kz = 1; kz <= Kz; ++kz)
        if (c.taskid() == c.task_coeff(0, c.mz(kz)))
            N += 2 * (Ny - 2);
    for (int kx = -Kx; kx <= Kx; ++kx) {
        if (kx == 0)
            continue;
        for (int kz = 1; kz <= Kz; ++kz)
            if (c.taskid() == c.task_coeff(c.mx(kx), c.mz(kz)))
                N += 2 * (Ny - 2);
    }
    return N;
}

void scalarToVector(const FlowField& c, VectorXd& x, int& n) {
    assert(c.xzstate() == Spectral && c.ystate() == Spectral && c.Nd() == 1);
    const int Kx = c.kxmaxDealiased(), Kz = c.kzmaxDealiased(), Ny = c.Ny();
    auto packmode = [&](int mx, int mz) {
        for (int ny = 2; ny < Ny; ++ny) {
            x(n++) = Re(c.cmplx(mx, ny, mz, 0));
            x(n++) = Im(c.cmplx(mx, ny, mz, 0));
        }
    };
    if (c.taskid() == c.task_coeff(0, 0))
        for (int ny = 3; ny < Ny; ++ny)
            x(n++) = Re(c.cmplx(0, ny, 0, 0));
    for (int kx = 1; kx <= Kx; ++kx)
        if (c.taskid() == c.task_coeff(c.mx(kx), 0))
            packmode(c.mx(kx), 0);
    for (int kz = 1; kz <= Kz; ++kz)
        if (c.taskid() == c.task_coeff(0, c.mz(kz)))
            packmode(0, c.mz(kz));
    for (int kx = -Kx; kx <= Kx; ++kx) {
        if (kx == 0)
            continue;
        for (int kz = 1; kz <= Kz; ++kz)
            if (c.taskid() == c.task_coeff(c.mx(kx), c.mz(kz)))
                packmode(c.mx(kx), c.mz(kz));
    }
}

void vectorToScalar(const VectorXd& x, int& n, FlowField& c, const SettlingScalar& sc) {
    assert(c.Nd() == 1);
    c.setToZero();
    c.setState(Spectral, Spectral);
    const int Kx = c.kxmaxDealiased(), Kz = c.kzmaxDealiased(), Ny = c.Ny();
    const BoundaryCond bc = sc.bc(c.Mx(), c.Mz());
    ComplexChebyCoeff f(Ny, c.a(), c.b(), Spectral);

    // unpack one non-mean mode and restore coefficients 0, 1 from the homogeneous no-flux BC
    auto unpackmode = [&](int mx, int mz) {
        f.setToZero();
        for (int ny = 2; ny < Ny; ++ny) {
            f.re[ny] = x(n++);
            f.im[ny] = x(n++);
        }
        fixBC(f, bc);
        for (int ny = 0; ny < Ny; ++ny)
            c.cmplx(mx, ny, mz, 0) = f[ny];
    };

    if (c.taskid() == c.task_coeff(0, 0)) {
        f.setToZero();
        for (int ny = 3; ny < Ny; ++ny)
            f.re[ny] = x(n++);
        fixNoFluxNormaliseMass(f.re, bc.ga_, bc.gb_, bc.alpha_);  // 0, 1 from the BC, 2 from mean = 1
        for (int ny = 0; ny < Ny; ++ny)
            c.cmplx(0, ny, 0, 0) = f[ny];
    }

    for (int kx = 1; kx <= Kx; ++kx) {
        const int mx = c.mx(kx), mxm = c.mx(-kx);
        const int send_id = c.task_coeff(mx, 0), rec_id = c.task_coeff(mxm, 0);
        if (c.taskid() == send_id)
            unpackmode(mx, 0);
        // the (-kx,0) mode is the conjugate of (kx,0), stored for FFTW
        for (int ny = 0; ny < Ny; ++ny) {
            if (c.taskid() == send_id && send_id == rec_id)
                c.cmplx(mxm, ny, 0, 0) = conj(c.cmplx(mx, ny, 0, 0));
#ifdef HAVE_MPI
            else {
                if (c.taskid() == send_id) {
                    Complex tmp = conj(c.cmplx(mx, ny, 0, 0));
                    MPI_Send(&tmp, 1, MPI_DOUBLE_COMPLEX, rec_id, 3, MPI_COMM_WORLD);
                }
                if (c.taskid() == rec_id) {
                    Complex tmp;
                    MPI_Status status;
                    MPI_Recv(&tmp, 1, MPI_DOUBLE_COMPLEX, send_id, 3, MPI_COMM_WORLD, &status);
                    c.cmplx(mxm, ny, 0, 0) = tmp;
                }
            }
#endif
        }
    }

    for (int kz = 1; kz <= Kz; ++kz)
        if (c.taskid() == c.task_coeff(0, c.mz(kz)))
            unpackmode(0, c.mz(kz));

    for (int kx = -Kx; kx <= Kx; ++kx) {
        if (kx == 0)
            continue;
        for (int kz = 1; kz <= Kz; ++kz)
            if (c.taskid() == c.task_coeff(c.mx(kx), c.mz(kz)))
                unpackmode(c.mx(kx), c.mz(kz));
    }
    c.setPadded(true);
}

SettlingDSI::SettlingDSI(DNSFlags& dnsflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch,
                         bool xrelative, bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u,
                         ostream* os)
    : cfDSI(dnsflags, sigma, h, dt, Tsearch, xrelative, zrelative, Tnormalize, Unormalize, u, os) {
    if (u.Nd() != 4)
        cferror("SettlingDSI : the state must be a 4-component field (velocity and settling scalar)");
    if (h != 0)
        cferror("SettlingDSI : Poincare sections are not supported by the settling module");
}

int SettlingDSI::fieldVectorSize(const FlowField& u4) {
    FlowField u, c;
    splitVelocityScalar(u4, u, c);
    return field2vector_size(u) + scalarVectorSize(c);
}

void SettlingDSI::fieldToVector(const FlowField& u4, VectorXd& x) {
    // As field2vector: grow x if needed, zero it (extra unknowns such as T come after the field part)
    FlowField u, c;
    splitVelocityScalar(u4, u, c);
    VectorXd xu;
    field2vector(u, xu);
    int n = field2vector_size(u);
    const int N = n + scalarVectorSize(c);
    if (x.size() < N)
        x.resize(N);
    x.setZero();
    x.head(n) = xu.head(n);
    scalarToVector(c, x, n);
}

void SettlingDSI::vectorToField(const VectorXd& x, FlowField& u4) {
    FlowField u, c;
    splitVelocityScalar(u4, u, c);
    vector2field(x, u);
    int n = field2vector_size(u);
    vectorToScalar(x, n, c, settlingScalars(dnsflags_)[0]);
    combineVelocityScalar(u, c, u4);
}

void SettlingDSI::integrate(const FlowField& u4, Real& T, PoincareCondition* h, FlowField& fu, const DNSFlags& flags_,
                            const TimeStep& dt_, int& fcount, Real& CFL, ostream& os) {
    // f^T(u) for the settling equations; follows the fixed-T branch of the free function f
    if (h != 0)
        cferror("SettlingDSI::integrate : Poincare sections are not supported by the settling module");
    if (!isfinite(L2Norm(u4))) {
        os << "error in f: u is not finite. exiting." << endl;
        exit(1);
    }
    DNSFlags flags(flags_);
    flags.logstream = &os;
    TimeStep dt(dt_);
    fu = u4;
    if (T < 0) {
        os << "f: negative integration time T == " << T << endl << "returning f(u,T) == (1+abs(T))*u" << endl;
        fu *= 1 + abs(T);
        return;
    }
    if (T == 0) {
        os << "f: T==0, no integration, returning u " << endl;
        return;
    }

    FlowField u, c;
    splitVelocityScalar(u4, u, c);
    vector<FlowField> fields = {u, c, FlowField(u.Nx(), u.Ny(), u.Nz(), 1, u.Lx(), u.Lz(), u.a(), u.b(), u.cfmpi())};

    dt.adjust_for_T(T, false);
    flags.dt = dt;
    SettlingDNS dns(fields, flags);
    if (dt.variable()) {
        dns.advance(fields, 1);
        dt.adjust(dns.CFL(fields[0]), false);
        dns.reset_dt(dt);
        fields[0] = u;
        fields[1] = c;
    }

    os << "f^T: " << flush;
    for (int s = 1; s <= dt.N(); ++s) {
        Real t = s * dt.dT();
        CFL = dns.CFL(fields[0]);
        if (s % 10 == 0)
            os << iround(t) << flush;
        else if (s % 2 == 0) {
            if (CFL < dt.CFLmin())
                os << '<' << flush;
            else if (CFL > dt.CFLmax())
                os << '>' << flush;
            else
                os << '.' << flush;
        }
        dns.advance(fields, dt.n());
        if (dt.variable() && dt.adjust(CFL, false))
            dns.reset_dt(dt);
    }

    combineVelocityScalar(fields[0], fields[1], fu);
    if (!isfinite(L2Norm(fu))) {
        os << "error in f: f(u,t) is not finite. exiting." << endl;
        exit(1);
    }
    ++fcount;
}

}  // namespace chflow
