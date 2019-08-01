/* JL complex dense matrix class with UL (not LU) decomposition solver.
 * We could use LU but keeping UL to match bandedtridiag and hopefully
 * forthcoming septdiagonal matrix class. */
// TODO a complex septdiagonal matrix class would be more efficient
 
#include "channelflow/realdensematrix.h"

namespace chflow {

RealDenseMatrix::~RealDenseMatrix() {
    delete[] data_;
}

RealDenseMatrix::RealDenseMatrix() : M_(0), UL_(false), data_(0) {}

RealDenseMatrix::RealDenseMatrix(int M)
    : M_(M), UL_(false), data_(new Real[2 * M * M]) {
    assert(M >= 0);
    for (int i = 0; i < 2 * M * M; i++) {
        data_[i] = 0.0;
    }
}

RealDenseMatrix::RealDenseMatrix(const RealDenseMatrix& A)
    : M_(A.M_), UL_(A.UL_), data_(new Real[2 * M_ * M_]) {
    for (int i = 0; i < 2 * M_ * M_; i++) {
        data_[i] = A.data_[i];
    }
}

RealDenseMatrix& RealDenseMatrix::operator=(const RealDenseMatrix& A) {
    if (this != &A) {
        if (M_ != A.M_) {
            delete[] data_;
            M_ = A.M_;
            UL_ = A.UL_;
            data_ = new Real[2 * M_ * M_];
        }
        for (int i = 0; i < 2 * M_ * M_; i++) {
            data_[i] = A.data_[i];
        }
    }

    return *this;
}

const Real& RealDenseMatrix::elem(int i, int j) const {
    assert(i >= 0 && j >= 0 && i <= M_ && j <= M_);
    assert(UL_ == false);

    return data_[i + j * M_];
}

Real& RealDenseMatrix::elem(int i, int j) {
    assert(i >= 0 && j >= 0 && i <= M_ && j <= M_);
    assert(UL_ == false);

    return data_[i + j * M_];
}

// Multiply matrix by x and store in b. How to do the multiplication depends on
// whether the UL decomposition has been done yet.
void RealDenseMatrix::multiply(const ChebyCoeff& x, ChebyCoeff&b) const {
    Real sum;
    assert(M_ == x.length());

    if (UL_ == false) {
        for (int i = 0; i < M_; i++) {
            sum = 0.0;
            for (int k = 0; k < M_; k++) {
                sum += elem(i, k) * x[k];
            }
            b[i] = sum;
        }
    } else {
        // Here we have to multiply by UL instead
        for (int i = 0; i < M_; i++) {
            sum = 0.0;
            for (int k = 0; k < M_; k++) {
                for (int l = 0; l < M_; l++) {
                    sum += U(i, l) * L(l, k) * x[k];
                }
            }
            b[i] = sum;
        }
    }
}

// UL decomposition without pivoting.
// We are free to choose L_ii = 1
void RealDenseMatrix::ULdecomp() {
    // Start from the bottom row and work up
    for (int i = M_ - 1; i >= 0; i--) {
        L(i, i) = 1.0;

        // solve for U entries
        for (int j = M_ - 1; j >= i; j--) {
            U(i, j) = elem(i, j);
            for (int k = M_ - 1; k > j; k--) {
                U(i, j) -= U(i, k) * L(k, j);
            }
        }

        // solve for L entries
        for (int j = i - 1; j >= 0; j--) {
            L(i, j) = elem(i, j);
            for (int k = M_ - 1; k > i; k--) {
                L(i, j) -= U(i, k) * L(k, j);
            }
            L(i, j) /= U(i, i);

            // This zeros U's lower triangle
            U(i, j) = 0.0;
        }
    }

    UL_ = true;
}

// Solve Ax = b given UL = A, via Uy = b, then Lx = y.
// N.B. In this case, b is a set of complex Chebyshev coefficients.
void RealDenseMatrix::ULsolve(ChebyCoeff& u) const {
    assert(UL_ == true);
    assert(M_ == u.N());

    Real b_i, y_i; // element of b, y

    // Solve Uy = b by back substitution, iterating last row to zeroth.
    for (int i = M_ - 1; i >= 0; i--) {
        b_i = u[i];
        for (int j = M_ - 1; j > i; j--) {
            b_i -= U(i, j) * u[j];
        }
        b_i /= U(i, i);
        u[i] = b_i;
    }

    // Solve Lx = y by forward substitution. First row needs no calculation
    // due to sparsity structure. N.B. different convention to BandedTridiag
    for (int i = 1; i < M_; i++) {
        y_i = u[i];
        for (int j = 0; j < i; j++) {
            y_i -= L(i, j) * u[j];
        }
        u[i] = y_i;
    }
}

void RealDenseMatrix::print() const {
    assert(UL_ == false);
    cout << "[\n";
    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j < M_; ++j) {
            // old format
            //cout << elem(i, j) << ' ';
            // julia format
            cout << elem(i, j) << " ";
        }
        cout << ";\n";
    }
    cout << "]\n";
}

void RealDenseMatrix::ULprint() const {
    assert(UL_ == true);

    cout << "U = [\n";
    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j < M_; ++j) {
            // old format
            //cout << U(i, j) << ' ';
            // julia format
            cout << U(i, j) << " ";
        }
        cout << ";\n";
    }
    cout << "]\n";

    cout << "L = [\n";
    for (int i = 0; i < M_; ++i) {
        for (int j = 0; j < M_; ++j) {
            // old format
            //cout << L(i, j) << ' ';
            // julia format
            cout << L(i, j) << " ";
        }
        cout << ";\n";
    }
    cout << "]\n";
}

} // namespace chflow
