#ifndef CHANNELFLOW_MATRIX_H
#define CHANNELFLOW_MATRIX_H

#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"

using namespace std;

namespace chflow {

class RealDenseMatrix {
public:
    RealDenseMatrix();
    RealDenseMatrix(int M);
    RealDenseMatrix(const RealDenseMatrix &A);
    RealDenseMatrix& operator=(const RealDenseMatrix& A);
    ~RealDenseMatrix();
    const Real& elem(int i, int j) const;
    Real& elem(int i, int j);
    void multiply(const ChebyCoeff&, ChebyCoeff&) const;
    inline Real& L(int i, int j);
    inline Real& U(int i, int j);
    inline const Real& L(int i, int j) const;
    inline const Real& U(int i, int j) const;
    void ULdecomp();
    void ULsolve(ChebyCoeff& u) const;
    void print() const;
    void ULprint() const;

private:
    int M_; // # rows/cols
    bool UL_; // has UL decomp been done?
    Real* data_; 
};

inline Real& RealDenseMatrix::L(int i, int j) {
    assert(i >= 0 && j >= 0 && i <= M_ && j <= M_);

    return data_[i + j * M_ + M_ * M_];
}
inline Real& RealDenseMatrix::U(int i, int j) {
    assert(i >= 0 && j >= 0 && i <= M_ && j <= M_);

    return data_[i + j * M_];
}
inline const Real& RealDenseMatrix::L(int i, int j) const {
    assert(i >= 0 && j >= 0 && i <= M_ && j <= M_);

    return data_[i + j * M_ + M_ * M_];
}
inline const Real& RealDenseMatrix::U(int i, int j) const {
    assert(i >= 0 && j >= 0 && i <= M_ && j <= M_);

    return data_[i + j * M_];
}

} // namespace chflow

#endif // CHANNELFLOW_MATRIX_H
