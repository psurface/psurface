#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include "Vector.h"
#include "StaticVector.h"
#include <stdexcept>
#include <sstream>

namespace psurface {

/** A template class for sparse matrices. 
 *
 * \tparam T one of float, double, complex<float> or complex<double>.
*/
template<class T> class SparseMatrix
{
protected:
    struct MatrixEntry{
    MatrixEntry() : value(T(0)), col(0){};

    MatrixEntry(const T& newVal, int column){
        value = newVal;
        col = column;
    }

    T value;
    int col;
    };

    ///
    std::vector<std::vector<MatrixEntry> > data;

    ///
    int numCols;

public:
    /// Default Constructor
    SparseMatrix() : numCols(0) {
    data.clear();
    }

    ///
    SparseMatrix(int n) : numCols(n) {
    data.resize(n);
    for (int i=0; i<n; i++){
        data[i].resize(1);
        data[i][0] = MatrixEntry(T(0), i);
    }
    }

    /// Multiplication with a scalar
    void operator*=(const T& scalar) {
    for (size_t i=0; i<data.size(); i++)
        for (size_t j=0; j<data[i].size(); j++)
        data[i][j].value *= scalar;
    }

    /// The number of rows of the matrix
    size_t nRows() const {
    return data.size();
    }

    /// The number of columns of the matrix
    size_t nCols() const {
    return numCols;
    }

    ///
    void setEntry(int i, int j, const T& newValue) {
    for (size_t k=0; k<data[i].size(); k++)
        if (data[i][k].col==j){
        data[i][k].value = newValue;
        return;
        }

    data[i].push_back(MatrixEntry(newValue, j));
    }

    ///
    void addToEntry(int i, int j, const T& newValue) {
    for (size_t k=0; k<data[i].size(); k++)
        if (data[i][k].col==j){
        data[i][k].value += newValue;
        return;
        }

    data[i].push_back(MatrixEntry(newValue, j));
    }

    ///
    Vector<T> multVec(const Vector<T>& v) const {
        assert(v.size()==nCols());

        Vector<T> result(v.size(), StaticVector<T, 2>(0));

        for (size_t i=0; i<nRows(); i++)
            for (size_t j=0; j<data[i].size(); j++)
                result[i] += data[i][j].value * v[data[i][j].col];

      return result;
    }

    /// another iterative solver for nonsymmetric matrices: BI-CGSTAB
    size_t BiCGSTAB(const Vector<T>& b, Vector<T>& x, Vector<T>& r,
                  const size_t& maxIter, const T& tolerance) const {
      const float EPSILON=1e-20;
      const size_t n = numCols;

      T rho, rho_old, alpha, beta, omega, h, norm0;
      Vector<T> r0(n), p(n, StaticVector<T, 2>(0)), v(n, StaticVector<T, 2>(0)), s(n), t(n);

      r = r0 = b - multVec(x);
      norm0 = r0.length();

      rho_old = alpha = omega = 1;

      if (r.length() < (tolerance * norm0) || r.length() < 1e-15)
        return 0;

      for (size_t k = 0; k < maxIter; ++k) {
        rho = r0 * r;

        if (std::abs(rho) <= EPSILON) {
          std::ostringstream os;
          os << "Breakdown in BiCGSTAB - rho "
             << rho << " <= EPSILON " << EPSILON
             << " after " << k << " iterations";
          throw std::runtime_error(os.str());
        }

        if (std::abs(omega) <= EPSILON) {
          std::ostringstream os;
          os << "Breakdown in BiCGSTAB - omega "
             << omega << " <= EPSILON " << EPSILON
             << " after " << k << " iterations";
          throw std::runtime_error(os.str());
        }

        beta = (rho/rho_old)*(alpha/omega);
        p = r + beta*(p - omega * v);

        v = multVec(p);
        h = r0 * v;

        if (std::abs(h) < EPSILON) {
          std::ostringstream os;
          os << "Breakdown in BiCGSTAB - h "
             << h << " < EPSILON " << EPSILON
             << " after " << k << " iterations";
          throw std::runtime_error(os.str());
        }

        alpha = rho/h;
        x += alpha * p;

        s = r - alpha * v;

        if (s.length() < (tolerance * norm0))
          return k;

        t = multVec(s);
        omega = (t*s)/(t*t);

        x += omega * s;
        r = s - omega*t;

        rho_old = rho;

        if (r.length() < (tolerance * norm0)  || r.length() < 1e-15)
          return k;
      }

      throw std::runtime_error("BiCGSTAB did not converge.");
    }
};

} // namespace psurface

#endif
