#ifndef VECTOR_H
#define VECTOR_H

#include "StaticVector.h"
#include <vector>
#include <cmath>
#include <assert.h>

namespace psurface {
  template <typename T>
  class Vector
    : public std::vector<StaticVector<T, 2> >
  {
  public:
    typedef StaticVector<T, 2> VType;

    /** \brief Default constructor.  Leaves the vector uninitialized. */
    Vector(const int& n)
      : std::vector<VType>(n)
    {}

    /** \brief Construction from a single scalar */
    explicit Vector(const int& n, const VType& s)
      : std::vector<VType>(n)
    {
      this->assign(n, s);
    }

    /** \brief Copy constructor */
    Vector(const Vector<T>& other)
      : std::vector<VType>(other.size())
    {
      for (int i=0; i<this->size(); i++)
        (*this)[i] = other[i];
    }

    /** \brief Addition */
    Vector<T>& operator+=(const Vector<T>& other) {
      assert(other.size() == this->size());

      for (size_t i=0; i<this->size(); i++)
        (*this)[i] += other[i];

      return *this;
    }

    /** \brief Division */
    Vector<T>& operator/=(const T& divisor) {
      for (size_t i=0; i<this->size(); i++)
        (*this)[i] /= divisor;

      return *this;
    }

    /** \brief Unary minus */
    friend Vector<T> operator-(const Vector<T>& a) {
      Vector<T> result(a.size());

      for (size_t i=0; i<a.size(); i++)
        result[i] = -a[i];

      return result;
    }

    /** \brief Addition */
    friend Vector<T> operator+(const Vector<T>& a, const Vector<T>& b) {
      assert(a.size() == b.size());

      Vector<T> result(a.size());

      for (size_t i=0; i<a.size(); i++)
        result[i] = a[i] + b[i];

      return result;
    }

    /** \brief Subtraction */
    friend Vector<T> operator-(const Vector<T>& a, const Vector<T>& b) {
      assert(a.size() == b.size());

      Vector<T> result(a.size());

      for (size_t i=0; i<a.size(); i++)
        result[i] = a[i] - b[i];

      return result;
    }

    /** \brief Vector Product */
    T operator*(const Vector<T>& other) const {
      assert(this->size() == other.size());

      T result = 0;

      for (size_t i=0; i<this->size(); i++)
        result += (*this)[i].dot(other[i]);

      return result;
    }

    /** \brief Scalar multiplication from the left */
    friend Vector<T> operator*(const T& s, const Vector<T>& a) {
      Vector<T> result(a.size());

      for (size_t i=0; i<a.size(); i++)
        result[i] = s * a[i];

      return result;
    }

    /** \brief Scalar multiplication from the right */
    friend Vector<T> operator*(const Vector<T>& a, const T& s) {
      return s * a;
    }

    /** \brief Scalar division */
    friend Vector<T> operator/(const Vector<T>& a, const T& s) {
      Vector<T> result(a.size());

      for (size_t i=0; i<a.size(); i++)
        result[i] = a[i] / s;

      return result;
    }

    /** \brief Vector product */
    T dot(const Vector<T>& a) const {
      return (*this) * a;
    }

    /** \brief Euclidean length */
    T length() const {
      return sqrt(dot(*this));
    }

    /** \brief Euclidean length squared */
    T length2() const {
      return dot(*this);
    }
  };
} // namespace psurface

#endif
