#ifndef STATIC_VECTOR_H
#define STATIC_VECTOR_H

#ifdef _MSC_VER
    // Required to make cmath define M_PI etc.
    #define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <cstring>  // for size_t

#include <array>
#include <iostream>
#include <assert.h>

namespace psurface {

template <class T, int N>
class StaticVector
    : public std::array<T,N>
{
public:

    /** \brief Default constructor.  Leaves the vector uninitialized. */
    StaticVector()
        : std::array<T,N>()
    {}

    /** \brief Construction from a single scalar */
    explicit StaticVector(const T& s) {
        this->fill(s);
    }

    /** \brief Construction from a two scalars */
    StaticVector(const T& a, const T& b) {
        assert(N==2);
        (*this)[0] = a;
        (*this)[1] = b;
    }

    /** \brief Construction from a three scalars */
    StaticVector(const T& a, const T& b, const T& c) {
        assert(N==3);
        (*this)[0] = a;
        (*this)[1] = b;
        (*this)[2] = c;
    }

    /** \brief Copy constructor */
    StaticVector(const StaticVector<T,N>& other) {
        for (int i=0; i<N; i++)
            (*this)[i] = other[i];
    }


    /** \brief Vector product */
    T dot(const StaticVector<T,N>& a) const {
        T result = 0;
        for (size_t i=0; i<N; i++)
            result += (*this)[i] * a[i];
        return result;
    }

    /// Cross Product.
    StaticVector<T,N> cross(const StaticVector<T,N>& o) const {
        assert(N==3);
        StaticVector<T,N> result;
        
        result[0] = (*this)[1]*o[2]-(*this)[2]*o[1];
        result[1] = (*this)[2]*o[0]-(*this)[0]*o[2];
        result[2] = (*this)[0]*o[1]-(*this)[1]*o[0];
        return result;
    }

    /** \brief Euclidean length */
    T length() const {
        return std::sqrt(dot(*this));
    }

    /** \brief Euclidean length squared */
    T length2() const {
        return dot(*this);
    }

    /** \brief Make unit vector */
    void normalize() {
        T l = length();
        for (size_t i=0; i<N; i++)
            (*this)[i] /= l;
    }

    /// returns the angle between two vectors
    T angle(const StaticVector<T,N>& other) const {
        assert(length()!=0 && other.length()!=0);
        const T ratio = dot(other) / (length()*other.length());
        if (ratio<-1) return M_PI;
        return (ratio>1) ? 0 : std::acos(ratio);
    }

    /** \brief Assignment */
    StaticVector<T,N>& operator=(const StaticVector<T,N>& other) {
        for (size_t i=0; i<N; i++)
            (*this)[i] = other[i];
        return *this;
    }

    /** \brief Addition */
    StaticVector<T,N>& operator+=(const StaticVector<T,N>& other) {
        for (size_t i=0; i<N; i++)
            (*this)[i] += other[i];
        return *this;
    }

    /** \brief Division */
    StaticVector<T,N>& operator/=(const T& divisor) {
        for (size_t i=0; i<N; i++)
            (*this)[i] /= divisor;
        return *this;
    }

    /** \brief Unary minus */
    friend StaticVector<T,N> operator-(const StaticVector<T,N>& a) {
        StaticVector<T,N> result;
        for (size_t i=0; i<N; i++)
            result[i] = -a[i];
        return result;
    }

    /** \brief Addition */
    friend StaticVector<T,N> operator+(const StaticVector<T,N>& a, const StaticVector<T,N>& b) {
        StaticVector<T,N> result;
        for (size_t i=0; i<N; i++)
            result[i] = a[i] + b[i];
        return result;
    }

    /** \brief Subtraction */
    friend StaticVector<T,N> operator-(const StaticVector<T,N>& a, const StaticVector<T,N>& b) {
        StaticVector<T,N> result;
        for (size_t i=0; i<N; i++)
            result[i] = a[i] - b[i];
        return result;
    }

    /** \brief Scalar multiplication from the left */
    friend StaticVector<T,N> operator*(const T& s, const StaticVector<T,N>& a) {
        StaticVector<T,N> result;
        for (size_t i=0; i<N; i++)
            result[i] = s * a[i];
        return result;
    }

    /** \brief Scalar multiplication from the right */
    friend StaticVector<T,N> operator*(const StaticVector<T,N>& a, const T& s) {
        StaticVector<T,N> result;
        for (size_t i=0; i<N; i++)
            result[i] = s * a[i];
        return result;
    }

    /** \brief Scalar division */
    friend StaticVector<T,N> operator/(const StaticVector<T,N>& a, const T& s) {
        StaticVector<T,N> result;
        for (size_t i=0; i<N; i++)
            result[i] = a[i] / s;
        return result;
    }
};
    template<class P, int K>
    std::ostream& operator<< (std::ostream& s, const StaticVector<P,K>& v)
    {
        for (int i=0; i<K; i++)
           s << ((i>0) ? " " : "") << v[i];
         return s;
    }


} // namespace psurface

#endif
