#ifndef STATIC_MATRIX_H
#define STATIC_MATRIX_H



template <class T, int N>
class StaticMatrix
    : public std::tr1::array<StaticVector<T,N>, N>
{
public:

    /** \brief Default constructor, doesn't initialize anything */
    StaticMatrix() {}

    /** \brief Construct matrix from column vectors */
    StaticMatrix(const StaticVector<T,N>& a, const StaticVector<T,N>& b, const StaticVector<T,N>& c)
    {
        (*this)[0][0]=a[0]; (*this)[0][1]=b[0]; (*this)[0][2]=c[0];
        (*this)[1][0]=a[1]; (*this)[1][1]=b[1]; (*this)[1][2]=c[1];
        (*this)[2][0]=a[2]; (*this)[2][1]=b[2]; (*this)[2][2]=c[2];

    }

    /** \brief Constructor for 3x3 matrices from separate scalars */
    StaticMatrix( const T& a00, const T& a01, const T& a02,
                  const T& a10, const T& a11, const T& a12,
                  const T& a20, const T& a21, const T& a22 )
    {
        assert(N==3);
        (*this)[0][0]=a00; (*this)[0][1]=a01; (*this)[0][2]=a02;
        (*this)[1][0]=a10; (*this)[1][1]=a11; (*this)[1][2]=a12;
        (*this)[2][0]=a20; (*this)[2][1]=a21; (*this)[2][2]=a22;
    }


    /// Computes the determinant of the matrix.
    T det() const {
        T ad1 = (*this)[0][0] * ((*this)[1][1]*(*this)[2][2] - (*this)[1][2]*(*this)[2][1]);
        T ad2 = (*this)[0][1] * ((*this)[1][0]*(*this)[2][2] - (*this)[1][2]*(*this)[2][0]);
        T ad3 = (*this)[0][2] * ((*this)[1][0]*(*this)[2][1] - (*this)[1][1]*(*this)[2][0]);

        return ad1 - ad2 + ad3;
    }

    StaticMatrix<T,N> inverse() const {

        StaticMatrix<T,N> result;
        assert(N==3);
        T d = det();

        result[0][0] =  ((*this)[1][1]*(*this)[2][2] - (*this)[1][2]*(*this)[2][1]) / d;
        result[0][1] = -((*this)[0][1]*(*this)[2][2] - (*this)[0][2]*(*this)[2][1]) / d;
        result[0][2] =  ((*this)[0][1]*(*this)[1][2] - (*this)[0][2]*(*this)[1][1]) / d;
        result[1][0] = -((*this)[1][0]*(*this)[2][2] - (*this)[1][2]*(*this)[2][0]) / d;
        result[1][1] =  ((*this)[0][0]*(*this)[2][2] - (*this)[0][2]*(*this)[2][0]) / d;
        result[1][2] = -((*this)[0][0]*(*this)[1][2] - (*this)[0][2]*(*this)[1][0]) / d;
        result[2][0] =  ((*this)[1][0]*(*this)[2][1] - (*this)[1][1]*(*this)[2][0]) / d;
        result[2][1] = -((*this)[0][0]*(*this)[2][1] - (*this)[0][1]*(*this)[2][0]) / d;
        result[2][2] =  ((*this)[0][0]*(*this)[1][1] - (*this)[0][1]*(*this)[1][0]) / d;

        return result;
    }

    void multMatrixVec(const StaticVector<T,N>& src, StaticVector<T,N>& dst) const {
        assert(N==3);
        dst[0] = src[0]*(*this)[0][0]+src[1]*(*this)[0][1]+src[2]*(*this)[0][2];
        dst[1] = src[0]*(*this)[1][0]+src[1]*(*this)[1][1]+src[2]*(*this)[1][2];
        dst[2] = src[0]*(*this)[2][0]+src[1]*(*this)[2][1]+src[2]*(*this)[2][2];
    }
};

#endif
