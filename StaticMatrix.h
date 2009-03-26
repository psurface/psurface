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
