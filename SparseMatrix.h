#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include <complex>

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
    int nRows() const {
	return data.size();
    }

    /// The number of columns of the matrix
    int nCols() const {
	return numCols;
    }

    ///
    void setEntry(int i, int j, const T& newValue) {
	for (int k=0; k<data[i].size(); k++)
	    if (data[i][k].col==j){
		data[i][k].value = newValue;
		return;
	    }

	data[i].push_back(MatrixEntry(newValue, j));
    }

    ///
    void addToEntry(int i, int j, const T& newValue) {
	for (int k=0; k<data[i].size(); k++)
	    if (data[i][k].col==j){
		data[i][k].value += newValue;
		return;
	    }

	data[i].push_back(MatrixEntry(newValue, j));
    }

    /// Multiply real sparse matrix with complex vector.
    void multVecC(const std::vector< std::complex<T> >& v,
		  std::vector< std::complex<T> >& result) const {

	assert(v.size()==nCols());
	result.resize(nRows());
	result.assign(result.size(), std::complex<T>(0.));

        for (int i=0; i<nRows(); i++)
            for (size_t j=0; j<data[i].size(); j++)
                result[i] += data[i][j].value * v[data[i][j].col];

    }

    /// another iterative solver for nonsymmetric matrices: BI-CGSTAB
    void BiCGSTABC(const std::vector<std::complex<T> >& rhs, std::vector<std::complex<T> >& result,
		  std::vector<std::complex<T> >& residuum,
		  int* maxIter, const double tolerance) const {

	int i, j;

	int countIter=0;

	const int N = rhs.size();

	const std::vector<std::complex<T> >& b = rhs;
	std::vector<std::complex<T> > xi = result;

	std::vector<std::complex<T> > tmp;
	multVecC(xi, tmp);

	std::vector<std::complex<T> > ri(N);
	for (j=0; j<N; j++)
	    ri[j] = b[j] - tmp[j];

	double normRes = sqrt(std::real(dotProductC(ri, ri)));
	double normRhs = sqrt(std::real(dotProductC(rhs, rhs)));

	if (normRes <= tolerance * normRhs) {
	    residuum = ri;
	    *maxIter = 0;
	    return;
	}

	std::vector<std::complex<T> > r0hat = ri;

startAgain:

	std::complex<double> tmp2 = dotProductCRI(r0hat, r0hat);

	double normR0hatr = sqrt(std::real(tmp2));
	double normR0hati = sqrt(std::imag(tmp2));
	double normRi_r   = normR0hatr;
	double normRi_i   = normR0hati;

	double alphar = 1;
	double alphai = 1;
	double omegai_r = 1;
	double omegai_i = 1;

	std::vector<std::complex<T> > pi(N);
	std::vector<std::complex<T> > vi(N);
	pi.assign(pi.size(), std::complex<T>(0));
	vi.assign(pi.size(), std::complex<T>(0));

	double rhoiMin1r    = 1;
	double rhoiMin1i    = 1;

	tmp2 = dotProductCRI(r0hat, ri);
	double rhoi_r = std::real(tmp2);
	double rhoi_i = std::imag(tmp2);


	i=0;
	while (i++<*maxIter) {

	    if (fabs(rhoi_r) < 1.e-10 * normR0hatr * normRi_r ||
		fabs(rhoi_i) < 1.e-10 * normR0hati * normRi_i) {
		//printf("iter %d  rhoi: %g %g \n",i,rhoi_r,rhoi_i);
		r0hat = ri;
		goto startAgain;
	    }

	    countIter++;

	    double betar = (rhoi_r / rhoiMin1r)*(alphar / omegai_r);
	    double betai = (rhoi_i / rhoiMin1i)*(alphai / omegai_i);

	    // compute omegaiMin1*viMin1

	    for (j=0; j<N; j++)
		pi[j] = std::complex<T>( std::real(ri[j]) + betar*(std::real(pi[j]) - omegai_r*std::real(vi[j])) ,
				       std::imag(ri[j]) + betai*(std::imag(pi[j]) - omegai_i*std::imag(vi[j])) );
	    multVecC(pi, vi);

	    tmp2 = dotProductCRI(r0hat, vi);
	    alphar = rhoi_r / std::real(tmp2);
	    alphai = rhoi_i / std::imag(tmp2);

	    std::vector<std::complex<T> > s(N);

	    for (j=0; j<N; j++)
		s[j] = std::complex<T>( std::real(ri[j]) - alphar*std::real(vi[j]),
				      std::imag(ri[j]) - alphai*std::imag(vi[j]) );

	    if (sqrt(std::real(dotProductC(s,s))) <= tolerance * normRhs) {
		for (j=0; j<N; j++)
		    xi[j] = std::complex<T>( std::real(xi[j]) + alphar*std::real(pi[j]),
			    		   std::imag(xi[j]) + alphai*std::imag(pi[j]) );
		ri = s;
		break;
	    }

	    std::vector<std::complex<T> > t;
	    multVecC(s, t);

	    tmp2 = dotProductCRI(t,t);
	    double denomr = std::real(tmp2);
            double denomi = std::imag(tmp2);

	    tmp2 = dotProductCRI(t,s);
	    omegai_r = std::real(tmp2) / denomr;
	    omegai_i = std::imag(tmp2) / denomi;

	    rhoiMin1r   = rhoi_r;
	    rhoiMin1i   = rhoi_i;
	    tmp2 = dotProductCRI(r0hat, t);
	    rhoi_r = -omegai_r * std::real(tmp2);
	    rhoi_i = -omegai_i * std::imag(tmp2);

	    for (j=0; j<N; j++)
		xi[j] = std::complex<T>( std::real(xi[j]) + alphar*std::real(pi[j]) + omegai_r*std::real(s[j]),
				       std::imag(xi[j]) + alphai*std::imag(pi[j]) + omegai_i*std::imag(s[j]) );

	    for (j=0; j<N; j++)
		ri[j] = std::complex<T>( std::real(s[j]) - omegai_r*std::real(t[j]),
				       std::imag(s[j]) - omegai_i*std::imag(t[j]) );

	    tmp2 = dotProductCRI(ri,ri);
	    normRi_r = sqrt(std::real(tmp2));
	    normRi_i = sqrt(std::imag(tmp2));
	    normRes = sqrt(std::real(tmp2) + std::imag(tmp2));

	    if (normRes <= tolerance * normRhs)
		break;

	    ////////////////////////
	}

	*maxIter = countIter;

	residuum = ri;

	result = xi;

    }

    static std::complex<T> dotProductC(const std::vector< std::complex<T> >& a,
			 const std::vector< std::complex<T> >& b){
	assert(a.size()==b.size());

	std::complex<T> sum = std::complex<T>(0);

	for (size_t i=0; i<a.size(); i++)
	    sum += a[i]*std::conj(b[i]);

	return sum;
    }

    static std::complex<double> dotProductCRI(const std::vector< std::complex<T> >& a,
			 const std::vector< std::complex<T> >& b){
	assert(a.size()==b.size());

	std::complex<double> sum = std::complex<double>(0);

	for (size_t i=0; i<a.size(); i++)
	    sum += std::complex<double>( (double)std::real(a[i])*std::real(b[i]) , (double)std::imag(a[i])*std::imag(b[i]) );

	return sum;
    }

};

/// @endif

} // namespace psurface

#endif
