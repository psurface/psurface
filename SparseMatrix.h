#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <vector>
#include <complex>

/** A template class for sparse matrices.  The first template parameter
    should be one of float, double, complex<float> or complex<double>.
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

    ///
    SparseMatrix(int rows, int columns) : numCols(columns) {
	data.resize(rows);
    }

    /// Index operator.  Only for element reading.
    T operator()(int i, int j) const {
	assert(i>=0 && i<data.size() && j>=0 && j<numCols);
	for (int k=0; k<data[i].size(); k++)
	    if (data[i][k].col==j)
		return data[i][k].value;

	return T(0);
    }

    /// Multiplication with a scalar
    void operator*=(const T& scalar) {
	int i, j;
	for (i=0; i<data.size(); i++)
	    for (j=0; j<data[i].size(); j++)
		data[i][j].value *= scalar;
    }

    ///
    void init(const int n) {
	numCols = n;
	data.resize(n);
	for (int i=0; i<n; i++){
	    data[i].resize(1);
	    data[i][0] = MatrixEntry(T(0), i);
	}
    }

    ///
    void resize(int m, int n) {
	assert(m>=0 && n>=0);

	data.resize(m);
	int i, j;
	if (numCols>n)
	    for (i=0; i<m; i++)
		for (j=data[i].size()-1; j>=0; j--)
		    if (data[i][j].col>=n)
			data[i].remove(j);

	numCols = n;
    }

    ///
    void resizeRows(int n){
	assert(n>=0);
	data.resize(n);
    }

    /// The number of rows of the matrix
    int nRows() const {
	return data.size();
    }

    /// The number of columns of the matrix
    int nCols() const {
	return numCols;
    }

    /// The number of elements
    int nElements() const {
	int nElem = 0;
	int i;
	for (i=0; i<nRows(); i++)
	    nElem += data[i].size();
	return nElem;
    }

    /// Quadratic matrix or not?
    bool isQuadratic() const {
	return nRows()==nCols();
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

    /// Multiply sparse matrix with vector.
    void multVec(const std::vector<T>& v, std::vector<T>& result) const {
	int i,j,k;

	assert(v.size()==nCols());
	result.resize(nRows());
	result.assign(result.size(), T(0));

        for (i=0; i<nRows(); i++)
            for (j=0; j<data[i].size(); j++)
                result[i] += data[i][j].value * v[data[i][j].col];

    }

    /// Multiply real sparse matrix with complex vector.
    void multVecC(const std::vector< std::complex<T> >& v,
		  std::vector< std::complex<T> >& result) const {
	int i,j,k;

	assert(v.size()==nCols());
	result.resize(nRows());
	result.assign(result.size(), std::complex<T>(0.));

        for (i=0; i<nRows(); i++)
            for (j=0; j<data[i].size(); j++)
                result[i] += data[i][j].value * v[data[i][j].col];

    }

    /// Multiply sparse matrix with multiple vectors.
    void multVecM(const std::vector< std::vector<T> >& v,
		  std::vector< std::vector<T> >& result) const {
	int i,j,k,l;

	int nRhs = v.size();
	if (nRhs == 0) return;

	result.resize(nRhs);
	for (i=0; i<nRhs; i++) {
	    assert(v[i].size()==nCols());
	    result[i].resize(nRows());
	    result[i].fill(T(0));
	}

        for (i=0; i<nRows(); i++)
            for (j=0; j<data[i].size(); j++)
                for (l=0; l<nRhs; l++)
                    result[l][i] += data[i][j].value * v[l][data[i][j].col];

    }

    /// another iterative solver for nonsymmetric matrices: BI-CGSTAB
    void BiCGSTAB(const std::vector<T>& rhs, std::vector<T>& result,
		  std::vector<T>& residuum,
		  int* maxIter, const double tolerance) const {

	int j;
	int countIter=0;
	const int N = rhs.size();

	const std::vector<T>& b = rhs;
	std::vector<T> xi = result;

	std::vector<T> tmp;
	multVec(xi, tmp);

	std::vector<T> ri(N);
	for (j=0; j<N; j++)
	    ri[j] = b[j] - tmp[j];

	double normRes = sqrt(realPart(dotProduct(ri, ri)));
	double normRhs = sqrt(realPart(dotProduct(rhs, rhs)));

	if (normRes <= tolerance * normRhs) {
	    residuum = ri;
	    *maxIter = 0;
	    return;
	}

	std::vector<T> r0hat = ri;

startAgain:

	double normR0hat = sqrt(realPart(dotProduct(r0hat,r0hat)));
	double normRi = normR0hat;

	double alpha = 1;
	double omegai = 1;

	std::vector<T> pi(N);
	std::vector<T> vi(N);
        std::fill(pi.begin(), pi.end(), T(0));
        std::fill(vi.begin(), vi.end(), T(0));

	double rhoiMin1    = 1;

	double rhoi = realPart(dotProduct(ri, ri));

	while (countIter++<*maxIter) {

    	    if (fabs(rhoi) < 1.e-10 * normR0hat * normRi) {
		//printf("iter %d  rhoi: %g \n",i,rhoi);
		r0hat = ri;
		goto startAgain;
	    }

	    double beta = (rhoi / rhoiMin1)*(alpha / omegai);


	    // compute omegaiMin1*viMin1

	    for (j=0; j<N; j++)
		pi[j] = ri[j] + beta*(pi[j] - omegai*vi[j]);

	    multVec(pi, vi);

	    double denom = realPart(dotProduct(r0hat, vi));

	    alpha = rhoi / denom;

	    std::vector<T> s(N);

	    for (j=0; j<N; j++)
		s[j] = ri[j] - alpha*vi[j];

	    if (sqrt(realPart(dotProduct(s,s))) <= tolerance * normRhs) {
		for (j=0; j<N; j++)
		    xi[j] = xi[j] + alpha*pi[j];

		ri = s;
		break;
	    }

	    std::vector<T> t;
	    multVec(s, t);

	    omegai = realPart(dotProduct(t,s)) / realPart(dotProduct(t,t));

	    rhoiMin1 = rhoi;
	    rhoi = -omegai*realPart(dotProduct(r0hat,t));

	    for (j=0; j<N; j++)
		xi[j] = xi[j] + alpha*pi[j] + omegai*s[j];

	    for (j=0; j<N; j++)
		ri[j] = s[j] - omegai*t[j];

	    normRes = sqrt(realPart(dotProduct(ri, ri)));
	    normRi = normRes;

	    if (normRes <= normRhs * tolerance)
		break;

	    ////////////////////////
	}

	*maxIter = countIter;

	residuum = ri;

	result = xi;

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


    ///
    static T dotProduct(const std::vector<T>& a, const std::vector<T>& b){
	assert(a.size()==b.size());

	T sum = T(0);
	int i;

	for (i=0; i<a.size(); i++)
	    sum += a[i]*conjug(b[i]);

	return sum;
    }

    static std::complex<T> dotProductC(const std::vector< std::complex<T> >& a,
			 const std::vector< std::complex<T> >& b){
	assert(a.size()==b.size());

	std::complex<T> sum = std::complex<T>(0);
	int i;

	for (i=0; i<a.size(); i++)
	    sum += a[i]*std::conj(b[i]);

	return sum;
    }

    static std::complex<double> dotProductCRI(const std::vector< std::complex<T> >& a,
			 const std::vector< std::complex<T> >& b){
	assert(a.size()==b.size());

	std::complex<double> sum = std::complex<double>(0);
	int i;

	for (i=0; i<a.size(); i++)
	    sum += std::complex<double>( (double)std::real(a[i])*std::real(b[i]) , (double)std::imag(a[i])*std::imag(b[i]) );

	return sum;
    }

    /// complex conjugate, if T is complex, nothing if T is not
    static T conjug(const T& z){
	return z;
    }

    static double realPart(const T& z){
	return z;
    }

};

/// @if EXCLUDETHIS

template<>
inline std::complex<float> SparseMatrix<std::complex<float> >::conjug(const std::complex<float>& z)
{
    return std::conj(z);
}

template<>
inline std::complex<double> SparseMatrix<std::complex<double> >::conjug(const std::complex<double>& z)
{
    return std::conj(z);
}

template<>
inline double SparseMatrix<std::complex<float> >::realPart(const std::complex<float>& z)
{
    return std::real(z);
}

template<>
inline double SparseMatrix<std::complex<double> >::realPart(const std::complex<double>& z)
{
    return std::real(z);
}

/// @endif

#endif
