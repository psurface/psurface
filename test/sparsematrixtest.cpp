#include "config.h"

#include "fenv.h"
#include "SparseMatrix.h"

using namespace std;
using namespace psurface;


template<typename ctype>
void check_result(const Vector<ctype>& result,
                  const Vector<ctype>& trueResult,
                  const ctype& tolerance, char const * const message) {
  const size_t n = trueResult.size();
  Vector<ctype> res(n);

  if ((result - trueResult).length() > tolerance)
    throw runtime_error(message);
}

template<typename ctype>
void check_result_by_mult(const SparseMatrix<ctype>& matrix,
                          const Vector<ctype>& result,
                          const Vector<ctype>& b,
                          const ctype& tolerance, char const * const message) {
  if ((matrix.multVec(result) - b).length() > tolerance)
    throw runtime_error(message);
}


template<typename ctype>
void test(const ctype tolerance, const int n, const int maxIter) {
  //// setup for n-dimensional tests
  // id matrix test
  Vector<ctype> residue(n);
  Vector<ctype> result(n);

  // setup first rhs
  Vector<ctype> b(n);
  for (size_t i = 0; i < n; ++i)
    b[i] = StaticVector<ctype, 2>(i+1, i+1+n);

  SparseMatrix<ctype> id(n);
  for (size_t i = 0; i < n; ++i)
    id.setEntry(i, i, 1);

  id.BiCGSTAB(b, result, residue, maxIter, tolerance);
  check_result<ctype>(result, b, tolerance, "id matrix check 1 failed");

  // setup another rhs
  for (size_t i = 0; i < n; ++i)
    b[i] = StaticVector<ctype, 2>(0, 0);

  id.BiCGSTAB(b, result, residue, maxIter, tolerance);
  check_result<ctype>(result, b, tolerance, "id matrix check 2 failed");


  // ladder matrix test
  SparseMatrix<ctype> ladder(n);
  for (size_t i = 0; i < n; ++i)
    ladder.setEntry(i, i, i+1);
  for (size_t i = 0; i < n; ++i)
    b[i] = StaticVector<ctype, 2>(i+1, (i+1)*2);

  // obtain solutions
  ladder.BiCGSTAB(b, result, residue, maxIter, tolerance);
  Vector<ctype> trueLadder(n);
  for (size_t i = 0; i < n; ++i)
    trueLadder[i] = StaticVector<ctype, 2>(1, 2);

  // check results
  check_result<ctype>(result, b, tolerance, "ladder matrix check failed");

  ///////////////////////////////////////////////////////////////////////
  //// setup specific tests
  ///////////////////////////////////////////////////////////////////////

  Vector<ctype> residue1(3);
  Vector<ctype> result1(3);
  Vector<ctype> b1(3);
  SparseMatrix<ctype> test1(3);
  /*
    1 2 3             0 5               -7 -1
    0 2 0  (x, y)  =  4 0 ==> (x, y) = ( 2, 0 )
    0 0 3             3 6                1  2
   */
  test1.setEntry(0, 0, 1);
  test1.setEntry(0, 1, 2);
  test1.setEntry(0, 2, 3);
  test1.setEntry(1, 1, 2);
  test1.setEntry(2, 2, 3);

  b1[0] = StaticVector<ctype, 2>(0, 5);
  b1[1] = StaticVector<ctype, 2>(4, 0);
  b1[2] = StaticVector<ctype, 2>(3, 6);

  // solutions
  test1.BiCGSTAB(b1, result1, residue1, maxIter, tolerance);
  Vector<ctype> trueTest1(3);
  trueTest1[0] = StaticVector<ctype, 2>(-7, -1);
  trueTest1[1] = StaticVector<ctype, 2>(2, 0);
  trueTest1[2] = StaticVector<ctype, 2>(1, 2);

  // check results
  check_result<ctype>(result1, trueTest1, tolerance, "specific matrix check 1 failed");

  // second specific matrix
  Vector<ctype> residue2(6);
  Vector<ctype> result2(6);
  Vector<ctype> b2(6);
  SparseMatrix<ctype> test2(6);

  test2.setEntry(0, 0, 1);
  test2.setEntry(1, 1, 1);
  test2.setEntry(2, 2, 1);
  test2.setEntry(3, 3, 1);
  test2.setEntry(4, 4, 1);
  test2.setEntry(5, 5, 1);
  test2.setEntry(1, 2, -0.0666667);
  test2.setEntry(1, 3, -0.312071);
  test2.setEntry(1, 4, -0.312071);
  test2.setEntry(1, 0, -0.154595);
  test2.setEntry(1, 5, -0.154595);

  b2[0] = StaticVector<ctype, 2>(0.5,0.866025);
  b2[1] = StaticVector<ctype, 2>(0,0);
  b2[2] = StaticVector<ctype, 2>(1,0);
  b2[3] = StaticVector<ctype, 2>(-0.5,0.866025);
  b2[4] = StaticVector<ctype, 2>(-0.5,-0.866025);
  b2[5] = StaticVector<ctype, 2>(0.5,-0.866025);

  // solutions
  test2.BiCGSTAB(b2, result2, residue2, maxIter, tolerance);

  // check results
  check_result_by_mult(test2, result2, b2, tolerance, "specific matrix check 2 failed");

  // third specific matrix -- the identity!
  Vector<ctype> residual3(2);
  Vector<ctype> result3(2);
  Vector<ctype> b3(2);
  SparseMatrix<ctype> test3(2);

  test3.setEntry(0, 0, 1);
  test3.setEntry(1, 1, 1);

  b3[0] = StaticVector<ctype, 2>(0,0);
  b3[1] = StaticVector<ctype, 2>(0,0);

  // solutions
  test3.BiCGSTAB(b3, result3, residual3, maxIter, tolerance);

  // check results
  check_result_by_mult(test3, result3, b3, tolerance, "specific matrix check 3 failed");
}

int main (int argc, char* argv[]) {

  feenableexcept(FE_INVALID);

  const int n = 1, maxIter = 3000;
  const double tolerance = 1e-3;

  try {
    test<float>(tolerance, n, maxIter);
    test<double>(tolerance, n, maxIter);
    test<long double>(tolerance, n, maxIter);
  } catch (const exception& e) {
    cout << e.what() << endl;

    return 1;
  }

  return 0;
}
