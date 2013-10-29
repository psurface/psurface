#include "config.h"

#include <fenv.h>
#include <cassert>

#include <algorithm>
#include <complex>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <vector>

#include "SparseMatrix.h"

using namespace std;
using namespace psurface;

template<typename ctype>
void print_vector(const vector<ctype>& vec) {
  for (size_t i = 0; i < vec.size(); ++i)
    cout << vec[i] << endl;
  cout << endl;
}

// !!! Test reasonable ?
template<typename ctype, typename rhstype>
void check_result(const vector<rhstype>& result,
                  const vector<rhstype>& trueResult,
                  const double& tolerance, char const * const message) {
  const size_t n = trueResult.size();
  vector<rhstype> res(n);

  for (size_t i = 0; i < n; ++i)
    res[i] = result[i] - trueResult[i];

  if (sqrt(std::real(SparseMatrix<ctype>::dotProductC(res, res))) > tolerance)
    throw runtime_error(message);
}


template<typename ctype, typename rhstype>
void test(const double tolerance, const int n, const int maxIter) {
  //// setup for n-dimensional tests
  // id matrix test
  int countIter;
  vector<rhstype> residue(n);
  vector<rhstype> result(n);

  // setup first rhs
  vector<rhstype> b(n);
  for (size_t i = 0; i < n; ++i)
    b[i] = rhstype(i+1, i+1+n);

  SparseMatrix<ctype> id(n);
  for (size_t i = 0; i < n; ++i)
    id.setEntry(i, i, 1);

  countIter = maxIter;
  id.BiCGSTABC(b, result, residue, &countIter, tolerance);

  print_vector(result);
  check_result<ctype, rhstype>(result, b, tolerance, "id matrix check 1 failed");

  // setup another rhs
  for (size_t i = 0; i < n; ++i)
    b[i] = rhstype(0, 0);

  countIter = maxIter;
  id.BiCGSTABC(b, result, residue, &countIter, tolerance);

  print_vector(result);
  check_result<ctype, rhstype>(result, b, tolerance, "id matrix check 2 failed");


  // ladder matrix test
  SparseMatrix<ctype> ladder(n);
  for (size_t i = 0; i < n; ++i)
    ladder.setEntry(i, i, i+1);
  for (size_t i = 0; i < n; ++i)
    b[i] = rhstype(i+1, (i+1)*2);

  // obtain solutions
  countIter = maxIter;
  ladder.BiCGSTABC(b, result, residue, &countIter, tolerance);
  vector<rhstype> trueLadder(n);
  for (size_t i = 0; i < n; ++i)
    trueLadder[i] = rhstype(1, 2);

  // check results
  print_vector(result);
  check_result<ctype, rhstype>(result, b, tolerance, "ladder matrix check failed");

  //// setup specific tests
  vector<rhstype> residue1(3);
  vector<rhstype> result1(3);
  vector<rhstype> b1(3);
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

  b1[0] = rhstype(0, 5);
  b1[1] = rhstype(4, 0);
  b1[2] = rhstype(3, 6);

  // solutions
  countIter = maxIter;
  test1.BiCGSTABC(b1, result1, residue1, &countIter, tolerance);
  vector<rhstype> trueTest1(3);
  trueTest1[0] = rhstype(-7, -1);
  trueTest1[1] = rhstype(2, 0);
  trueTest1[2] = rhstype(1, 2);

  // check results
  print_vector(result1);
  check_result<ctype, rhstype>(result1, trueTest1, tolerance, "specific matrix check 1 failed");

  // second specific matrix
  vector<rhstype> residue2(5);
  vector<rhstype> result2(5);
  vector<rhstype> b2(5);
  SparseMatrix<ctype> test2(5);

  test2.setEntry(0, 0, 1);
  test2.setEntry(1, 1, 1);
  test2.setEntry(2, 2, 1);
  test2.setEntry(3, 3, 1);
  test2.setEntry(4, 4, 1);
  test2.setEntry(1, 2, -0.473423);
  test2.setEntry(1, 0, -0.166667);
  test2.setEntry(1, 3, -0.193244);
  test2.setEntry(1, 4, -0.166667);

  b2[0] = rhstype(-.5, 0.866025);
  b2[1] = rhstype(0, 0);
  b2[2] = rhstype(1, 0);
  b2[3] = rhstype(-1, -1.01577e-07);
  b2[4] = rhstype(-0.5, -0.866026);

  // solutions
  countIter = maxIter;
  test2.BiCGSTABC(b2, result2, residue2, &countIter, tolerance);
  vector<rhstype> trueTest2(5);
  trueTest2[0] = rhstype(-.5, 0.866);
  trueTest2[1] = rhstype(0.1135, 0);
  trueTest2[2] = rhstype(1, 0);
  trueTest2[3] = rhstype(-1, 0);
  trueTest2[4] = rhstype(-.5, -0.866);

  // check results
  print_vector(result2);
  check_result<ctype, rhstype>(result2, trueTest2, 1e-3, "specific matrix check 2 failed");
}

int main (int argc, char* argv[]) {

  feenableexcept(FE_INVALID);

  const int n = 1, maxIter = 3000;
  const double tolerance = 1e-6;

  try {
    test<float, complex<float> >(tolerance, n, maxIter);
    test<double, complex<double> >(tolerance, n, maxIter);
    test<long double, complex<long double> >(tolerance, n, maxIter);
  } catch (const exception& e) {
    cout << e.what() << endl;

    return 1;
  }

  return 0;
}
