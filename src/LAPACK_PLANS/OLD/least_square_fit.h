//-*-C++-*-

/*
 * least_square_plan.h
 *
 *      Author: peter staar
 */

#ifndef LEAST_SQUARE_PLAN_H_
#define LEAST_SQUARE_PLAN_H_


template<typename scalartype>
class least_square_plan
{

 public:

  least_square_plan(int m, int n, int nrhs, scalartype* a, scalartype* b);
  ~least_square_plan();

  int execute_plan();

 public:

  int M;
  int N;

  int NRHS;

  scalartype* A;
  int LDA;

  scalartype* B;
  int LDB;

  double* S;

  double RCOND;
  int RANK;

  scalartype* WORK;
  int LWORK;

  double* RWORK;

  int* IWORK;

  int INFO;

};

//       SUBROUTINE ZGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
//      $                   WORK, LWORK, RWORK, INFO )
// *
// *  -- LAPACK driver routine (version 3.2) --
// *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
// *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
// *     November 2006
// *
// *     .. Scalar Arguments ..
//       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
//       DOUBLE PRECISION   RCOND
// *     ..
// *     .. Array Arguments ..
//       DOUBLE PRECISION   RWORK( * ), S( * )
//       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
// *     ..
// *
// *  Purpose
// *  =======
// *
// *  ZGELSS computes the minimum norm solution to a complex linear
// *  least squares problem:
// *
// *  Minimize 2-norm(| b - A*x |).
// *
// *  using the singular value decomposition (SVD) of A. A is an M-by-N
// *  matrix which may be rank-deficient.
// *
// *  Several right hand side vectors b and solution vectors x can be
// *  handled in a single call; they are stored as the columns of the
// *  M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
// *  X.
// *
// *  The effective rank of A is determined by treating as zero those
// *  singular values which are less than RCOND times the largest singular
// *  value.
// *
// *  Arguments
// *  =========
// *
// *  M       (input) INTEGER
// *          The number of rows of the matrix A. M >= 0.
// *
// *  N       (input) INTEGER
// *          The number of columns of the matrix A. N >= 0.
// *
// *  NRHS    (input) INTEGER
// *          The number of right hand sides, i.e., the number of columns
// *          of the matrices B and X. NRHS >= 0.
// *
// *  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
// *          On entry, the M-by-N matrix A.
// *          On exit, the first min(m,n) rows of A are overwritten with
// *          its right singular vectors, stored rowwise.
// *
// *  LDA     (input) INTEGER
// *          The leading dimension of the array A. LDA >= max(1,M).
// *
// *  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
// *          On entry, the M-by-NRHS right hand side matrix B.
// *          On exit, B is overwritten by the N-by-NRHS solution matrix X.
// *          If m >= n and RANK = n, the residual sum-of-squares for
// *          the solution in the i-th column is given by the sum of
// *          squares of the modulus of elements n+1:m in that column.
// *
// *  LDB     (input) INTEGER
// *          The leading dimension of the array B.  LDB >= max(1,M,N).
// *
// *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
// *          The singular values of A in decreasing order.
// *          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
// *
// *  RCOND   (input) DOUBLE PRECISION
// *          RCOND is used to determine the effective rank of A.
// *          Singular values S(i) <= RCOND*S(1) are treated as zero.
// *          If RCOND < 0, machine precision is used instead.
// *
// *  RANK    (output) INTEGER
// *          The effective rank of A, i.e., the number of singular values
// *          which are greater than RCOND*S(1).
// *
// *  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
// *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
// *
// *  LWORK   (input) INTEGER
// *          The dimension of the array WORK. LWORK >= 1, and also:
// *          LWORK >=  2*min(M,N) + max(M,N,NRHS)
// *          For good performance, LWORK should generally be larger.
// *
// *          If LWORK = -1, then a workspace query is assumed; the routine
// *          only calculates the optimal size of the WORK array, returns
// *          this value as the first entry of the WORK array, and no error
// *          message related to LWORK is issued by XERBLA.
// *
// *  RWORK   (workspace) DOUBLE PRECISION array, dimension (5*min(M,N))
// *
// *  INFO    (output) INTEGER
// *          = 0:  successful exit
// *          < 0:  if INFO = -i, the i-th argument had an illegal value.
// *          > 0:  the algorithm for computing the SVD failed to converge;
// *                if INFO = i, i off-diagonal elements of an intermediate
// *                bidiagonal form did not converge to zero.


extern "C" void zgelss_(const int* M, const int* N, const int* NRHS, 
			std::complex<double>* A, const int* LDA, 
			std::complex<double>* B, const int* LDB, 
			double* S, 
			const double* RCOND, int* RANK,
			std::complex<double>* WORK, const int* LWORK, 
			double* RWORK, int* INFO);


template<typename scalartype>
least_square_plan<scalartype>::least_square_plan(int m, int n, int nrhs, scalartype* a, scalartype* b):
  M(m),
  N(n),
  NRHS(nrhs),

  A(a),
  LDA(m),

  B(b),
  LDB(m),

  RCOND(-1),

  LWORK(2*(2*min(M,N) + max(M,max(N,NRHS))))
{
  //cout << __FUNCTION__ << endl;

  S    = new double[min(N,M)];
  WORK = new scalartype[LWORK];
  RWORK = new double[5*(5*min(N,M))];
}


template<typename scalartype>
least_square_plan<scalartype>::~least_square_plan()
{
  //cout << __FUNCTION__ << endl;

  delete [] S;
  delete [] WORK;
  delete [] RWORK;
}

template<typename scalartype>
int least_square_plan<scalartype>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

// template<>
// void least_square_plan<float>::execute_plan()
// {
//   sgelsd_( &M, &N, &NRHS, A, &LDA, B, &LDB, S, &RCOND, &RANK,
// 	   WORK, &LWORK, IWORK, &INFO);
// }

// template<>
// void least_square_plan<double>::execute_plan()
// {
//   dgelsd_( &M, &N, &NRHS, A, &LDA, B, &LDB, S, &RCOND, &RANK,
// 	   WORK, &LWORK, IWORK, &INFO);
// }

// template<>
// void least_square_plan<std::complex<float> >::execute_plan()
// {
//   cgelsd_( &M, &N, &NRHS, A, &LDA, B, &LDB, S, &RCOND, &RANK,
// 	   WORK, &LWORK, IWORK, &INFO);
//}

template<>
int least_square_plan<std::complex<double> >::execute_plan()
{
//   cout << __FUNCTION__ << endl;
//   cout << LWORK << endl;

  zgelss_( &M, &N, &NRHS, 
	   A, &LDA, 
	   B, &LDB, 
	   S, &RCOND, &RANK,
	   WORK, &LWORK, 
	   RWORK,
	   &INFO);

  return INFO;
}

#endif
