//-*-C++-*-

#ifndef SINGULAR_VALUE_DECOMPOSITION_PLAN_H
#define SINGULAR_VALUE_DECOMPOSITION_PLAN_H

/*!
 *  \class   singular_value_decomposition_plan
 *  \author  Peter Staar
 *  \version 1.0
 */
template<typename real_scalartype, matrix_form_type matrix>
class singular_value_decomposition_plan
{};

template<typename scalartype>
class singular_value_decomposition_plan<std::complex<scalartype>, GENERAL>
{
 public:
  
  singular_value_decomposition_plan(int n, char jobu = 'A', char jobvt = 'A');
  singular_value_decomposition_plan(int m, int n, char jobu = 'N', char jobvt = 'N');

  ~singular_value_decomposition_plan();

  void reset(int n, char jobu, char jobvt);

  void execute_plan();

  char JOBU;
  char JOBVT;

  int M;
  int N;

  std::complex<scalartype>* A;
  int                       LDA;

  scalartype* S;

  std::complex<scalartype>* U;
  int                       LDU;

  std::complex<scalartype>* VT;
  int                       LDVT;

  int LWORK;
  std::complex<scalartype>*  WORK;
  scalartype*               RWORK;

  int INFO;
};

template<typename scalartype>
singular_value_decomposition_plan<std::complex<scalartype>, GENERAL>::singular_value_decomposition_plan(int n, char jobu, char jobvt):
  JOBU(jobu),
  JOBVT(jobvt),

  M(n),
  N(n),

  A(NULL),
  LDA(n),

  S(NULL),
  
  U(NULL),
  LDU(n),

  VT(NULL),
  LDVT(n),

  LWORK(256*(2*min(M,N)+max(M,N))),
  WORK(NULL),
  RWORK(NULL),

  INFO(0)

{
  A  = new std::complex<scalartype>[N*N];
  U  = new std::complex<scalartype>[N*N];
  S  = new scalartype[N];
  VT = new std::complex<scalartype>[N*N];

  WORK  = new std::complex<scalartype>[LWORK];
  RWORK = new scalartype[5*min(M,N)];
}

template<typename scalartype>
singular_value_decomposition_plan<std::complex<scalartype>, GENERAL>::singular_value_decomposition_plan(int m, int n, char jobu, char jobvt):
  JOBU(jobu),
  JOBVT(jobvt),

  M(m),
  N(n),

  A(NULL),
  LDA(m),

  S(NULL),
  
  U(NULL),
  LDU(M),

  VT(NULL),
  LDVT(n),

  LWORK(256*(2*min(M,N)+max(M,N))),
  WORK(NULL),
  RWORK(NULL),

  INFO(0)

{
  A  = new std::complex<scalartype>[M*N];
  U  = new std::complex<scalartype>[M*M];
  S  = new scalartype[std::min(M,N)];
  VT = new std::complex<scalartype>[N*N];

  WORK  = new std::complex<scalartype>[LWORK];
  RWORK = new scalartype[5*min(M,N)];
}

template<typename scalartype>
singular_value_decomposition_plan<std::complex<scalartype>, GENERAL>::~singular_value_decomposition_plan()
{
  delete [] A;
  delete [] U;
  delete [] S;
  delete [] VT;

  delete [] WORK;
  delete [] RWORK;
}

template<typename real_scalartype>
void singular_value_decomposition_plan<std::complex<real_scalartype>, GENERAL>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

template<>
void singular_value_decomposition_plan<std::complex<double>, GENERAL>::execute_plan()
{
  LAPACK::zgesvd_( &JOBU, &JOBVT, &M, &N, A, &LDA, S, U, &LDU, VT, &LDVT, WORK, &LWORK, RWORK, &INFO);

  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

#endif
