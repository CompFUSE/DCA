//-*-C++-*-

#ifndef QR_PLAN_REAL_H
#define QR_PLAN_REAL_H

/*!
 *   \author: peter staar
 */
template<typename scalartype>
class qr_plan
{

 public:

  qr_plan(int m, int n);
  ~qr_plan();

  int execute_plan();

  void compute_Q_and_R(scalartype* Q, scalartype* R);

 public:

  int M;
  int N;

  scalartype* A;
  int LDA;

  int*    JPVT;
  scalartype* TAU;

  scalartype* WORK;
  int LWORK;

  int INFO;
};


template<typename scalartype>
qr_plan<scalartype>::qr_plan(int m, int n):
  M(m),
  N(n),

  A(NULL),
  LDA(m),

  JPVT(NULL),
  TAU(NULL),

  WORK(NULL),
  LWORK(-1),
  
  INFO(0)
{
  A = new scalartype[M*N];

  JPVT = new int[N];
  TAU  = new scalartype[std::min(M,N)];

  WORK = new scalartype[1];
  
  execute_plan();

  LWORK = WORK[0];
  
  delete [] WORK;
  WORK = new scalartype[LWORK];
}

template<typename scalartype>
qr_plan<scalartype>::~qr_plan()
{
  delete [] A;
  delete [] JPVT;
  delete [] TAU;
  delete [] WORK;
}


template<typename scalartype>
int qr_plan<scalartype>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

template<>
int qr_plan<float>::execute_plan()
{
  for(int i=0; i<N; ++i)
    JPVT[i] = 0;

  LAPACK::sgeqp3_(&M,&N, A, &LDA, JPVT, TAU, WORK, &LWORK, &INFO );
  return INFO;
}

template<>
int qr_plan<double>::execute_plan()
{
  for(int i=0; i<N; ++i)
    JPVT[i] = 0;

  LAPACK::dgeqp3_(&M,&N, A, &LDA, JPVT, TAU, WORK, &LWORK, &INFO );
  return INFO;
}

template<typename scalartype>
void qr_plan<scalartype>::compute_Q_and_R(scalartype* /*Q*/, scalartype* /*R*/)
{
  throw std::logic_error(__FUNCTION__);
}

template<>
void qr_plan<float>::compute_Q_and_R(float* Q, float* R)
{
  for(int i=0; i<std::min(M,N); ++i)
    for(int j=0; j<std::min(M,N); ++j)
      if(i>j)
	R[i+j*std::min(M,N)] = 0;
      else
	R[i+j*std::min(M,N)] = A[i+j*M];

  memcpy(Q, A, M*N*sizeof(float));

  LAPACK::sorgqr_(&M, &N, &N, Q, &LDA, TAU, WORK, &LWORK, &INFO);
}

template<>
void qr_plan<double>::compute_Q_and_R(double* Q, double* R)
{
  for(int i=0; i<std::min(M,N); ++i)
    for(int j=0; j<std::min(M,N); ++j)
      if(i>j)
	R[i+j*std::min(M,N)] = 0;
      else
	R[i+j*std::min(M,N)] = A[i+j*M];

  memcpy(Q, A, M*N*sizeof(double));

  LAPACK::dorgqr_(&M, &N, &N, Q, &LDA, TAU, WORK, &LWORK, &INFO);
}

#endif
