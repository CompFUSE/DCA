//-*-C++-*-

#ifndef EIGENSYSTEM_PLAN_SYMMETRIC_REAL_H
#define EIGENSYSTEM_PLAN_SYMMETRIC_REAL_H

/*!
 *  \author  Peter Staar
 *  \version 1.0
 */
template<typename real_scalartype>
class eigensystem_plan<real_scalartype, HERMITIAN>
{

 public:
  
  eigensystem_plan(int n, char jobz = 'V', char uplo = 'U');
  ~eigensystem_plan();

  void find_optimum_workspace();

  void execute_plan();

  int execute_plan(real_scalartype VL, real_scalartype VU);

 private:

  int N;
  char JOBZ;
  char UPLO;
  int LDA;
  int LWORK;
  int INFO;

 public:

  real_scalartype*  A;
  real_scalartype*  W;

 private:

  real_scalartype*  WORK;
};

template<typename real_scalartype>
eigensystem_plan<real_scalartype, HERMITIAN>::eigensystem_plan(int n, char jobz, char uplo):
  N(n),
  JOBZ(jobz),
  UPLO(uplo),
  LDA(N)
{
  A      = new real_scalartype[N*N];
  memset(A, 0, sizeof(real_scalartype)*N*N);

  W = new real_scalartype[N];
  
  WORK        = new real_scalartype[1];

  LWORK = -1;

  find_optimum_workspace();
}

template<typename real_scalartype>
eigensystem_plan<real_scalartype, HERMITIAN>::~eigensystem_plan()
{
  delete [] A;
  delete [] W;
  delete [] WORK;
}

template<typename real_scalartype>
void eigensystem_plan<real_scalartype, HERMITIAN>::find_optimum_workspace()
{
  assert(false);
}

template<>
void eigensystem_plan<float, HERMITIAN>::find_optimum_workspace()
{
  LAPACK::ssyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO );
  assert(INFO == 0);

  LWORK = int(WORK[0]);

  delete [] WORK;
  WORK = new float[LWORK];
}

template<>
void eigensystem_plan<double, HERMITIAN>::find_optimum_workspace()
{
  LAPACK::dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO );
  assert(INFO == 0);

  LWORK = int(WORK[0]);

  delete [] WORK;
  WORK = new double[LWORK];
}

template<typename real_scalartype>
void eigensystem_plan<real_scalartype, HERMITIAN>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<float, HERMITIAN>::execute_plan()
{
  LAPACK::ssyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);
  if(INFO != 0) throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<double, HERMITIAN>::execute_plan()
{
  LAPACK::dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);
  if(INFO != 0) throw std::logic_error(__FUNCTION__);
}

template<typename real_scalartype>
int eigensystem_plan<real_scalartype, HERMITIAN>::execute_plan(real_scalartype VL, real_scalartype VU)
{
  throw std::logic_error(__FUNCTION__);
  return -1;
}

template<>
int eigensystem_plan<float, HERMITIAN>::execute_plan(float VL, float VU)
{
  char RANGE = 'V';
  
  int IL=0, IU=0;

  char tmp='S';
  float ABSTOL=2*LAPACK::slamch_(&tmp);

  int M=0;

  int LDZ   = N;
  float* Z = new float[N*N];

  int* IWORK = new int[5*N];
  int* IFAIL = new int[N];


  {
    LWORK = -1;
    LAPACK::ssyevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU,
		    &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, IWORK,
		    IFAIL, &INFO);
    
    LWORK = int(WORK[0]);

    delete [] WORK;
    WORK = new float[LWORK];
  }


  LAPACK::ssyevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU,
		  &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, IWORK,
		  IFAIL, &INFO);

  if(INFO != 0) throw std::logic_error(__FUNCTION__);

  for(int j=0; j<M; ++j)
    for(int i=0; i<N; ++i)
      A[i+N*j] = Z[i+N*j];

  delete [] Z;
  delete [] IWORK;
  delete [] IFAIL;

  return M;
}

template<>
int eigensystem_plan<double, HERMITIAN>::execute_plan(double VL, double VU)
{
  char RANGE = 'V';
  
  int IL=0, IU=0;

  char tmp='S';
  double ABSTOL=2*LAPACK::dlamch_(&tmp);

  int M=0;

  int LDZ   = N;
  double* Z = new double[N*N];

  int* IWORK = new int[5*N];
  int* IFAIL = new int[N];


  {
    LWORK = -1;
    LAPACK::dsyevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU,
		    &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, IWORK,
		    IFAIL, &INFO);
    
    LWORK = int(WORK[0]);

    delete [] WORK;
    WORK = new double[LWORK];
  }


  LAPACK::dsyevx_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU,
		  &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, IWORK,
		  IFAIL, &INFO);

  if(INFO != 0) throw std::logic_error(__FUNCTION__);

  for(int j=0; j<M; ++j)
    for(int i=0; i<N; ++i)
      A[i+N*j] = Z[i+N*j];

  delete [] Z;
  delete [] IWORK;
  delete [] IFAIL;

  return M;
}

#endif
