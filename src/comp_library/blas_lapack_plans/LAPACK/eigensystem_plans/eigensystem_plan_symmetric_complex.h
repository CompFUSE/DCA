//-*-C++-*-

#ifndef EIGENSYSTEM_PLAN_SYMMETRIC_COMPLEX_H
#define EIGENSYSTEM_PLAN_SYMMETRIC_COMPLEX_H

/*!
 *  \author  Peter Staar
 *  \version 1.0
 */
template<typename real_scalartype>
class eigensystem_plan<std::complex<real_scalartype>, HERMITIAN >
{

 public:
  
  eigensystem_plan(int n, char jobz = 'V', char uplo = 'U');
  ~eigensystem_plan();

  void set_N(int n);

  void find_optimum_workspace();
  void execute_plan();

 private:

  int N;
  char JOBZ;
  char UPLO;
  int LDA;
  int LWORK;
  int INFO;

 public:
  std::complex<real_scalartype>*  Matrix;
  real_scalartype              *  eigenvalues;

 private:

  real_scalartype*                RWORK;
  std::complex<real_scalartype>*  WORK;
};

template<typename real_scalartype>
eigensystem_plan<std::complex<real_scalartype>, HERMITIAN >::eigensystem_plan(int n, char jobz, char uplo):
  N(n),
  JOBZ(jobz),
  UPLO(uplo),
  LDA(N)
{
  Matrix      = new std::complex<real_scalartype>[N*N];
  memset(Matrix, 0, sizeof(real_scalartype)*N*N);

  eigenvalues = new real_scalartype[N];
  
  RWORK       = new real_scalartype[3*N-2]; 
  WORK        = new std::complex<real_scalartype>[1];

  LWORK = -1;

  find_optimum_workspace();
}

template<typename real_scalartype>
eigensystem_plan<std::complex<real_scalartype>, HERMITIAN >::~eigensystem_plan()
{
  delete [] Matrix;
  delete [] eigenvalues;

  delete [] RWORK;
  delete [] WORK;
}

template<typename real_scalartype>
void eigensystem_plan<std::complex<real_scalartype>, HERMITIAN>::set_N(int n)
{
  N   = n;
  LDA = n;

  Matrix      = new std::complex<real_scalartype>[N*N];
  memset(Matrix, 0, sizeof(std::complex<real_scalartype>)*N*N);

  eigenvalues = new real_scalartype[N];
  
  RWORK       = new real_scalartype[3*N-2]; 
  WORK        = new std::complex<real_scalartype>[1];

  LWORK = -1;

  find_optimum_workspace();
}

template<typename real_scalartype>
void eigensystem_plan<std::complex<real_scalartype>, HERMITIAN >::find_optimum_workspace()
{
  assert(false);
}

template<>
void eigensystem_plan<std::complex<float>, HERMITIAN >::find_optimum_workspace()
{
  LAPACK::cheev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, RWORK, &INFO );

  LWORK = int(real(WORK[0]));

  delete [] WORK;
  WORK = new std::complex<float>[LWORK];
}

template<>
void eigensystem_plan<std::complex<double>, HERMITIAN >::find_optimum_workspace()
{
  LAPACK::zheev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, RWORK, &INFO );

  LWORK = int(real(WORK[0]));

  delete [] WORK;
  WORK = new std::complex<double>[LWORK];
}

template<typename real_scalartype>
void eigensystem_plan<std::complex<real_scalartype>, HERMITIAN >::execute_plan()
{
  assert(false);
}

template<>
void eigensystem_plan<std::complex<float>, HERMITIAN >::execute_plan()
{
  LAPACK::cheev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, RWORK, &INFO );
}

template<>
void eigensystem_plan<std::complex<double>, HERMITIAN >::execute_plan()
{
  LAPACK::zheev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, RWORK, &INFO );
}

#endif
