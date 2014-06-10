//-*-C++-*-

/*
 * eigensystem_plan.h
 *      Author: peter staar
 */


#ifndef EIGENSYSTEM_PLAN_H_
#define EIGENSYSTEM_PLAN_H_



template<typename real_scalartype, matrix_form_type matrix>
class eigensystem_plan
{

};

template<typename scalartype>
class eigensystem_plan<scalartype, GENERAL>
{
 public:
  
  //  eigensystem_plan(int n);
  eigensystem_plan(int n, char jobvl = 'N', char jobvr = 'V');
  ~eigensystem_plan();

  void execute_plan();

  char JOBVL;
  char JOBVR;

  int N;
  int LDA;

  int LDVL;
  int LDVR;

  int LWORK;
  int INFO;

  scalartype*  A;
  scalartype*  W;

  scalartype*  VL;
  scalartype*  VR;

  scalartype*  WORK;
  scalartype*  RWORK;

};

extern "C" void sgeev_(char* JOBVL,  char* JOBVR, 
		       int* N, float* A,  int* LDA, 
		       float* W, 
		       float* VL, int* LDVL,
		       float* VR, int* LDVR,
		       float* WORK, int* LWORK, float* RWORK,
		       int*   INFO );

extern "C" void dgeev_(char* JOBVL,  char* JOBVR, 
		       int* N, double* A,  int* LDA, 
		       double* W, 
		       double* VL, int* LDVL,
		       double* VR, int* LDVR,
		       double* WORK, int* LWORK, double* RWORK,
		       int*   INFO );

extern "C" void cgeev_(char* JOBVL,  char* JOBVR, 
		       int* N, std::complex<float>* A,  int* LDA, 
		       std::complex<float>* W, 
		       std::complex<float>* VL, int* LDVL,
		       std::complex<float>* VR, int* LDVR,
		       std::complex<float>* WORK, int* LWORK, std::complex<float>* RWORK,
		       int*   INFO );

extern "C" void zgeev_(char* JOBVL,  char* JOBVR, 
		       int* N, std::complex<double>* A,  int* LDA, 
		       std::complex<double>* W, 
		       std::complex<double>* VL, int* LDVL,
		       std::complex<double>* VR, int* LDVR,
		       std::complex<double>* WORK, int* LWORK, std::complex<double>* RWORK,
		       int*   INFO );

template<typename scalartype>
eigensystem_plan<scalartype, GENERAL>::eigensystem_plan(int n, char jobvl, char jobvr):
  JOBVL(jobvl),
  JOBVR(jobvr),

  N(n),
  LDA(n),

  LDVL(1),
  LDVR(n),

  LWORK(4*n),
  INFO(1)
{
  A     = new scalartype[N*N];
  W     = new scalartype[N];

  if(jobvl == 'V')
    {
      LDVL  = N;
      VL    = new scalartype[N*N];
    }
  else
    {
      LDVL  = 1;
      VL    = new scalartype[1];
    }
  
  if(jobvr == 'V')
    {
      LDVR  = N;
      VR    = new scalartype[N*N];
    }
  else
    {
      LDVR  = 1;
      VR    = new scalartype[1];
    }
  
  VR    = new scalartype[N*N];

  WORK  = new scalartype[LWORK];
  RWORK = new scalartype[2*N];
}

template<typename scalartype>
eigensystem_plan<scalartype, GENERAL>::~eigensystem_plan()
{
  delete []  A     ;
  delete []  W     ;
  
  delete []  VL    ;
  delete []  VR    ;
  
  delete []  WORK  ;
  delete []  RWORK ;
}

template<typename real_scalartype>
void eigensystem_plan<real_scalartype, GENERAL>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<float, GENERAL>::execute_plan()
{
  sgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, &INFO);
  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<double, GENERAL>::execute_plan()
{
  dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, &INFO);
  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<std::complex<float>, GENERAL>::execute_plan()
{
  cgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, &INFO);
  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<std::complex<double>, GENERAL>::execute_plan()
{
  zgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, &INFO);
  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}


/* 
*******************************************
***       REAL SYMMETRIC MATRICES       ***             
*******************************************
*/

template<typename real_scalartype>
class eigensystem_plan<real_scalartype, HERMITIAN>
{

 public:
  
  //  eigensystem_plan(int n);
  eigensystem_plan(int n, char jobz = 'V', char uplo = 'U');
  ~eigensystem_plan();

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
  real_scalartype*  Matrix;
  real_scalartype*  eigenvalues;

 private:

  real_scalartype*  WORK;
};

extern "C" void ssyev_(const char* JOBZ, const char* UPLO, const int* N, float* Matrix, const int* LDA, 
		      float* eigenvalues, float* WORK, const int* LWORK, int* INFO );

extern "C" void dsyev_(const char* JOBZ, const char* UPLO, const int* N, double* Matrix, const int* LDA, 
		      double* eigenvalues, double* WORK, const int* LWORK, int* INFO );

template<typename real_scalartype>
eigensystem_plan<real_scalartype, HERMITIAN>::eigensystem_plan(int n, char jobz, char uplo):
N(n),
  JOBZ(jobz),
  UPLO(uplo),
  LDA(N)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  Matrix      = new real_scalartype[N*N];
  memset(Matrix, 0, sizeof(real_scalartype)*N*N);

  eigenvalues = new real_scalartype[N];
  
  WORK        = new real_scalartype[1];

  LWORK = -1;

  find_optimum_workspace();
}

template<typename real_scalartype>
eigensystem_plan<real_scalartype, HERMITIAN>::~eigensystem_plan()
{
  //cout << __PRETTY_FUNCTION__ << " is terminated" << endl;

  delete [] Matrix;
  delete [] eigenvalues;

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
  //cout << __PRETTY_FUNCTION__ << endl;

  ssyev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, &INFO );
  assert(INFO == 0);

  LWORK = int(WORK[0]);

  delete [] WORK;
  WORK = new float[LWORK];
}

template<>
void eigensystem_plan<double, HERMITIAN>::find_optimum_workspace()
{
  //cout << __PRETTY_FUNCTION__ << endl;

  dsyev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, &INFO );
  assert(INFO == 0);

  LWORK = int(WORK[0]);

  delete [] WORK;
  WORK = new double[LWORK];
}

template<typename real_scalartype>
void eigensystem_plan<real_scalartype, HERMITIAN>::execute_plan()
{
  assert(false);
}

template<>
void eigensystem_plan<float, HERMITIAN>::execute_plan()
{
  //cout << __PRETTY_FUNCTION__ << endl;

  ssyev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, &INFO);
  assert(INFO == 0);
}

template<>
void eigensystem_plan<double, HERMITIAN>::execute_plan()
{
  //cout << __PRETTY_FUNCTION__ << endl;

  dsyev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, &INFO);
  assert(INFO == 0);
}



/* 
**********************************************
***       COMPLEX HERMITIAN MATRICES       ***             
**********************************************
*/

template<typename real_scalartype>
class eigensystem_plan<std::complex<real_scalartype>, HERMITIAN >
{

 public:
  
  //  eigensystem_plan(int n);
  eigensystem_plan(int n, char jobz = 'V', char uplo = 'U');
  ~eigensystem_plan();

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

extern "C" void cheev_( const char* jobz, const char* uplo, const int* n, std::complex<float>* a, const int* lda,
		       float* w, std::complex<float>* work, const int* lwork, float* rwork, int* info );


extern "C" void zheev_(const char* jobz, const char* uplo, const int* n, std::complex<double>* a, const int* lda,
		      double* w, std::complex<double>* work, const int* lwork, double* rwork, int* info );


template<typename real_scalartype>
eigensystem_plan<std::complex<real_scalartype>, HERMITIAN >::eigensystem_plan(int n, char jobz, char uplo):
  N(n),
  JOBZ(jobz),
  UPLO(uplo),
  LDA(N)
{
  //cout << __PRETTY_FUNCTION__ << endl;

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
  //cout << __PRETTY_FUNCTION__ << " is terminated" << endl;

  delete [] Matrix;
  delete [] eigenvalues;

  delete [] RWORK;
  delete [] WORK;
}

template<typename real_scalartype>
void eigensystem_plan<std::complex<real_scalartype>, HERMITIAN >::find_optimum_workspace()
{
  assert(false);
}

template<>
void eigensystem_plan<std::complex<float>, HERMITIAN >::find_optimum_workspace()
{
  //cout << __PRETTY_FUNCTION__ << endl;

  cheev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, RWORK, &INFO );

  LWORK = int(real(WORK[0]));

  delete [] WORK;
  WORK = new std::complex<float>[LWORK];
}

template<>
void eigensystem_plan<std::complex<double>, HERMITIAN >::find_optimum_workspace()
{
  //cout << __PRETTY_FUNCTION__ << endl;

  zheev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, RWORK, &INFO );

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
  //cout << __PRETTY_FUNCTION__ << endl;

  cheev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, RWORK, &INFO );
}

template<>
void eigensystem_plan<std::complex<double>, HERMITIAN >::execute_plan()
{
  //cout << __PRETTY_FUNCTION__ << endl;

  zheev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, RWORK, &INFO );
}

#endif
