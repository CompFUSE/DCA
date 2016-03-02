//-*-C++-*-

#ifndef EIGENSYSTEM_PLAN_H_
#define EIGENSYSTEM_PLAN_H_

/*!
 *  \class   eigensystem_plan
 *  \author  Peter Staar
 *  \version 1.0
 */
/*
template<typename real_scalartype, matrix_form_type matrix>
class eigensystem_plan
{};

template<typename scalartype>
class eigensystem_plan<scalartype, GENERAL>
{
 public:
  
  eigensystem_plan(int n, char jobvl = 'N', char jobvr = 'V');
  ~eigensystem_plan();

  void reset(int n, char jobvl = 'N', char jobvr = 'V');
//   void set_N(int n);

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

  scalartype*  WR;
  scalartype*  WI;

  scalartype*  VL;
  scalartype*  VR;

  scalartype*  WORK;
  scalartype*  RWORK;

};

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

  W  = new scalartype[N];
  WR = new scalartype[N];
  WI = new scalartype[N];

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
  
  //VR    = new scalartype[N*N];

  WORK  = new scalartype[LWORK];
  RWORK = new scalartype[2*N];
}

template<typename scalartype>
eigensystem_plan<scalartype, GENERAL>::~eigensystem_plan()
{
  delete []  A     ;
  delete []  W     ;
  delete []  WR     ;
  delete []  WI     ;
  
  delete []  VL    ;
  delete []  VR    ;
  
  delete []  WORK  ;
  delete []  RWORK ;
}

template<typename scalartype>
void eigensystem_plan<scalartype, GENERAL>::reset(int n, char jobvl, char jobvr)
{
  {
    delete []  A     ;
    delete []  W     ;
    delete []  WR     ;
    delete []  WI     ;
  
    
    delete []  VL    ;
    delete []  VR    ;
    
    delete []  WORK  ;
    delete []  RWORK ;
  }

  {
    JOBVL = jobvl;
    JOBVR = jobvr;

    N   = n;
    LDA = n;

    LDVL = 1;
    LDVR = n;

    LWORK = 4*n;

    A     = new scalartype[N*N];
    W     = new scalartype[N];
    
    WR = new scalartype[N];
    WI = new scalartype[N];

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
  
    //VR    = new scalartype[N*N];

    WORK  = new scalartype[LWORK];
    RWORK = new scalartype[2*N];
  }
}

template<typename real_scalartype>
void eigensystem_plan<real_scalartype, GENERAL>::execute_plan()
{
  throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<float, GENERAL>::execute_plan()
{
  LAPACK::sgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);
  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<double, GENERAL>::execute_plan()
{
  LAPACK::dgeev_(&JOBVL, &JOBVR, &N, A, &LDA, WR, WI, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &INFO);
  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<std::complex<float>, GENERAL>::execute_plan()
{
  LAPACK::cgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, &INFO);
  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<std::complex<double>, GENERAL>::execute_plan()
{
  LAPACK::zgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, RWORK, &INFO);
  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}
*/





/* 
*******************************************
***       REAL SYMMETRIC MATRICES       ***             
*******************************************
*/
/*
template<typename real_scalartype>
class eigensystem_plan<real_scalartype, HERMITIAN>
{

 public:
  
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

template<typename real_scalartype>
eigensystem_plan<real_scalartype, HERMITIAN>::eigensystem_plan(int n, char jobz, char uplo):
  N(n),
  JOBZ(jobz),
  UPLO(uplo),
  LDA(N)
{
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
  LAPACK::ssyev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, &INFO );
  assert(INFO == 0);

  LWORK = int(WORK[0]);

  delete [] WORK;
  WORK = new float[LWORK];
}

template<>
void eigensystem_plan<double, HERMITIAN>::find_optimum_workspace()
{
  LAPACK::dsyev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, &INFO );
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
  LAPACK::ssyev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, &INFO);
  assert(INFO == 0);
}

template<>
void eigensystem_plan<double, HERMITIAN>::execute_plan()
{
  LAPACK::dsyev_(&JOBZ, &UPLO, &N, Matrix, &LDA, eigenvalues, WORK, &LWORK, &INFO);
  assert(INFO == 0);
}
*/


/* 
**********************************************
***       COMPLEX HERMITIAN MATRICES       ***             
**********************************************
*/

/*
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
  N = n;
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
*/

#endif
