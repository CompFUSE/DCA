//-*-C++-*-

#ifndef EIGENSYSTEM_PLAN_H
#define EIGENSYSTEM_PLAN_H

/*!
 *  \author  Peter Staar
 *  \version 1.0
 */
template<typename scalartype>
class eigensystem_plan<scalartype, GENERAL>
{
 public:
  
  eigensystem_plan(int n, char jobvl = 'N', char jobvr = 'V');
  ~eigensystem_plan();

  void reset(int n, char jobvl = 'N', char jobvr = 'V');

  void execute_plan();

public:

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
  A  = new scalartype[N*N];

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
  LAPACK::cgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &real(RWORK[0]), &INFO);

  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

template<>
void eigensystem_plan<std::complex<double>, GENERAL>::execute_plan()
{
  LAPACK::zgeev_(&JOBVL, &JOBVR, &N, A, &LDA, W, VL, &LDVL, VR, &LDVR, WORK, &LWORK, &real(RWORK[0]), &INFO);

  if(INFO != 0)
    throw std::logic_error(__FUNCTION__);
}

#endif
