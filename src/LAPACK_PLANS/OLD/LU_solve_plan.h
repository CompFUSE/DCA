//-*-C++-*-

/*
 * LU_solve_plan.h
 *
 *      Author: peter staar
 */

#ifndef LU_SOLVE_PLAN_H_
#define LU_SOLVE_PLAN_H_

template<typename scalartype>
class LU_solve_plan
{

 public:

  LU_solve_plan(int n, int nrhs, int lda, int ldb);
  ~LU_solve_plan();

  void execute_plan();

  static void check();
  
public:

  char TRANS;
  int N;
  int NRHS;
  int LDA;
  int LDB;
  int INFO;

  scalartype*  Matrix_A;
  scalartype*  Matrix_B;
  int*             IPIV;

};

// extern "C" void sgetrs_(const char* TRANS, const int *N, const int *NRHS, float *A, const int *LDA, int *IPIV, float *B, const int *LDB, int *INFO);
// extern "C" void dgetrs_(const char* TRANS, const int *N, const int *NRHS, double *A, const int *LDA, int *IPIV, double *B, const int *LDB, int *INFO);
// extern "C" void cgetrs_(const char* TRANS, const int *N, const int *NRHS, std::complex<float> *A, const int *LDA, int *IPIV, std::complex<float> *B, const int *LDB, int *INFO);
// extern "C" void zgetrs_(const char* TRANS, const int *N, const int *NRHS, std::complex<double> *A, const int *LDA, int *IPIV, std::complex<double> *B, const int *LDB, int *INFO);

template<typename scalartype>
LU_solve_plan<scalartype>::LU_solve_plan(int n, int nrhs, int lda, int ldb):
  TRANS('N'),
  N(n),
  NRHS(nrhs),
  LDA(lda),
  LDB(ldb)
{
  //cout << __FUNCTION__ << endl;

  IPIV = new int[N];

  for(int i=0; i<N; i++)
    IPIV[i] = i+1;
}
  
template<typename scalartype>
LU_solve_plan<scalartype>::~LU_solve_plan()
{
  delete [] IPIV;
}
  
template<typename scalartype>
void LU_solve_plan<scalartype>::execute_plan()
{
  throw std::logic_error(__PRETTY_FUNCTION__);
}

template<>
void LU_solve_plan<float>::execute_plan()
{
  sgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);

  if(INFO!=0)
    throw std::logic_error(__PRETTY_FUNCTION__); 
}

template<>
void LU_solve_plan<double>::execute_plan()
{

  dgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);

  if(INFO!=0)
    throw std::logic_error(__PRETTY_FUNCTION__);
}

template<>
void LU_solve_plan<std::complex<float> >::execute_plan()
{
  cgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);
  
  if(INFO!=0)
    throw std::logic_error(__PRETTY_FUNCTION__);
}

template<>
void LU_solve_plan<std::complex<double> >::execute_plan()
{
  zgetrs_(&TRANS, &N, &NRHS, Matrix_A, &LDA, IPIV, Matrix_B, &LDB, &INFO);
  
  if(INFO!=0)
    throw std::logic_error(__PRETTY_FUNCTION__); 
}


template<>
void LU_solve_plan<double>::check()
{
  cout << __FUNCTION__ << endl;
  /*double A[9*9] = {1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6,  
    1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4, 
    1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  
    1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  
    1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  
    2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  
    2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  
    2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  
    2.6,  2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0};*/

  double A_LU[9*9] = {2.6,   2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  
		      0.4,   0.3,  0.6,  0.8,  1.1,  1.4,  1.7,  1.9,  2.2,  
		      0.5,  -0.4,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  2.8,  
		      0.5,  -0.3,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  2.4,  
		      0.6,  -0.3,  0.0,  0.0,  0.4,  0.8,  1.2,  1.6,  2.0,  
		      0.7,  -0.2,  0.0,  0.0,  0.0,  0.4,  0.8,  1.2,  1.6,  
		      0.8,  -0.2,  0.0,  0.0,  0.0,  0.0,  0.4,  0.8,  1.2,  
		      0.8,  -0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.4,  0.8,  
		      0.9,  -0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.4};

  double BX[9*5] = {93.0,  186.0,  279.0,  372.0,  465.0,  
		    84.4,  168.8,  253.2,  337.6,  422.0,  
		    76.6,  153.2,  229.8,  306.4,  383.0,  
		    70.0,  140.0,  210.0,  280.0,  350.0,  
		    65.0,  130.0,  195.0,  260.0,  325.0,  
		    62.0,  124.0,  186.0,  248.0,  310.0,  
		    61.4,  122.8,  184.2,  245.6,  307.0,  
		    63.6,  127.2,  190.8,  254.4,  318.0,  
		    69.0,  138.0,  207.0,  276.0,  345.0 };


  double BX_output[9*5] = {1.0,   2.0,   3.0,   4.0,   5.0,  
			   2.0,   4.0,   6.0,   8.0,  10.0,  
			   3.0,   6.0,   9.0,  12.0,  15.0,  
			   4.0,   8.0,  12.0,  16.0,  20.0,  
			   5.0,  10.0,  15.0,  20.0,  25.0,  
			   6.0,  12.0,  18.0,  24.0,  30.0,  
			   7.0,  14.0,  21.0,  28.0,  35.0,  
			   8.0,  16.0,  24.0,  32.0,  40.0,  
			   9.0,  18.0,  27.0,  36.0,  45.0};

  LU_solve_plan<double> pln(9, 5, 9, 9);

  pln.TRANS = 'T';
  pln.Matrix_A = A_LU;
  pln.Matrix_B = BX;

  pln.execute_plan();

  cout << scientific;
  for(int i=0; i<9; i++){
    for(int j=0; j<5; j++){
      cout << "\t" << BX[i+9*j];
    }
    cout << endl;
  }

  cout << endl;

  for(int i=0; i<9; i++){
    for(int j=0; j<5; j++){
      cout << "\t" << BX_output[i+9*j];
    }
    cout << endl;
  }

  cout << endl;

}
#endif
