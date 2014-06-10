//-*-C++-*-

/*
 * invert_plan.h
 *
 *      Author: peter staar
 */

#ifndef INVERT_PLAN_H_
#define INVERT_PLAN_H_

/* 
*******************************************
***       REAL MATRICES                 ***             
*******************************************
*/

template<typename real_scalartype>
class invert_plan
{

 public:

  invert_plan(int n);
  invert_plan(int n, int lda);
  ~invert_plan();

  void execute_plan();
//   void check_plan();

 private:

  void reset_inverted_matrix();

public:

  int N;
  int LDA;
  int LWORK;
  int INFO;

  real_scalartype*  Matrix;
  real_scalartype*  inverted_matrix;

private:

  int*              IPIV;

};

// extern "C" void sgesv_(const int *N, const int *NRHS, float *A, const int *LDA, int *IPIV, float *B, const int *LDB, int *INFO);
// extern "C" void dgesv_(const int *N, const int *NRHS, double *A, const int *LDA, int *IPIV, double *B, const int *LDB, int *INFO);

template<typename real_scalartype>
invert_plan<real_scalartype>::invert_plan(int n):
  N(n),
  LDA(N)
{
  Matrix           = new real_scalartype[N*N];
  memset(Matrix, 0, sizeof(real_scalartype)*N*N);

  inverted_matrix  = new real_scalartype[N*N];
  memset(inverted_matrix, 0, sizeof(real_scalartype)*N*N);

  for(int i=0; i<N; i++)
    inverted_matrix[i+N*i] = real_scalartype(1);
  
  IPIV = new int[N];
}

template<typename real_scalartype>
invert_plan<real_scalartype>::invert_plan(int n, int lda):
  N(n),
  LDA(lda)
{
  Matrix      = new real_scalartype[LDA*LDA];
  memset(Matrix, 0, sizeof(real_scalartype)*LDA*LDA);

  inverted_matrix      = new real_scalartype[LDA*LDA];
  memset(inverted_matrix, 0, sizeof(real_scalartype)*LDA*LDA);

  for(int i=0; i<LDA; i++)
    inverted_matrix[i+LDA*i] = real_scalartype(1);
  
  IPIV = new int[N];
}

template<typename real_scalartype>
invert_plan<real_scalartype>::~invert_plan()
{
  delete [] Matrix;
  delete [] inverted_matrix;

  delete [] IPIV;
}

template<typename real_scalartype>
void invert_plan<real_scalartype>::reset_inverted_matrix()
{
  int lda = std::max(N,LDA);

  memset(inverted_matrix, 0, sizeof(real_scalartype)*lda*lda);
  
  for(int i=0; i<lda; i++)
    inverted_matrix[i+lda*i] = real_scalartype(1.);
}

template<typename real_scalartype>
void invert_plan<real_scalartype>::execute_plan()
{
  throw std::logic_error(__PRETTY_FUNCTION__);
  assert(false);
}

template<>
void invert_plan<float>::execute_plan()
{
  reset_inverted_matrix();
  LAPACK::sgesv_(&N, &N, Matrix, &LDA, IPIV, inverted_matrix, &LDA, &INFO);

  if(INFO != 0) 
    throw std::logic_error(__FUNCTION__);
}

template<>
void invert_plan<double>::execute_plan()
{
  reset_inverted_matrix();
  LAPACK::dgesv_(&N, &N, Matrix, &LDA, IPIV, inverted_matrix, &LDA, &INFO);

  if(INFO != 0) 
    throw std::logic_error(__FUNCTION__);
}

// template<typename real_scalartype>
// void invert_plan<real_scalartype>::check_plan()
// {
//   cout << endl << endl <<  __PRETTY_FUNCTION__ << endl << endl;

//   cout << scientific;
//   cout.precision(6);
//   cout.width(6);

//   int lda = std::max(N,LDA);

//   for(int i=0; i<N; i++){
//       for(int j=0; j<N; j++){
	
// 	cout << "\t" << Matrix[i+j*lda];
//       }
//       cout << endl;
//   }
//   cout << endl;

//   for(int i=0; i<N; i++){
//       for(int j=0; j<N; j++){
	
// 	cout << "\t" << inverted_matrix[i+j*lda];
//       }
//       cout << endl;
//   }
//   cout << endl;


//   // makes no sense to check back --> on exit we have Matrix = P*L*U !!!!

//   /*for(int i=0; i<N; i++){
//       for(int j=0; j<N; j++){
// 	real_scalartype c(0);
// 	for(int l=0; l<N; l++)
// 	  c += Matrix[i+l*lda]*inverted_matrix[l+j*lda];

// 	//i==j? c-=1 : c;

// 	cout << "\t" << c;
// 	//assert(abs(c)<1.e-6);
//       }
//       cout << endl;
//   }
//   cout << endl;*/
// }


/* 
**********************************************
***       COMPLEX  MATRICES                ***             
**********************************************
*/

template<typename real_scalartype>
class invert_plan<std::complex<real_scalartype> >
{

 public:

  invert_plan(int n);
  invert_plan(int n, int LDA);

  ~invert_plan();

  void execute_plan();
//   void check_plan();

 private:

  void reset_inverted_matrix();

  int N;
  int LDA;
  int LWORK;
  int INFO;

 public:

  std::complex<real_scalartype>*  Matrix;
  std::complex<real_scalartype>*  inverted_matrix;
  int*                            IPIV;

  const static bool check = true;
};

// extern "C" void cgesv_(const int *N, const int *NRHS, std::complex<float> *A, const int *LDA, int *IPIV, std::complex<float> *B, const int *LDB, int *INFO);
// extern "C" void zgesv_(const int *N, const int *NRHS, std::complex<double> *A, const int *LDA, int *IPIV, std::complex<double> *B, const int *LDB, int *INFO);

template<typename real_scalartype>
invert_plan<std::complex<real_scalartype> >::invert_plan(int n):
  N(n),
  LDA(N)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  Matrix          = new std::complex<real_scalartype>[N*N];
  memset(Matrix, 0, sizeof(std::complex<real_scalartype>)*N*N);

  inverted_matrix = new std::complex<real_scalartype>[N*N];
  memset(inverted_matrix, 0, sizeof(std::complex<real_scalartype>)*N*N);
  
  for(int i=0; i<N; i++)
    inverted_matrix[i+N*i] = std::complex<real_scalartype>(1.);

  IPIV = new int[N];
}

template<typename real_scalartype>
invert_plan<std::complex<real_scalartype> >::invert_plan(int n, int lda):
  N(n),
  LDA(lda)
{
  //cout << __PRETTY_FUNCTION__ << endl;

  Matrix      = new std::complex<real_scalartype>[LDA*LDA];
  memset(Matrix, 0, sizeof(std::complex<real_scalartype>)*LDA*LDA);

  inverted_matrix      = new std::complex<real_scalartype>[LDA*LDA];
  memset(inverted_matrix, 0, sizeof(std::complex<real_scalartype>)*LDA*LDA);
  
  for(int i=0; i<LDA; i++)
    inverted_matrix[i+LDA*i] = std::complex<real_scalartype>(1.);

  IPIV = new int[LDA];
}

template<typename real_scalartype>
invert_plan<std::complex<real_scalartype> >::~invert_plan()
{
  //cout << __PRETTY_FUNCTION__ << " is terminated" << endl;

  delete [] Matrix;
  delete [] inverted_matrix;
  delete [] IPIV;
}

template<typename real_scalartype>
void invert_plan<std::complex<real_scalartype> >::reset_inverted_matrix()
{
  int lda = std::max(N,LDA);

  memset(inverted_matrix, 0, sizeof(std::complex<real_scalartype>)*lda*lda);
  
  for(int i=0; i<lda; i++)
    inverted_matrix[i+lda*i] = std::complex<real_scalartype>(1.);
}

template<typename real_scalartype>
void invert_plan<std::complex<real_scalartype> >::execute_plan()
{
  //cout << __PRETTY_FUNCTION__ << endl;
  throw std::logic_error(__PRETTY_FUNCTION__);
}

template<>
void invert_plan<std::complex<float> >::execute_plan()
{
  reset_inverted_matrix();

  LAPACK::cgesv_(&N, &N, Matrix, &LDA, IPIV, inverted_matrix, &LDA, &INFO);

//   if(check && INFO != 0)
//     check_plan();

  assert(INFO == 0);
}

template<>
void invert_plan<std::complex<double> >::execute_plan()
{
  reset_inverted_matrix();

  LAPACK::zgesv_(&N, &N, Matrix, &LDA, IPIV, inverted_matrix, &LDA, &INFO);

//   if(check && INFO != 0)
//     check_plan();

  assert(INFO == 0);
}

// template<typename real_scalartype>
// void invert_plan<std::complex<real_scalartype> >::check_plan()
// {
//   cout << endl << endl <<  __PRETTY_FUNCTION__ << endl << endl;

//   int lda = std::max(N,LDA);

//   for(int i=0; i<N; i++){
//       for(int j=0; j<N; j++){
	
// 	cout << "\t" << Matrix[i+j*lda];
//       }
//       cout << endl;
//   }
//   cout << endl;

//   // makes no sense to check back --> on exit we have Matrix = P*L*U !!!!

//   /*for(int i=0; i<N; i++){
//       for(int j=0; j<N; j++){
// 	std::complex<real_scalartype> c(0,0);
// 	for(int l=0; l<N; l++)
// 	  c += Matrix[i+l*N]*inverted_matrix[l+j*lda];

// 	i==j? c-=1 : c;

// 	cout << "\t" << c;
// 	assert(abs(c)<1.e-6);
//       }
//       cout << endl;
//   }
//   cout << endl;*/
// }

#endif

