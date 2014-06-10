//-*-C++-*-

/*
 * invert_plan.h
 *
 *      Author: peter staar
 */

#ifndef SOLVE_PLAN_H_
#define SOLVE_PLAN_H_

template<typename real_scalartype>
class solve_plan
{
public:
  
  solve_plan(int n, int nrhs);
  solve_plan(int n, int lda, int nrhs);
  ~solve_plan();
  
  int execute_plan();
  
public:

  int N;
  int NRHS;
  int LDA;
  int LWORK;
  int INFO;

  real_scalartype*  matrix;
  real_scalartype*  solved_matrix;

private:

  int*              IPIV;
};

template<typename real_scalartype>
solve_plan<real_scalartype>::solve_plan(int n, int nrhs):
  N(n),
  NRHS(nrhs),
  LDA(N)
{
  matrix           = new real_scalartype[N*N];
  memset(matrix, 0, sizeof(real_scalartype)*N*N);

  solved_matrix  = new real_scalartype[N*NRHS];
  memset(solved_matrix, 0, sizeof(real_scalartype)*N*NRHS);
  
  IPIV = new int[N];
}

template<typename real_scalartype>
solve_plan<real_scalartype>::solve_plan(int n, int lda, int nrhs):
  N(n),
  LDA(lda)
{
  matrix      = new real_scalartype[LDA*LDA];
  memset(matrix, 0, sizeof(real_scalartype)*LDA*LDA);

  solved_matrix      = new real_scalartype[LDA*NRHS];
  memset(solved_matrix, 0, sizeof(real_scalartype)*LDA*NRHS);
  
  IPIV = new int[N];
}

template<typename real_scalartype>
solve_plan<real_scalartype>::~solve_plan()
{
  delete [] matrix;
  delete [] solved_matrix;

  delete [] IPIV;
}

template<typename real_scalartype>
int solve_plan<real_scalartype>::execute_plan()
{
  throw std::logic_error(__PRETTY_FUNCTION__);
}

template<>
int solve_plan<float>::execute_plan()
{
  LAPACK::sgesv_(&N, &NRHS, matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);
  return INFO;
}

template<>
int solve_plan<double>::execute_plan()
{
  LAPACK::dgesv_(&N, &NRHS, matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);
  return INFO;
}

template<>
int solve_plan<std::complex<float> >::execute_plan()
{
  LAPACK::cgesv_(&N, &NRHS, matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);
  return INFO;
}

template<>
int solve_plan<std::complex<double> >::execute_plan()
{
  LAPACK::zgesv_(&N, &NRHS, matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);
  return INFO;
}

#endif













// /* 
// **********************************************
// ***       COMPLEX  MATRICES                ***             
// **********************************************
// */




// template<typename real_scalartype>
// class solve_plan<std::complex<real_scalartype> >
// {

//  public:

//   solve_plan(int n);
//   solve_plan(int n, int LDA);

//   ~solve_plan();

//   void execute_plan();
//   void check_plan();

//  private:

//   void reset_solved_matrix();

//   int N;
//   int LDA;
//   int LWORK;
//   int INFO;

//  public:

//   std::complex<real_scalartype>*  Matrix;
//   std::complex<real_scalartype>*  solved_matrix;
//   int*                            IPIV;

//   const static bool check = true;
// };


// extern "C" void cgesv_(const int *N, const int *NRHS, std::complex<float> *A, const int *LDA, int *IPIV, std::complex<float> *B, const int *LDB, int *INFO);
// extern "C" void zgesv_(const int *N, const int *NRHS, std::complex<double> *A, const int *LDA, int *IPIV, std::complex<double> *B, const int *LDB, int *INFO);

// template<typename real_scalartype>
// solve_plan<std::complex<real_scalartype> >::solve_plan(int n):
//   N(n),
//   LDA(N)
// {
//   //cout << __PRETTY_FUNCTION__ << endl;

//   Matrix          = new std::complex<real_scalartype>[N*N];
//   memset(Matrix, 0, sizeof(std::complex<real_scalartype>)*N*N);

//   solved_matrix = new std::complex<real_scalartype>[N*N];
//   memset(solved_matrix, 0, sizeof(std::complex<real_scalartype>)*N*N);
  
//   for(int i=0; i<N; i++)
//     solved_matrix[i+N*i] = std::complex<real_scalartype>(1.);

//   IPIV = new int[N];
// }

// template<typename real_scalartype>
// solve_plan<std::complex<real_scalartype> >::solve_plan(int n, int lda):
//   N(n),
//   LDA(lda)
// {
//   //cout << __PRETTY_FUNCTION__ << endl;

//   Matrix      = new std::complex<real_scalartype>[LDA*LDA];
//   memset(Matrix, 0, sizeof(std::complex<real_scalartype>)*LDA*LDA);

//   solved_matrix      = new std::complex<real_scalartype>[LDA*LDA];
//   memset(solved_matrix, 0, sizeof(std::complex<real_scalartype>)*LDA*LDA);
  
//   for(int i=0; i<LDA; i++)
//     solved_matrix[i+LDA*i] = std::complex<real_scalartype>(1.);

//   IPIV = new int[LDA];
// }

// template<typename real_scalartype>
// solve_plan<std::complex<real_scalartype> >::~solve_plan()
// {
//   //cout << __PRETTY_FUNCTION__ << " is terminated" << endl;

//   delete [] Matrix;
//   delete [] solved_matrix;
//   delete [] IPIV;
// }

// template<typename real_scalartype>
// void solve_plan<std::complex<real_scalartype> >::reset_solved_matrix()
// {
//   int lda = std::max(N,LDA);

//   memset(solved_matrix, 0, sizeof(std::complex<real_scalartype>)*lda*lda);
  
//   for(int i=0; i<lda; i++)
//     solved_matrix[i+lda*i] = std::complex<real_scalartype>(1.);
// }

// template<typename real_scalartype>
// void solve_plan<std::complex<real_scalartype> >::execute_plan()
// {
//   //cout << __PRETTY_FUNCTION__ << endl;
//   throw std::logic_error(__PRETTY_FUNCTION__);
// }

// template<>
// void solve_plan<std::complex<float> >::execute_plan()
// {
//   reset_solved_matrix();

//   cgesv_(&N, &N, Matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);

//   if(check && INFO != 0)
//     check_plan();

//   assert(INFO == 0);
// }

// template<>
// void solve_plan<std::complex<double> >::execute_plan()
// {
//   reset_solved_matrix();

//   zgesv_(&N, &N, Matrix, &LDA, IPIV, solved_matrix, &LDA, &INFO);

//   if(check && INFO != 0)
//     check_plan();

//   assert(INFO == 0);
// }

// template<typename real_scalartype>
// void solve_plan<std::complex<real_scalartype> >::check_plan()
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
// 	  c += Matrix[i+l*N]*solved_matrix[l+j*lda];

// 	i==j? c-=1 : c;

// 	cout << "\t" << c;
// 	assert(abs(c)<1.e-6);
//       }
//       cout << endl;
//   }
//   cout << endl;*/
// }

