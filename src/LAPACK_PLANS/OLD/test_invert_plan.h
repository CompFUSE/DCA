#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "invert_plan.h"

using namespace std;

void test_invert_plan()
{
  srand ( time(NULL) );

  int N = 2000;

  float* mf = new float[N*N];
  double* md = new double[N*N];

  std::complex<float>* mcf = new std::complex<float>[N*N];
  //cout << "OK" << endl;
  //std::complex<double>* mcd = new std:complex<double>[N*N];

  invert_plan<float> inv_pln_f(N);
  invert_plan<double> inv_pln_d(N);

  invert_plan<std::complex<float> > inv_pln_cf(N);
  invert_plan<std::complex<double> > inv_pln_cd(N);

  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      mf[i+N*j] = rand() % 100;
      md[i+N*j] = rand() % 100;
      real(mcf[i+N*j]) = rand() % 100; imag(mcf[i+N*j]) = rand() % 100;
      //real(mcd[i+N*j]) = rand() % 100; imag(mcd[i+N*j]) = rand() % 100;
    }
  }

  {
    profiler p2("float");
    memcpy(inv_pln_f.Matrix, mf, sizeof(float)*N*N);
    inv_pln_f.execute_plan();
  }

  {
    profiler p2("double");
    memcpy(inv_pln_d.Matrix, md, sizeof(double)*N*N);
    inv_pln_d.execute_plan();
  }

  {
    profiler p2("complex float");
    memcpy(inv_pln_cf.Matrix, mcf, sizeof(std::complex<float>)*N*N);
    inv_pln_cf.execute_plan();
  }
  /*
  {
    profiler p2("complex double");
    memcpy(inv_pln_cd.Matrix, mcd, sizeof(std::complex<double>)*N*N);
    inv_pln_cd.execute_plan();
  }
  */
  /*
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      std::complex<float> c(0,0);
      for(int l=0; l<N; l++)
	c += mcf[i+N*l]*inv_pln_cf.identity_Matrix[l+j*N];

      if(i==j)
	c -= 1;

      assert( fabs(real(c)*real(c)+imag(c)*imag(c)) < 10e-6);
    }
    }*/
}
