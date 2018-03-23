//-*-C++-*-

#ifndef LIN_ALG_PERFORMANCE_INSPECTOR_CPU_H
#define LIN_ALG_PERFORMANCE_INSPECTOR_CPU_H

#include <sys/time.h>
#include <thread>
#include <vector>

namespace LIN_ALG {
template <typename scalartype>
class PERFORMANCE_INSPECTOR<CPU, scalartype> {
public:
  static double execute(int M, int K, int N, scalartype a = 1., scalartype b = 0.);

  static double execute_threadded(int N_th, int M, int K, int N, scalartype a = 1.,
                                  scalartype b = 0.);

private:
  struct matrix_data {
    int thread_id;

    int M;
    int K;
    int N;

    scalartype a;
    scalartype b;

    double GFLOPS;
  };

  static void* execute_gemm(void* ptr);
};

template <typename scalartype>
double PERFORMANCE_INSPECTOR<CPU, scalartype>::execute(int M, int K, int N, scalartype a,
                                                       scalartype b) {
  matrix_data data;

  {
    data.thread_id = 0;

    data.M = M;
    data.K = K;
    data.N = N;

    data.a = a;
    data.b = b;

    data.GFLOPS = 0.;

    execute_gemm(&data);
  }

  return data.GFLOPS;
}

template <typename scalartype>
double PERFORMANCE_INSPECTOR<CPU, scalartype>::execute_threadded(int N_th, int M, int K, int N,
                                                                 scalartype a, scalartype b) {
  std::vector<std::thread> threads;
  std::vector<matrix_data> data(N_th);

  for (int l = 0; l < N_th; ++l) {
    data[l].thread_id = l;

    data[l].M = M;
    data[l].K = K;
    data[l].N = N;

    data[l].a = a;
    data[l].b = b;

    data[l].GFLOPS = 0.;

    threads.push_back(execute_gemm, &data[l]);
  }

  for (int l = 0; l < N_th; ++l)
    threads[l].join();

  double total = 0;
  for (int l = 0; l < N_th; ++l)
    total += data[l].GFLOPS;

  return total;
}

template <typename scalartype>
void* PERFORMANCE_INSPECTOR<CPU, scalartype>::execute_gemm(void* ptr) {
  double N_ITERATIONS = 10.;

  matrix_data* data_ptr = reinterpret_cast<matrix_data*>(ptr);

  matrix<scalartype, CPU> A(std::pair<int, int>(data_ptr->M, data_ptr->K));
  matrix<scalartype, CPU> B(std::pair<int, int>(data_ptr->K, data_ptr->N));
  matrix<scalartype, CPU> C(std::pair<int, int>(data_ptr->M, data_ptr->N));

  for (int j = 0; j < data_ptr->K; ++j)
    for (int i = 0; i < data_ptr->M; ++i)
      A(i, j) = drand48();

  for (int j = 0; j < data_ptr->N; ++j)
    for (int i = 0; i < data_ptr->K; ++i)
      B(i, j) = drand48();

  timeval start;
  gettimeofday(&start, NULL);

  for (double i = 0; i < N_ITERATIONS; ++i)
    dca::linalg::matrixop::gemm(data_ptr->a, A, B, data_ptr->b, C, 0, 0);

  timeval end;
  gettimeofday(&end, NULL);

  double time_start = double(start.tv_sec) + (1.e-6) * double(start.tv_usec);
  double time_end = double(end.tv_sec) + (1.e-6) * double(end.tv_usec);

  double time = time_end - time_start;

  data_ptr->GFLOPS =
      double(2 * (data_ptr->M) * (data_ptr->K) * (data_ptr->N) * N_ITERATIONS) / time * (1.e-9);

  // cout << "\n\t" << data_ptr->thread_id << "\t" <<
  // 2*(data_ptr->M)*(data_ptr->K)*(data_ptr->N)*100 << "\t" << time << "\t" <<  data_ptr->GFLOPS;

  return 0;
}
}

#endif
