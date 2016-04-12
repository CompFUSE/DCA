//-*-C++-*-

#ifndef DCA_QMCI_N_MATRIX_ROUTINES_CPU_H
#define DCA_QMCI_N_MATRIX_ROUTINES_CPU_H

#include "ctaux_N_matrix_routines_TEM.h"

namespace DCA {
namespace QMCI {

template <typename parameters_type>
class N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type> {
  const static int MAX_VERTEX_SINGLETS = 4;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  N_MATRIX_TOOLS(int id, parameters_type& parameters_ref);
  ~N_MATRIX_TOOLS();

  double* get_device_ptr(LIN_ALG::vector<double, LIN_ALG::CPU>& v);

  int* get_permutation();
  void set_permutation(std::vector<int>& p);

  void set_d_vector(std::vector<int>& d_index, LIN_ALG::matrix<double, LIN_ALG::CPU>& N,
                    LIN_ALG::vector<double, LIN_ALG::CPU>& d_inv);

  void set_d_vector(LIN_ALG::vector<double, LIN_ALG::CPU>& d_inv);

  void scale_rows(LIN_ALG::matrix<double, LIN_ALG::CPU>& N);

  void copy_rows(LIN_ALG::matrix<double, LIN_ALG::CPU>& N,
                 LIN_ALG::matrix<double, LIN_ALG::CPU>& N_new_spins);

  void compute_G_cols(std::vector<double>& exp_V, LIN_ALG::matrix<double, LIN_ALG::CPU>& N,
                      LIN_ALG::matrix<double, LIN_ALG::CPU>& G,
                      LIN_ALG::matrix<double, LIN_ALG::CPU>& G_cols);

private:
  int thread_id;
  int stream_id;

  parameters_type& parameters;
  concurrency_type& concurrency;

  LIN_ALG::vector<int, LIN_ALG::CPU> identity;
  LIN_ALG::vector<int, LIN_ALG::CPU> permutation;

  LIN_ALG::vector<double, LIN_ALG::CPU> exp_V;
  LIN_ALG::vector<double, LIN_ALG::CPU> d_vec;

  // LIN_ALG::matrix<double, LIN_ALG::CPU> data;
};

template <typename parameters_type>
N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::N_MATRIX_TOOLS(int id, parameters_type& parameters_ref)
    : thread_id(id),
      stream_id(0),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      identity(MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),
      permutation(MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),

      exp_V(MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()),
      d_vec(MAX_VERTEX_SINGLETS * parameters.get_K_PHANI()) {
  for (int l = 0; l < MAX_VERTEX_SINGLETS * parameters.get_K_PHANI(); ++l)
    identity[l] = l;
}

template <typename parameters_type>
N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::~N_MATRIX_TOOLS() {}

template <typename parameters_type>
double* N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::get_device_ptr(
    LIN_ALG::vector<double, LIN_ALG::CPU>& v) {
  return v.get_ptr();
}

template <typename parameters_type>
int* N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::get_permutation() {
  return permutation.get_ptr();
}

template <typename parameters_type>
void N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::set_permutation(std::vector<int>& p) {
  permutation.set(p);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::set_d_vector(
    LIN_ALG::vector<double, LIN_ALG::CPU>& d_inv) {
  d_vec.set(d_inv);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::scale_rows(LIN_ALG::matrix<double, LIN_ALG::CPU>& N) {
  assert(permutation.size() == d_vec.size());

  int N_i = permutation.size();
  int N_c = N.get_number_of_cols();

  int N_LD = N.get_leading_dimension();

  LIN_ALG::SCALE<LIN_ALG::CPU>::many_rows(N_c, N_i, permutation.get_ptr(), d_vec.get_ptr(),
                                          N.get_ptr(), N_LD, thread_id, stream_id);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::copy_rows(
    LIN_ALG::matrix<double, LIN_ALG::CPU>& N, LIN_ALG::matrix<double, LIN_ALG::CPU>& N_new_spins) {
  assert(N_new_spins.get_number_of_cols() == N.get_number_of_cols());
  assert(N_new_spins.get_number_of_rows() == permutation.size());

  int N_i = permutation.size();
  int N_c = N.get_number_of_cols();

  assert(N_i <= identity.size());

  LIN_ALG::COPY<LIN_ALG::CPU>::many_rows(
      N_c, N_i, permutation.get_ptr(), N.get_ptr(), N.get_leading_dimension(), identity.get_ptr(),
      N_new_spins.get_ptr(), N_new_spins.get_leading_dimension(), thread_id, stream_id);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<LIN_ALG::CPU, parameters_type>::compute_G_cols(
    std::vector<double>& exp_V, LIN_ALG::matrix<double, LIN_ALG::CPU>& N,
    LIN_ALG::matrix<double, LIN_ALG::CPU>& G, LIN_ALG::matrix<double, LIN_ALG::CPU>& G_cols) {
  assert(N.get_number_of_rows() == G.get_number_of_rows());
  assert(N.get_number_of_rows() == G_cols.get_number_of_rows());
  assert(permutation.size() == G_cols.get_number_of_cols());
  assert(int(exp_V.size()) == permutation.size());

  int N_r = N.get_number_of_rows();

  int N_ind = N.get_number_of_cols() - G.get_number_of_cols();

  for (int l = 0; l < permutation.size(); ++l) {
    if (permutation[l] >= N_ind) {
      for (int i = 0; i < N_r; ++i)
        G_cols(i, l) = G(i, permutation[l] - N_ind);
    }
    else {
      double alpha = exp_V[l] / (exp_V[l] - 1.);

      for (int i = 0; i < N_r; ++i)
        G_cols(i, l) = alpha * N(i, permutation[l]);

      G_cols(permutation[l], l) -= 1. / (exp_V[l] - 1.);
    }
  }
}
}
}

#endif
