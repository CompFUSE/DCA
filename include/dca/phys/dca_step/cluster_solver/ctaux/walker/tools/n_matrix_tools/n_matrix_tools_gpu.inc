// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// N-matrix tools class.
// Template specialization for GPU.

template <typename parameters_type>
class N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type> {
  const static int MAX_VERTEX_SINGLETS = 4;

  typedef typename parameters_type::concurrency_type concurrency_type;
  typedef typename parameters_type::profiler_type profiler_t;

public:
  N_MATRIX_TOOLS(int id, parameters_type& parameters_ref);

  int* get_permutation();
  void set_permutation(std::vector<int>& p);

  void set_d_vector(dca::linalg::Vector<double, dca::linalg::CPU>& d_inv);

  void scale_rows(dca::linalg::Matrix<double, dca::linalg::GPU>& N);

  double* get_device_ptr(dca::linalg::Vector<double, dca::linalg::CPU>& v);

  void copy_rows(dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                 dca::linalg::Matrix<double, dca::linalg::GPU>& N_new_spins);

  void compute_G_cols(std::vector<double>& exp_V, dca::linalg::Matrix<double, dca::linalg::GPU>& N,
                      dca::linalg::Matrix<double, dca::linalg::GPU>& G,
                      dca::linalg::Matrix<double, dca::linalg::GPU>& G_cols);

private:
  int thread_id;
  int stream_id;

  parameters_type& parameters;
  concurrency_type& concurrency;

  dca::linalg::Vector<int, dca::linalg::GPU> identity;
  dca::linalg::Vector<int, dca::linalg::GPU> permutation;

  dca::linalg::Vector<double, dca::linalg::GPU> tmp;

  dca::linalg::Vector<double, dca::linalg::GPU> exp_V;

  dca::linalg::Vector<int, dca::linalg::GPU> d_ind;
  dca::linalg::Vector<double, dca::linalg::GPU> d_vec;
};

template <typename parameters_type>
N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::N_MATRIX_TOOLS(int id,
                                                                  parameters_type& parameters_ref)
    : thread_id(id),
      stream_id(0),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      identity("identity    N_MATRIX_TOOLS<dca::linalg::GPU>",
               MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),
      permutation("permutation N_MATRIX_TOOLS<dca::linalg::GPU>",
                  MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),

      tmp("tmp   N_MATRIX_TOOLS<dca::linalg::GPU>",
          MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),
      exp_V("exp_V N_MATRIX_TOOLS<dca::linalg::GPU>",
            MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),

      d_ind("d_ind N_MATRIX_TOOLS<dca::linalg::GPU>",
            MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),
      d_vec("d_vec N_MATRIX_TOOLS<dca::linalg::GPU>",
            MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()) {
  std::vector<int> id_tmp(MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size());

  for (int l = 0; l < MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size(); ++l)
    id_tmp[l] = l;

  identity = id_tmp;
}

template <typename parameters_type>
int* N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::get_permutation() {
  return permutation.ptr();
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::set_permutation(std::vector<int>& p) {
  permutation.setAsync(p, linalg::util::getStream(thread_id, stream_id));
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::set_d_vector(
    dca::linalg::Vector<double, dca::linalg::CPU>& d_inv) {
  d_vec.set(d_inv, thread_id, stream_id);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::scale_rows(
    dca::linalg::Matrix<double, dca::linalg::GPU>& N) {
  assert(permutation.size() == d_vec.size());

  dca::linalg::matrixop::scaleRows(N, permutation, d_vec, thread_id, stream_id);
}

template <typename parameters_type>
double* N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::get_device_ptr(
    dca::linalg::Vector<double, dca::linalg::CPU>& v) {
  tmp.set(v, thread_id, stream_id);

  return tmp.ptr();
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::copy_rows(
    dca::linalg::Matrix<double, dca::linalg::GPU>& N,
    dca::linalg::Matrix<double, dca::linalg::GPU>& N_new_spins) {
  assert(N_new_spins.nrCols() == N.nrCols());
  assert(N_new_spins.nrRows() == permutation.size());
  assert(permutation.size() <= identity.size());

  dca::linalg::matrixop::copyRows(N, permutation, N_new_spins, identity, thread_id, stream_id);
}

template <typename parameters_type>
void N_MATRIX_TOOLS<dca::linalg::GPU, parameters_type>::compute_G_cols(
    std::vector<double>& exp_V_CPU, dca::linalg::Matrix<double, dca::linalg::GPU>& N,
    dca::linalg::Matrix<double, dca::linalg::GPU>& G,
    dca::linalg::Matrix<double, dca::linalg::GPU>& G_cols) {
  exp_V.setAsync(exp_V_CPU, linalg::util::getStream(thread_id, stream_id));

  assert(N.nrRows() == G.nrRows());
  assert(N.nrRows() == G_cols.nrRows());
  assert(exp_V.size() == permutation.size());

  int N_i = permutation.size();
  int N_r = N.nrRows();

  int N_c = N.nrCols() - G.nrCols();

  nkernels::compute_G_cols(N_i, N_r, N_c, permutation.ptr(), exp_V.ptr(), N.ptr(),
                           N.leadingDimension(), G.ptr(), G.leadingDimension(), G_cols.ptr(),
                           G_cols.leadingDimension(), thread_id, stream_id);
}
