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
// Template specialization for CPU.

 template <class Parameters, typename Real>
class N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real> {
  const static int MAX_VERTEX_SINGLETS = 4;

  typedef typename Parameters::concurrency_type concurrency_type;
  typedef typename Parameters::profiler_type profiler_t;

public:
  N_MATRIX_TOOLS(int id,Parameters& parameters_ref);

  Real* get_device_ptr(dca::linalg::Vector<Real, dca::linalg::CPU>& v);

  int* get_permutation();
  void set_permutation(std::vector<int>& p);

  void set_d_vector(dca::linalg::Vector<Real, dca::linalg::CPU>& d_inv);

  void scale_rows(dca::linalg::Matrix<Real, dca::linalg::CPU>& N);

  void copy_rows(dca::linalg::Matrix<Real, dca::linalg::CPU>& N,
                 dca::linalg::Matrix<Real, dca::linalg::CPU>& N_new_spins);

  void compute_G_cols(std::vector<Real>& exp_V, dca::linalg::Matrix<Real, dca::linalg::CPU>& N,
                      dca::linalg::Matrix<Real, dca::linalg::CPU>& G,
                      dca::linalg::Matrix<Real, dca::linalg::CPU>& G_cols);

private:
  int thread_id;
  int stream_id;

 Parameters& parameters;
  concurrency_type& concurrency;

  dca::linalg::Vector<int, dca::linalg::CPU> identity;
  dca::linalg::Vector<int, dca::linalg::CPU> permutation;

  dca::linalg::Vector<Real, dca::linalg::CPU> exp_V;
  dca::linalg::Vector<Real, dca::linalg::CPU> d_vec;

  // dca::linalg::Matrix<Real, dca::linalg::CPU> data;
};

 template <class Parameters, typename Real>
N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real>::N_MATRIX_TOOLS(int id,
                                                                 Parameters& parameters_ref)
    : thread_id(id),
      stream_id(0),

      parameters(parameters_ref),
      concurrency(parameters.get_concurrency()),

      identity(MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),
      permutation(MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),

      exp_V(MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()),
      d_vec(MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size()) {
  for (int l = 0; l < MAX_VERTEX_SINGLETS * parameters.get_max_submatrix_size(); ++l)
    identity[l] = l;
}

 template <class Parameters, typename Real>
Real* N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real>::get_device_ptr(
    dca::linalg::Vector<Real, dca::linalg::CPU>& v) {
  return v.ptr();
}

 template <class Parameters, typename Real>
int* N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real>::get_permutation() {
  return permutation.ptr();
}

 template <class Parameters, typename Real>
void N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real>::set_permutation(std::vector<int>& p) {
  permutation = p;
}

 template <class Parameters, typename Real>
void N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real>::set_d_vector(
    dca::linalg::Vector<Real, dca::linalg::CPU>& d_inv) {
  d_vec = d_inv;
}

 template <class Parameters, typename Real>
void N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real>::scale_rows(
    dca::linalg::Matrix<Real, dca::linalg::CPU>& N) {
  assert(permutation.size() == d_vec.size());

  dca::linalg::matrixop::scaleRows(N, permutation, d_vec, thread_id, stream_id);
}

 template <class Parameters, typename Real>
void N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real>::copy_rows(
    dca::linalg::Matrix<Real, dca::linalg::CPU>& N,
    dca::linalg::Matrix<Real, dca::linalg::CPU>& N_new_spins) {
  assert(N_new_spins.nrCols() == N.nrCols());
  assert(N_new_spins.nrRows() == permutation.size());
  assert(permutation.size() <= identity.size());

  dca::linalg::matrixop::copyRows(N, permutation, N_new_spins, identity, thread_id, stream_id);
}

 template <class Parameters, typename Real>
void N_MATRIX_TOOLS<dca::linalg::CPU,Parameters, Real>::compute_G_cols(
    std::vector<Real>& exp_V, dca::linalg::Matrix<Real, dca::linalg::CPU>& N,
    dca::linalg::Matrix<Real, dca::linalg::CPU>& G,
    dca::linalg::Matrix<Real, dca::linalg::CPU>& G_cols) {
  assert(N.nrRows() == G.nrRows());
  assert(N.nrRows() == G_cols.nrRows());
  assert(permutation.size() == G_cols.nrCols());
  assert(int(exp_V.size()) == permutation.size());

  int N_r = N.nrRows();

  int N_ind = N.nrCols() - G.nrCols();

  for (int l = 0; l < permutation.size(); ++l) {
    if (permutation[l] >= N_ind) {
      for (int i = 0; i < N_r; ++i)
        G_cols(i, l) = G(i, permutation[l] - N_ind);
    }
    else {
      Real alpha = exp_V[l] / (exp_V[l] - 1.);

      for (int i = 0; i < N_r; ++i)
        G_cols(i, l) = alpha * N(i, permutation[l]);

      G_cols(permutation[l], l) -= 1. / (exp_V[l] - 1.);
    }
  }
}
