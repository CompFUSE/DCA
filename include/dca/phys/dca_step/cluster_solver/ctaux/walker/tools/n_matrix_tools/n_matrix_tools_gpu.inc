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

template <class Parameters>
class N_MATRIX_TOOLS<dca::linalg::GPU, Parameters> {
  const static int MAX_VERTEX_SINGLETS = 4;

  using Real = typename Parameters::Real;
  using Scalar = typename Parameters::Scalar;

  typedef typename Parameters::concurrency_type concurrency_type;
  typedef typename Parameters::profiler_type profiler_t;

  template<class T>
  using HostVector = linalg::util::HostVector<T>;

public:
  N_MATRIX_TOOLS(int id, const Parameters& parameters_ref);

  const int* get_permutation() const;
  void set_permutation(const HostVector<int>& p);

  void set_d_vector(dca::linalg::Vector<Scalar, dca::linalg::CPU>& d_inv);

  void scale_rows(dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N);

  auto get_device_ptr(dca::linalg::Vector<Scalar, dca::linalg::CPU>& v);

  void copy_rows(dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N,
                 dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N_new_spins);

  void compute_G_cols(HostVector<Scalar>& exp_V, dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N,
                      dca::linalg::Matrix<Scalar, dca::linalg::GPU>& G,
                      dca::linalg::Matrix<Scalar, dca::linalg::GPU>& G_cols);

private:
  int thread_id;
  int stream_id;

  const Parameters& parameters;
  const concurrency_type& concurrency;

  dca::linalg::Vector<int, dca::linalg::GPU> identity;
  dca::linalg::Vector<int, dca::linalg::GPU> permutation;

  dca::linalg::Vector<Scalar, dca::linalg::GPU> tmp;

  dca::linalg::Vector<Scalar, dca::linalg::GPU> exp_V;

  dca::linalg::Vector<int, dca::linalg::GPU> d_ind;
  dca::linalg::Vector<Scalar, dca::linalg::GPU> d_vec;
};

template <class Parameters>
N_MATRIX_TOOLS<dca::linalg::GPU, Parameters>::N_MATRIX_TOOLS(int id, const Parameters& parameters_ref)
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

template <class Parameters>
const int* N_MATRIX_TOOLS<dca::linalg::GPU, Parameters>::get_permutation() const {
  return permutation.ptr();
}

template <class Parameters>
void N_MATRIX_TOOLS<dca::linalg::GPU, Parameters>::set_permutation(const HostVector<int>& p) {
  permutation.setAsync(p, thread_id, stream_id);
}

template <class Parameters>
void N_MATRIX_TOOLS<dca::linalg::GPU, Parameters>::set_d_vector(
    dca::linalg::Vector<Scalar, dca::linalg::CPU>& d_inv) {
  d_vec.setAsync(d_inv, thread_id, stream_id);
}

template <class Parameters>
void N_MATRIX_TOOLS<dca::linalg::GPU, Parameters>::scale_rows(
    dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N) {
  assert(permutation.size() == d_vec.size());

  dca::linalg::matrixop::scaleRows(N, permutation, d_vec, thread_id, stream_id);
}

template <class Parameters>
auto N_MATRIX_TOOLS<dca::linalg::GPU, Parameters>::get_device_ptr(
    dca::linalg::Vector<Scalar, dca::linalg::CPU>& v) {
  tmp.setAsync(v, thread_id, stream_id);

  return tmp.ptr();
}

template <class Parameters>
void N_MATRIX_TOOLS<dca::linalg::GPU, Parameters>::copy_rows(
    dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N,
    dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N_new_spins) {
  assert(N_new_spins.nrCols() == N.nrCols());
  assert(N_new_spins.nrRows() == permutation.size());
  assert(permutation.size() <= identity.size());
  dca::linalg::util::syncStream(thread_id, stream_id);
  dca::linalg::matrixop::copyRows(N, permutation, N_new_spins, identity, thread_id, stream_id);
}

template <class Parameters>
void N_MATRIX_TOOLS<dca::linalg::GPU, Parameters>::compute_G_cols(
    HostVector<Scalar>& exp_V_CPU, dca::linalg::Matrix<Scalar, dca::linalg::GPU>& N,
    dca::linalg::Matrix<Scalar, dca::linalg::GPU>& G,
    dca::linalg::Matrix<Scalar, dca::linalg::GPU>& G_cols) {
  exp_V.setAsync(exp_V_CPU, linalg::util::getStream(thread_id, stream_id));
  dca::linalg::util::syncStream(thread_id, stream_id);

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
