// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// G-matrix tools class.
// Template specialization for GPU.

template <class Parameters>
class G_MATRIX_TOOLS<dca::linalg::GPU, Parameters> {
  const static int MAX_VERTEX_SINGLETS = 4;

  using Real = typename Parameters::Real;
  using Scalar = typename Parameters::Scalar;
  using Concurrency = typename Parameters::concurrency_type;
  using Profiler = typename Parameters::profiler_type;

public:
  G_MATRIX_TOOLS(int id, const Parameters& parameters_ref)
      : thread_id(id),

        parameters(parameters_ref),
        concurrency(parameters.get_concurrency()),

        i_index_gpu(0, MAX_VERTEX_SINGLETS * parameters_ref.get_max_submatrix_size()),
        j_index_gpu(0, MAX_VERTEX_SINGLETS * parameters_ref.get_max_submatrix_size()),
        is_Bennett_gpu(0, MAX_VERTEX_SINGLETS * parameters_ref.get_max_submatrix_size()),
        exp_Vj_gpu(0, MAX_VERTEX_SINGLETS * parameters_ref.get_max_submatrix_size()) {}

  void read_G_matrix_elements(dca::linalg::Vector<int, dca::linalg::CPU>& i_index,
                              dca::linalg::Vector<int, dca::linalg::CPU>& j_index,
                              dca::linalg::Vector<bool, dca::linalg::CPU>& is_Bennett,
                              dca::linalg::Vector<Real, dca::linalg::CPU>& exp_Vj,
                              dca::linalg::Matrix<Real, dca::linalg::GPU>& N,
                              dca::linalg::Matrix<Real, dca::linalg::GPU>& G_precomputed,
                              Real* result_ptr, int incr) {
    assert(i_index.size() == j_index.size());
    assert(i_index.size() == is_Bennett.size());
    assert(i_index.size() == exp_Vj.size());

    i_index_gpu = i_index;
    j_index_gpu = j_index;
    is_Bennett_gpu = is_Bennett;
    exp_Vj_gpu = exp_Vj;

    int vertex_index = N.nrCols() - G_precomputed.nrCols();

    gkernels::read_G_matrix_elements(i_index_gpu.size(), vertex_index, i_index_gpu.ptr(),
                                     j_index_gpu.ptr(), is_Bennett_gpu.ptr(), exp_Vj_gpu.ptr(),
                                     N.ptr(), N.leadingDimension(), G_precomputed.ptr(),
                                     G_precomputed.leadingDimension(), result_ptr, incr);
  }

  void compute_row_on_Gamma_matrix(int row_index, dca::linalg::Vector<int, dca::linalg::GPU>& indices,
                                   dca::linalg::Vector<Real, dca::linalg::GPU>& exp_V,
                                   dca::linalg::Matrix<Real, dca::linalg::GPU>& N,
                                   dca::linalg::Matrix<Real, dca::linalg::GPU>& G_precomputed,
                                   Real* row_ptr, int incr) {
    int vertex_index = N.nrCols() - G_precomputed.nrCols();

    gkernels::compute_row_on_Gamma_matrix(row_index, indices.size(), vertex_index,

                                          indices.ptr(), exp_V.ptr(),

                                          N.ptr(), N.leadingDimension(), G_precomputed.ptr(),
                                          G_precomputed.leadingDimension(), row_ptr, incr);
  }

  void compute_col_on_Gamma_matrix(int col_index, dca::linalg::Vector<int, dca::linalg::GPU>& indices,
                                   dca::linalg::Vector<Real, dca::linalg::GPU>& exp_V,
                                   dca::linalg::Matrix<Real, dca::linalg::GPU>& N,
                                   dca::linalg::Matrix<Real, dca::linalg::GPU>& G_precomputed,
                                   Real* col_ptr, int incr) {
    int vertex_index = N.nrCols() - G_precomputed.nrCols();

    gkernels::compute_col_on_Gamma_matrix(col_index, indices.size(), vertex_index,

                                          indices.ptr(), exp_V.ptr(),

                                          N.ptr(), N.leadingDimension(), G_precomputed.ptr(),
                                          G_precomputed.leadingDimension(), col_ptr, incr);
  }

  int deviceFingerprint() const {
    return (i_index_gpu.capacity() + j_index_gpu.capacity()) * sizeof(int) +
           is_Bennett_gpu.capacity() * sizeof(bool) + exp_Vj_gpu.capacity() * sizeof(Real);
  }

private:
  int thread_id;

  const Parameters& parameters;
  const Concurrency& concurrency;

  dca::linalg::Vector<int, dca::linalg::GPU> i_index_gpu, j_index_gpu;
  dca::linalg::Vector<bool, dca::linalg::GPU> is_Bennett_gpu;
  dca::linalg::Vector<Real, dca::linalg::GPU> exp_Vj_gpu;
};
