// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Shrink tools algorithms class.
// Template specialization for CPU.

template <typename Real>
class SHRINK_TOOLS_ALGORITHMS<dca::linalg::CPU, Real> {
  using HostVector = linalg::util::HostVector<int>;

public:
  SHRINK_TOOLS_ALGORITHMS(int id) : thread_id(id), stream_id(0) {}

  void execute(HostVector& source_index, HostVector& target_index,
               dca::linalg::Matrix<Real, dca::linalg::CPU>& N,
               dca::linalg::Matrix<Real, dca::linalg::CPU>& G0) {
    assert(source_index.size() == target_index.size());

    assert(N.nrRows() == G0.nrRows());
    assert(N.nrCols() == G0.nrCols());

    for (size_t l = 0; l < source_index.size(); ++l) {
      dca::linalg::matrixop::copyRow(N, source_index[l], N, target_index[l], thread_id, stream_id);
      dca::linalg::matrixop::copyCol(N, source_index[l], N, target_index[l], thread_id, stream_id);
      dca::linalg::matrixop::copyRow(G0, source_index[l], G0, target_index[l], thread_id, stream_id);
      dca::linalg::matrixop::copyCol(G0, source_index[l], G0, target_index[l], thread_id, stream_id);
    }
  }

  void execute(HostVector& source_index_up, HostVector& target_index_up,
               dca::linalg::Matrix<Real, dca::linalg::CPU>& N_up,
               dca::linalg::Matrix<Real, dca::linalg::CPU>& G0_up, HostVector& source_index_dn,
               HostVector& target_index_dn, dca::linalg::Matrix<Real, dca::linalg::CPU>& N_dn,
               dca::linalg::Matrix<Real, dca::linalg::CPU>& G0_dn) {
    execute(source_index_up, target_index_up, N_up, G0_up);
    execute(source_index_dn, target_index_dn, N_dn, G0_dn);
  }

  std::size_t deviceFingerprint() const {
    return 0;
  }

private:
  int thread_id;
  int stream_id;
};

template class SHRINK_TOOLS_ALGORITHMS<dca::linalg::CPU, float>;
template class SHRINK_TOOLS_ALGORITHMS<dca::linalg::CPU, double>;
