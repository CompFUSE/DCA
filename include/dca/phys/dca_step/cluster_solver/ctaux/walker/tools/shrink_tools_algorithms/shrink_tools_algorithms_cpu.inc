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

template <>
class SHRINK_TOOLS_ALGORITHMS<dca::linalg::CPU> {
public:
  SHRINK_TOOLS_ALGORITHMS(int id) : thread_id(id), stream_id(0) {}

  void execute(std::vector<int>& source_index, std::vector<int>& target_index,
               dca::linalg::Matrix<double, dca::linalg::CPU>& N,
               dca::linalg::Matrix<double, dca::linalg::CPU>& G0) {
    assert(source_index.size() == target_index.size());

    assert(N.nrRows() == G0.nrRows());
    assert(N.nrCols() == G0.nrCols());

    for (size_t l = 0; l < source_index.size(); ++l) {
      dca::linalg::matrixop::swapRowAndCol(N, source_index[l], target_index[l], thread_id, stream_id);
      dca::linalg::matrixop::swapRowAndCol(G0, source_index[l], target_index[l], thread_id,
                                           stream_id);
    }
  }

  void execute(std::vector<int>& source_index_up, std::vector<int>& target_index_up,
               dca::linalg::Matrix<double, dca::linalg::CPU>& N_up,
               dca::linalg::Matrix<double, dca::linalg::CPU>& G0_up,
               std::vector<int>& source_index_dn, std::vector<int>& target_index_dn,
               dca::linalg::Matrix<double, dca::linalg::CPU>& N_dn,
               dca::linalg::Matrix<double, dca::linalg::CPU>& G0_dn) {
    assert(source_index_up.size() == target_index_up.size());
    assert(source_index_dn.size() == target_index_dn.size());

    assert(N_up.size() == G0_up.size());
    assert(N_dn.size() == G0_dn.size());

    for (size_t l = 0; l < source_index_up.size(); ++l)
      dca::linalg::matrixop::swapRowAndCol(N_up, source_index_up[l], target_index_up[l], thread_id,
                                           stream_id);

    for (size_t l = 0; l < source_index_up.size(); ++l)
      dca::linalg::matrixop::swapRowAndCol(G0_up, source_index_up[l], target_index_up[l], thread_id,
                                           stream_id);

    for (size_t l = 0; l < source_index_dn.size(); ++l)
      dca::linalg::matrixop::swapRowAndCol(N_dn, source_index_dn[l], target_index_dn[l], thread_id,
                                           stream_id);

    for (size_t l = 0; l < source_index_dn.size(); ++l)
      dca::linalg::matrixop::swapRowAndCol(G0_dn, source_index_dn[l], target_index_dn[l], thread_id,
                                           stream_id);
  }

  std::size_t deviceFingerprint() const {
    return 0;
  }

private:
  int thread_id;
  int stream_id;
};
