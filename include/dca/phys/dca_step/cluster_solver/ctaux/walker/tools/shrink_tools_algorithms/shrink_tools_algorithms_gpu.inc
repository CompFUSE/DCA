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
// Template specialization for GPU.

template <typename Real>
class SHRINK_TOOLS_ALGORITHMS<dca::linalg::GPU, Real> {
public:
  SHRINK_TOOLS_ALGORITHMS(int id)
      : thread_id(id),
        stream_id(0),

        i_s_dn("i_s_dn SHRINK_TOOLS_ALGORITHMS<dca::linalg::GPU>", 512),
        i_t_dn("i_t_dn SHRINK_TOOLS_ALGORITHMS<dca::linalg::GPU>", 512),

        i_s_up("i_s_up SHRINK_TOOLS_ALGORITHMS<dca::linalg::GPU>", 512),
        i_t_up("i_t_up SHRINK_TOOLS_ALGORITHMS<dca::linalg::GPU>", 512),

        SHRINK_TOOLS_ALGORITHMS_CPU_obj(id) {}

  void execute(std::vector<int>& source_index_up, std::vector<int>& target_index_up,
               dca::linalg::Matrix<Real, dca::linalg::GPU>& N_up,
               dca::linalg::Matrix<Real, dca::linalg::GPU>& G0_up, std::vector<int>& source_index_dn,
               std::vector<int>& target_index_dn, dca::linalg::Matrix<Real, dca::linalg::GPU>& N_dn,
               dca::linalg::Matrix<Real, dca::linalg::GPU>& G0_dn) {
    assert(source_index_up.size() == target_index_up.size());
    assert(source_index_dn.size() == target_index_dn.size());

    assert(N_up.size() == G0_up.size());
    assert(N_dn.size() == G0_dn.size());

#ifdef DCA_WITH_QMC_BIT
    N_dn_CPU = N_dn;
    N_up_CPU = N_up;

    G0_dn_CPU = G0_dn;
    G0_up_CPU = G0_up;
#endif  // DCA_WITH_QMC_BIT

    if (true) {
      i_s_up.setAsync(source_index_up, thread_id, stream_id);
      i_t_up.setAsync(target_index_up, thread_id, stream_id);

      i_s_dn.setAsync(source_index_dn, thread_id, stream_id);
      i_t_dn.setAsync(target_index_dn, thread_id, stream_id);

      linalg::util::syncStream(thread_id, stream_id);

      dca::linalg::matrixop::swapRows(N_up, i_s_up, i_t_up, thread_id, stream_id);
      dca::linalg::matrixop::swapCols(N_up, i_s_up, i_t_up, thread_id, stream_id);

      dca::linalg::matrixop::swapRows(G0_up, i_s_up, i_t_up, thread_id, stream_id);
      dca::linalg::matrixop::swapCols(G0_up, i_s_up, i_t_up, thread_id, stream_id);

      dca::linalg::matrixop::swapRows(N_dn, i_s_dn, i_t_dn, thread_id, stream_id);
      dca::linalg::matrixop::swapCols(N_dn, i_s_dn, i_t_dn, thread_id, stream_id);

      dca::linalg::matrixop::swapRows(G0_dn, i_s_dn, i_t_dn, thread_id, stream_id);
      dca::linalg::matrixop::swapCols(G0_dn, i_s_dn, i_t_dn, thread_id, stream_id);
    }
    else {
      for (size_t l = 0; l < source_index_up.size(); ++l) {
        dca::linalg::matrixop::swapRowAndCol(N_up, source_index_up[l], target_index_up[l],
                                             thread_id, stream_id);
        dca::linalg::matrixop::swapRowAndCol(G0_up, source_index_up[l], target_index_up[l],
                                             thread_id, stream_id);
      }

      for (size_t l = 0; l < source_index_dn.size(); ++l) {
        dca::linalg::matrixop::swapRowAndCol(N_dn, source_index_dn[l], target_index_dn[l],
                                             thread_id, stream_id);
        dca::linalg::matrixop::swapRowAndCol(G0_dn, source_index_dn[l], target_index_dn[l],
                                             thread_id, stream_id);
      }
    }

#ifdef DCA_WITH_QMC_BIT
    SHRINK_TOOLS_ALGORITHMS_CPU_obj.execute(source_index_up, target_index_up, N_up_CPU, G0_up_CPU,
                                            source_index_dn, target_index_dn, N_dn_CPU, G0_dn_CPU);

    dca::linalg::matrixop::difference(N_dn_CPU, N_dn);
    dca::linalg::matrixop::difference(N_up_CPU, N_up);

    dca::linalg::matrixop::difference(G0_dn_CPU, G0_dn);
    dca::linalg::matrixop::difference(G0_up_CPU, G0_up);
#endif  // DCA_WITH_QMC_BIT
  }

  int deviceFingerprint() const {
    return (i_s_dn.capacity() + i_t_up.capacity()) * 2 * sizeof(int);
  }

private:
  bool test_swap_vectors(std::vector<int>& source_index, std::vector<int>& target_index) {
    if (source_index.size() != target_index.size())
      throw std::logic_error("source_index.size() != target_index.size()");

    for (size_t i = 0; i < source_index.size(); ++i)
      for (size_t j = i + 1; j < source_index.size(); ++j)
        if (source_index[i] == source_index[j])
          throw std::logic_error("source_index[i] == source_index[j]");

    for (size_t i = 0; i < target_index.size(); ++i)
      for (size_t j = i + 1; j < target_index.size(); ++j)
        if (target_index[i] == target_index[j])
          throw std::logic_error("target_index[i] == target_index[j]");

    for (size_t i = 0; i < source_index.size(); ++i)
      for (size_t j = 0; j < target_index.size(); ++j)
        if (source_index[i] == target_index[j])
          throw std::logic_error("source_index[i] == target_index[j]");

    return true;
  }

private:
  int thread_id;
  int stream_id;

  dca::linalg::Vector<int, dca::linalg::GPU> i_s_dn;
  dca::linalg::Vector<int, dca::linalg::GPU> i_t_dn;

  dca::linalg::Vector<int, dca::linalg::GPU> i_s_up;
  dca::linalg::Vector<int, dca::linalg::GPU> i_t_up;

  dca::linalg::Matrix<Real, dca::linalg::CPU> N_dn_CPU;
  dca::linalg::Matrix<Real, dca::linalg::CPU> G0_dn_CPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> N_up_CPU;
  dca::linalg::Matrix<Real, dca::linalg::CPU> G0_up_CPU;

  SHRINK_TOOLS_ALGORITHMS<dca::linalg::CPU, Real> SHRINK_TOOLS_ALGORITHMS_CPU_obj;
};

template class SHRINK_TOOLS_ALGORITHMS<dca::linalg::GPU, float>;
template class SHRINK_TOOLS_ALGORITHMS<dca::linalg::GPU, double>;
