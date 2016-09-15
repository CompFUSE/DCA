// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization for CPU.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_SHRINK_ROUTINES_CTAUX_SHRINK_ROUTINES_CPU_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_SHRINK_ROUTINES_CTAUX_SHRINK_ROUTINES_CPU_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_shrink_routines/ctaux_shrink_routines_TEM.h"

#include <cassert>
#include <vector>

#include "comp_library/linalg/linalg.hpp"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

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

private:
  int thread_id;
  int stream_id;
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_SHRINK_ROUTINES_CTAUX_SHRINK_ROUTINES_CPU_H
