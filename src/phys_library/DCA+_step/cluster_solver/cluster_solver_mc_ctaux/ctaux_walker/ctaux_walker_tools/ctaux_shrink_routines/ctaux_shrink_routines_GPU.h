// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Template specialization for GPU.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_SHRINK_ROUTINES_CTAUX_SHRINK_ROUTINES_GPU_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_SHRINK_ROUTINES_CTAUX_SHRINK_ROUTINES_GPU_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_shrink_routines/ctaux_shrink_routines_TEM.h"

#include <cassert>
#include <stdexcept>
#include <vector>

#include "comp_library/linalg/linalg.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_shrink_routines/ctaux_shrink_routines_CPU.h"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <>
class SHRINK_TOOLS_ALGORITHMS<LIN_ALG::GPU> {
public:
  SHRINK_TOOLS_ALGORITHMS(int id)
      : thread_id(id),
        stream_id(0),

        i_s_dn("i_s_dn SHRINK_TOOLS_ALGORITHMS<LIN_ALG::GPU>", 512),
        i_t_dn("i_t_dn SHRINK_TOOLS_ALGORITHMS<LIN_ALG::GPU>", 512),

        i_s_up("i_s_up SHRINK_TOOLS_ALGORITHMS<LIN_ALG::GPU>", 512),
        i_t_up("i_t_up SHRINK_TOOLS_ALGORITHMS<LIN_ALG::GPU>", 512),

        SHRINK_TOOLS_ALGORITHMS_CPU_obj(id) {
    i_s_dn.set_thread_and_stream_id(thread_id, stream_id);
    i_t_dn.set_thread_and_stream_id(thread_id, stream_id);

    i_s_up.set_thread_and_stream_id(thread_id, stream_id);
    i_t_up.set_thread_and_stream_id(thread_id, stream_id);
  }

  void execute(std::vector<int>& source_index_up, std::vector<int>& target_index_up,
               LIN_ALG::matrix<double, LIN_ALG::GPU>& N_up,
               LIN_ALG::matrix<double, LIN_ALG::GPU>& G0_up, std::vector<int>& source_index_dn,
               std::vector<int>& target_index_dn, LIN_ALG::matrix<double, LIN_ALG::GPU>& N_dn,
               LIN_ALG::matrix<double, LIN_ALG::GPU>& G0_dn) {
    assert(source_index_up.size() == target_index_up.size());
    assert(source_index_dn.size() == target_index_dn.size());

    assert(N_up.get_current_size() == G0_up.get_current_size());
    assert(N_dn.get_current_size() == G0_dn.get_current_size());

#ifdef DCA_WITH_QMC_BIT
    N_dn_CPU.copy_from(N_dn);
    N_up_CPU.copy_from(N_up);

    G0_dn_CPU.copy_from(G0_dn);
    G0_up_CPU.copy_from(G0_up);
#endif  // DCA_WITH_QMC_BIT

    if (true) {
      i_s_up.set(source_index_up, LIN_ALG::ASYNCHRONOUS);
      i_t_up.set(target_index_up, LIN_ALG::ASYNCHRONOUS);

      i_s_dn.set(source_index_dn, LIN_ALG::ASYNCHRONOUS);
      i_t_dn.set(target_index_dn, LIN_ALG::ASYNCHRONOUS);

      LIN_ALG::SWAP<LIN_ALG::GPU>::many_rows(N_up, i_s_up, i_t_up, thread_id, stream_id);
      LIN_ALG::SWAP<LIN_ALG::GPU>::many_cols(N_up, i_s_up, i_t_up, thread_id, stream_id);

      LIN_ALG::SWAP<LIN_ALG::GPU>::many_rows(G0_up, i_s_up, i_t_up, thread_id, stream_id);
      LIN_ALG::SWAP<LIN_ALG::GPU>::many_cols(G0_up, i_s_up, i_t_up, thread_id, stream_id);

      LIN_ALG::SWAP<LIN_ALG::GPU>::many_rows(N_dn, i_s_dn, i_t_dn, thread_id, stream_id);
      LIN_ALG::SWAP<LIN_ALG::GPU>::many_cols(N_dn, i_s_dn, i_t_dn, thread_id, stream_id);

      LIN_ALG::SWAP<LIN_ALG::GPU>::many_rows(G0_dn, i_s_dn, i_t_dn, thread_id, stream_id);
      LIN_ALG::SWAP<LIN_ALG::GPU>::many_cols(G0_dn, i_s_dn, i_t_dn, thread_id, stream_id);
    }
    else {
      for (size_t l = 0; l < source_index_up.size(); ++l) {
        LIN_ALG::SWAP<LIN_ALG::GPU>::row_and_column(N_up, source_index_up[l], target_index_up[l],
                                                    thread_id, stream_id);
        LIN_ALG::SWAP<LIN_ALG::GPU>::row_and_column(G0_up, source_index_up[l], target_index_up[l],
                                                    thread_id, stream_id);
      }

      for (size_t l = 0; l < source_index_dn.size(); ++l) {
        LIN_ALG::SWAP<LIN_ALG::GPU>::row_and_column(N_dn, source_index_dn[l], target_index_dn[l],
                                                    thread_id, stream_id);
        LIN_ALG::SWAP<LIN_ALG::GPU>::row_and_column(G0_dn, source_index_dn[l], target_index_dn[l],
                                                    thread_id, stream_id);
      }
    }

#ifdef DCA_WITH_QMC_BIT
    SHRINK_TOOLS_ALGORITHMS_CPU_obj.execute(source_index_up, target_index_up, N_up_CPU, G0_up_CPU,
                                            source_index_dn, target_index_dn, N_dn_CPU, G0_dn_CPU);

    N_dn_CPU.difference(N_dn);
    N_up_CPU.difference(N_up);

    G0_dn_CPU.difference(G0_dn);
    G0_up_CPU.difference(G0_up);
  }
#endif  // DCA_WITH_QMC_BIT

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

LIN_ALG::vector<int, LIN_ALG::GPU> i_s_dn;
LIN_ALG::vector<int, LIN_ALG::GPU> i_t_dn;

LIN_ALG::vector<int, LIN_ALG::GPU> i_s_up;
LIN_ALG::vector<int, LIN_ALG::GPU> i_t_up;

LIN_ALG::matrix<double, LIN_ALG::CPU> N_dn_CPU;
LIN_ALG::matrix<double, LIN_ALG::CPU> G0_dn_CPU;

LIN_ALG::matrix<double, LIN_ALG::CPU> N_up_CPU;
LIN_ALG::matrix<double, LIN_ALG::CPU> G0_up_CPU;

SHRINK_TOOLS_ALGORITHMS<LIN_ALG::CPU> SHRINK_TOOLS_ALGORITHMS_CPU_obj;
};

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_TOOLS_CTAUX_SHRINK_ROUTINES_CTAUX_SHRINK_ROUTINES_GPU_H
