// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class computes the nonlocal single-particle Greens-function \f$G^{I}(k_1,k_2)\f$.

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_CTAUX_ACCUMULATOR_NONLOCAL_G_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_CTAUX_ACCUMULATOR_NONLOCAL_G_H

#include <complex>
#include <vector>

#include "comp_library/function_library/include_function_library.h"
#include "comp_library/linalg/src/matrix.h"
#include "math_library/functional_transforms/function_transforms/function_transforms.hpp"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_accumulator/tp_accumulator/ctaux_tp_nft.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_vertex_singleton.h"
#include "phys_library/domains/cluster/cluster_domain.h"
#include "phys_library/domains/Quantum_domain/electron_band_domain.h"
#include "phys_library/domains/Quantum_domain/electron_spin_domain.h"
#include "phys_library/domains/time_and_frequency/frequency_domain_compact.h"

namespace DCA {
namespace QMCI {
namespace CT_AUX_ACCUMULATION {
// DCA::QMCI::CT_AUX_ACCUMULATION::

template <class parameters_type, class MOMS_type>
class accumulator_nonlocal_G {
public:
  typedef typename parameters_type::profiler_type profiler_t;
  typedef typename parameters_type::concurrency_type concurrency_type;

  typedef typename parameters_type::MC_measurement_scalar_type scalar_type;

  typedef typename parameters_type::G4_w1_dmn_t w1_dmn_t;
  typedef typename parameters_type::G4_w2_dmn_t w2_dmn_t;

  using w_VERTEX_EXTENDED = dmn_0<DCA::vertex_frequency_domain<DCA::EXTENDED>>;
  using w_VERTEX_EXTENDED_POS = dmn_0<DCA::vertex_frequency_domain<DCA::EXTENDED_POSITIVE>>;

  using b = dmn_0<electron_band_domain>;
  using s = dmn_0<electron_spin_domain>;

  using r_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     REAL_SPACE, BRILLOUIN_ZONE>>;
  using k_DCA = dmn_0<cluster_domain<double, parameters_type::lattice_type::DIMENSION, CLUSTER,
                                     MOMENTUM_SPACE, BRILLOUIN_ZONE>>;
  using r_dmn_t = r_DCA;
  using k_dmn_t = k_DCA;

  typedef dmn_6<b, b, r_dmn_t, r_dmn_t, w1_dmn_t, w2_dmn_t> b_b_r_r_w_w_dmn_t;
  typedef dmn_6<b, b, k_dmn_t, r_dmn_t, w1_dmn_t, w2_dmn_t> b_b_k_r_w_w_dmn_t;
  typedef dmn_6<b, b, k_dmn_t, k_dmn_t, w1_dmn_t, w2_dmn_t> b_b_k_k_w_w_dmn_t;

  typedef vertex_singleton vertex_singleton_type;

  typedef LIN_ALG::matrix<double, LIN_ALG::CPU> vertex_vertex_matrix_type;

public:
  accumulator_nonlocal_G(parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id);

  double get_GFLOP();

  void initialize();

  void finalize();

  /*!
   * \brief Compute the nonlocal single particle Greens-function \f$G^{I}_{\sigma}(k_1, k_2)\f$,
   * from the M-matrices.
   */
  void execute(std::vector<vertex_singleton_type>& HS_configuration_e_UP,
               vertex_vertex_matrix_type& M_e_UP,
               std::vector<vertex_singleton_type>& HS_configuration_e_DN,
               vertex_vertex_matrix_type& M_e_DN);

  /*!
   * \brief Get a reference to \f$G^{I}_{\uparrow}(k_1, k_2)\f$.
   */
  FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& get_G_k_k_w_w_e_UP() {
    return G2_k_k_w_w_e_UP;
  }

  /*!
   * \brief Get a reference to \f$G^{I}_{\downarrow}(k_1, k_2)\f$.
   */
  FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& get_G_k_k_w_w_e_DN() {
    return G2_k_k_w_w_e_DN;
  }

  /*!
   * \brief print help-functions.
   */
  template <class stream_type>
  void to_JSON(stream_type& ss);

private:
  void FT_M_v_v_2_M_k_k_w_w(std::vector<vertex_singleton_type>& HS_configuration_e_spin,
                            FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& M_k_k_w_w,
                            vertex_vertex_matrix_type& M_e_spin);

  void FT_M_v_v_2_M_k_k_w_w_test();

  void compute_G2_k_k_w_w(FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& M_k_k_w_w,
                          FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w,
                          e_spin_states_type e_spin);

  void compute_G2_k_k_w_w_atomic(std::complex<double>* G0_k_w_l, std::complex<scalar_type>* M_k_w,
                                 std::complex<double>* G0_k_w_r,
                                 std::complex<scalar_type>* G2_k_k_w_w);

  void compute_G2_k_k_w_w_atomic_1(std::complex<double>* G0_k_w_l, std::complex<scalar_type>* M_k_w,
                                   std::complex<double>* G0_k_w_r,
                                   std::complex<scalar_type>* G2_k_k_w_w);

  void compute_G2_k_k_w_w_atomic_default(std::complex<double>* G0_k_w_l,
                                         std::complex<scalar_type>* M_k_w,
                                         std::complex<double>* G0_k_w_r,
                                         std::complex<scalar_type>* G2_k_k_w_w);

private:
  parameters_type& parameters;
  MOMS_type& MOMS;
  concurrency_type& concurrency;

  int thread_id;

  double GFLOP;

  scalar_type beta;

  int b_dmn_size;
  int H_dmn_size;

  std::vector<int> corresponding_w1_index;
  std::vector<int> corresponding_w2_index;

  // QMC::cached_nft<2, scalar_type, r_dmn_t, w1_dmn_t, w2_dmn_t> nft_obj;
  // QMC::cached_nft<2, scalar_type, r_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED_POS> nft_obj;
  cached_nft<2, scalar_type, r_dmn_t, w_VERTEX_EXTENDED, w_VERTEX_EXTENDED_POS> nft_obj;

  FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t> G2_k_k_w_w_e_UP;
  FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t> G2_k_k_w_w_e_DN;

  FUNC_LIB::function<std::complex<scalar_type>, b_b_r_r_w_w_dmn_t> M_r_r_w_w;
  FUNC_LIB::function<std::complex<scalar_type>, b_b_k_r_w_w_dmn_t> M_k_r_w_w;
  FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t> M_k_k_w_w;
};

template <class parameters_type, class MOMS_type>
accumulator_nonlocal_G<parameters_type, MOMS_type>::accumulator_nonlocal_G(
    parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      thread_id(id),

      GFLOP(0),

      beta(parameters.get_beta()),

      b_dmn_size(b::dmn_size()),
      H_dmn_size(b::dmn_size() * s::dmn_size()),

      corresponding_w1_index(w1_dmn_t::parameter_type::get_corresponding_frequency_domain_index()),
      corresponding_w2_index(w2_dmn_t::parameter_type::get_corresponding_frequency_domain_index()),

      nft_obj(w_VERTEX_EXTENDED::parameter_type::get_elements()),

      G2_k_k_w_w_e_UP("G2_k_k_w_w_e_UP"),
      G2_k_k_w_w_e_DN("G2_k_k_w_w_e_DN"),

      M_r_r_w_w("M_r_r_w_w_tmp"),
      M_k_r_w_w("M_k_r_w_w_tmp"),
      M_k_k_w_w("M_k_k_w_w_tmp") {}

template <class parameters_type, class MOMS_type>
double accumulator_nonlocal_G<parameters_type, MOMS_type>::get_GFLOP() {
  double result = GFLOP;
  GFLOP = 0;

  return result;
}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_G<parameters_type, MOMS_type>::initialize() {
  beta = parameters.get_beta();

  b_dmn_size = b::dmn_size();
  H_dmn_size = b::dmn_size() * s::dmn_size();
}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_G<parameters_type, MOMS_type>::finalize() {}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_G<parameters_type, MOMS_type>::execute(
    std::vector<vertex_singleton_type>& HS_configuration_e_UP, vertex_vertex_matrix_type& M_e_UP,
    std::vector<vertex_singleton_type>& HS_configuration_e_DN, vertex_vertex_matrix_type& M_e_DN) {
  profiler_t profiler("compute nonlocal-G from M-matrix", "CT-AUX accumulator", __LINE__, thread_id);

  { // e_UP
   {// profiler_t profiler("nonlocal-G-FT-M          ", __FILE__, __LINE__);
    FT_M_v_v_2_M_k_k_w_w(HS_configuration_e_UP, M_k_k_w_w, M_e_UP);
}

{
  // profiler_t profiler("(B)\tb. c. --> nonlocal-G-compute       ", __FILE__, __LINE__);
  compute_G2_k_k_w_w(M_k_k_w_w, G2_k_k_w_w_e_UP, e_UP);
}
}

{  // e_DN
  {
    // profiler_t profiler("(B)\tb. c. --> nonlocal-G-FT-M          ", __FILE__, __LINE__);
    FT_M_v_v_2_M_k_k_w_w(HS_configuration_e_DN, M_k_k_w_w, M_e_DN);
  }

  {
    // profiler_t profiler("(B)\tb. c. --> nonlocal-G-compute       ", __FILE__, __LINE__);
    compute_G2_k_k_w_w(M_k_k_w_w, G2_k_k_w_w_e_DN, e_DN);
  }
}
}

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_G<parameters_type, MOMS_type>::FT_M_v_v_2_M_k_k_w_w(
    std::vector<vertex_singleton_type>& configuration_e_spin,
    FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& M_k_k_w_w_e_spin,
    vertex_vertex_matrix_type& M_e_spin) {
  {
    profiler_t profiler("nonlocal-G-FT t --> w", "CT-AUX accumulator", __LINE__, thread_id);
    GFLOP += nft_obj.execute(configuration_e_spin, M_e_spin, M_r_r_w_w);
  }

  // FT_M_v_v_2_M_k_k_w_w_test();

  {
    profiler_t profiler("nonlocal-G-FT r --> k", "CT-AUX accumulator", __LINE__, thread_id);

    // FT<r_dmn_t, k_dmn_t>::execute(M_r_r_w_w, M_k_k_w_w_e_spin);
    math_algorithms::functional_transforms::TRANSFORM<r_dmn_t, k_dmn_t>::execute(M_r_r_w_w,
                                                                                 M_k_r_w_w);
    math_algorithms::functional_transforms::TRANSFORM<r_dmn_t, k_dmn_t>::execute(M_k_r_w_w,
                                                                                 M_k_k_w_w_e_spin);

    // scalar_type  Nc        = scalar_type(base_cluster_type::get_cluster_size());
    scalar_type Nc = scalar_type(k_DCA::dmn_size());
    scalar_type one_div_Nc = 1. / Nc;  // scalar_type(base_cluster_type::get_cluster_size());
    M_k_k_w_w_e_spin *= one_div_Nc;

    GFLOP += M_r_r_w_w.size() * double(Nc) * (1.e-9);
  }
}

/*
template<class parameters_type, class MOMS_type>
void accumulator_nonlocal_G<parameters_type, MOMS_type>::FT_M_v_v_2_M_k_k_w_w_test()
{
  cout << __FUNCTION__ << endl;

  {
    typedef dmn_6<b,b,r_DCA,r_DCA,w1_dmn_t,w2_dmn_t> b_b_r_DCA_r_DCA_w_w_dmn_t;
    typedef dmn_6<b,b,r_PCM,r_PCM,w1_dmn_t,w2_dmn_t> b_b_r_PCM_r_PCM_w_w_dmn_t;

    typedef dmn_6<b,b,k_PCM,r_PCM,w1_dmn_t,w2_dmn_t> b_b_k_PCM_r_PCM_w_w_dmn_t;

    typedef dmn_6<b,b,k_DCA,k_DCA,w1_dmn_t,w2_dmn_t> b_b_k_DCA_k_DCA_w_w_dmn_t;
    typedef dmn_6<b,b,k_PCM,k_PCM,w1_dmn_t,w2_dmn_t> b_b_k_PCM_k_PCM_w_w_dmn_t;

    FUNC_LIB::function<std::complex<scalar_type>, b_b_r_DCA_r_DCA_w_w_dmn_t> tmp_RR_1("tmp_RR_1");
    FUNC_LIB::function<std::complex<scalar_type>, b_b_r_PCM_r_PCM_w_w_dmn_t> tmp_RR_2("tmp_RR_2");

    FUNC_LIB::function<std::complex<scalar_type>, b_b_k_PCM_r_PCM_w_w_dmn_t> tmp_KR_2("tmp_KR_2");

    FUNC_LIB::function<std::complex<scalar_type>, b_b_k_DCA_k_DCA_w_w_dmn_t> tmp_KK_1("tmp_KK_1");
    FUNC_LIB::function<std::complex<scalar_type>, b_b_k_PCM_k_PCM_w_w_dmn_t> tmp_KK_2("tmp_KK_2");

    for(int l=0; l<M_r_r_w_w.size(); l++)
      tmp_RR_1(l) = M_r_r_w_w(l);

    for(int l=0; l<M_r_r_w_w.size(); l++)
      tmp_RR_2(l) = M_r_r_w_w(l);

    double time_1=0;
    {
      clock_t start = clock();

      FT<r_dmn_t, k_dmn_t>::execute(tmp_RR_1, tmp_KK_1);

      clock_t end = clock();

      //      cout << "\n\n\t time (1) : " << (end-start) << endl;

      time_1 = (end-start);
    }

    double time_2=0;
    {
      clock_t start = clock();

      math_algorithms::functional_transforms::TRANSFORM<r_PCM, k_PCM>::execute(tmp_RR_2, tmp_KR_2);
      math_algorithms::functional_transforms::TRANSFORM<r_PCM, k_PCM>::execute(tmp_KR_2, tmp_KK_2);

      //      math_algorithms::functional_transforms::TRANSFORM<r_PCM,
k_PCM>::execute_on_all(tmp_RR_2, tmp_KK_2);

      clock_t end = clock();

      //      cout << "\n\n\t time (2) : " << (end-start) << endl;

      time_2 = (end-start);
    }

    for(int l=0; l<M_r_r_w_w.size(); l++)
      if(abs(tmp_KK_1(l)-tmp_KK_2(l))>1.e-12)
        cout << abs(tmp_KK_1(l)-tmp_KK_2(l)) << endl;

    {
      static double N = 0;
      static double R = 0;

      N += 1;
      R += time_1/time_2;

      cout << "\n\tspeed-up : " << R/N << endl;
    }

//             for(int i=0; i<16; i++){
//             for(int j=0; j<16; j++)
//             cout << tmp_KR_1(0,0,i,j,15,16) << "\t";
//             cout << "\n";
//             }
//             cout << "\n";

//             for(int i=0; i<16; i++){
//             for(int j=0; j<16; j++)
//             cout << tmp_KK(0,0,i,j,15,16) << "\t";
//             cout << "\n";
//             }
//             cout << "\n";

//             for(int i=0; i<16; i++){
//             for(int j=0; j<16; j++)
//             cout << tmp_KR_1(0,0,i,j,15,16)-tmp_KK(0,0,i,j,15,16) << "\t";
//             cout << "\n";
//             }
//             cout << "\n";

//             for(int i=0; i<16; i++){
//             for(int j=0; j<16; j++)
//             cout << tmp_KR_1(0,0,i,j,15,16)/tmp_KK(0,0,i,j,15,16) << "\t";
//             cout << "\n";
//             }
//             cout << "\n";

//             assert(false);

  }
}
*/

template <class parameters_type, class MOMS_type>
void accumulator_nonlocal_G<parameters_type, MOMS_type>::compute_G2_k_k_w_w(
    FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& M_k_k_w_w,
    FUNC_LIB::function<std::complex<scalar_type>, b_b_k_k_w_w_dmn_t>& G2_k_k_w_w,
    e_spin_states_type e_spin) {
  std::complex<double> *G0_k_w_l_ptr, *G0_k_w_r_ptr;
  std::complex<scalar_type> *M_k_k_w_w_ptr, *G_k_k_w_w_ptr;

  // int spin_index = do_cast<int>::execute(e_spin);
  int spin_index = electron_spin_domain::to_coordinate(e_spin);

  for (int i = 0; i < G2_k_k_w_w.size(); i++)
    G2_k_k_w_w(i) = 0;

  for (int w2 = 0; w2 < w2_dmn_t::dmn_size(); w2++) {
    int corr_w2 = corresponding_w2_index[w2];

    for (int w1 = 0; w1 < w1_dmn_t::dmn_size(); w1++) {
      int corr_w1 = corresponding_w1_index[w1];

      for (int k2 = 0; k2 < k_dmn_t::dmn_size(); k2++) {
        for (int k1 = 0; k1 < k_dmn_t::dmn_size(); k1++) {
          G0_k_w_l_ptr = &MOMS.G0_k_w_cluster_excluded(0, spin_index, 0, spin_index, k1, corr_w1);
          M_k_k_w_w_ptr = &M_k_k_w_w(0, 0, k1, k2, w1, w2);
          G0_k_w_r_ptr = &MOMS.G0_k_w_cluster_excluded(0, spin_index, 0, spin_index, k2, corr_w2);

          G_k_k_w_w_ptr = &G2_k_k_w_w(0, 0, k1, k2, w1, w2);

          compute_G2_k_k_w_w_atomic(G0_k_w_l_ptr, M_k_k_w_w_ptr, G0_k_w_r_ptr, G_k_k_w_w_ptr);
        }
      }
    }
  }
}

template <class parameters_type, class MOMS_type>
inline void accumulator_nonlocal_G<parameters_type, MOMS_type>::compute_G2_k_k_w_w_atomic(
    std::complex<double>* G0_k_w_l_ptr, std::complex<scalar_type>* M_k_k_w_w_ptr,
    std::complex<double>* G0_k_w_r_ptr, std::complex<scalar_type>* G2_k_k_w_w_ptr) {
  switch (b_dmn_size) {
    case 1:
      compute_G2_k_k_w_w_atomic_1(G0_k_w_l_ptr, M_k_k_w_w_ptr, G0_k_w_r_ptr, G2_k_k_w_w_ptr);
      break;

    default:
      compute_G2_k_k_w_w_atomic_default(G0_k_w_l_ptr, M_k_k_w_w_ptr, G0_k_w_r_ptr, G2_k_k_w_w_ptr);
  }
}

template <class parameters_type, class MOMS_type>
inline void accumulator_nonlocal_G<parameters_type, MOMS_type>::compute_G2_k_k_w_w_atomic_1(
    std::complex<double>* G0_k_w_l_ptr, std::complex<scalar_type>* M_k_k_w_w_ptr,
    std::complex<double>* G0_k_w_r_ptr, std::complex<scalar_type>* G2_k_k_w_w_ptr) {
  G2_k_k_w_w_ptr[0] -= (std::complex<scalar_type>(G0_k_w_l_ptr[0]) * M_k_k_w_w_ptr[0] *
                        std::complex<scalar_type>(G0_k_w_r_ptr[0]));

  if (G0_k_w_l_ptr == G0_k_w_r_ptr)
    G2_k_k_w_w_ptr[0] += std::complex<scalar_type>(G0_k_w_l_ptr[0]) * beta;
}

template <class parameters_type, class MOMS_type>
inline void accumulator_nonlocal_G<parameters_type, MOMS_type>::compute_G2_k_k_w_w_atomic_default(
    std::complex<double>* G0_k_w_l_ptr, std::complex<scalar_type>* M_k_k_w_w_ptr,
    std::complex<double>* G0_k_w_r_ptr, std::complex<scalar_type>* G2_k_k_w_w_ptr) {
  for (int b1 = 0; b1 < b_dmn_size; b1++) {
    for (int b2 = 0; b2 < b_dmn_size; b2++) {
      for (int l1 = 0; l1 < b_dmn_size; l1++)
        for (int l2 = 0; l2 < b_dmn_size; l2++)
          G2_k_k_w_w_ptr[b1 + b2 * b_dmn_size] -=
              (std::complex<scalar_type>(G0_k_w_l_ptr[b1 + l1 * H_dmn_size]) *
               M_k_k_w_w_ptr[l1 + l2 * b_dmn_size] *
               std::complex<scalar_type>(G0_k_w_r_ptr[l2 + b2 * H_dmn_size]));

      if (G0_k_w_l_ptr == G0_k_w_r_ptr)
        G2_k_k_w_w_ptr[b1 + b2 * b_dmn_size] +=
            std::complex<scalar_type>(G0_k_w_l_ptr[b1 + b2 * H_dmn_size]) * beta;
    }
  }
}

}  // CT_AUX_ACCUMULATION
}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_ACCUMULATOR_TP_ACCUMULATOR_CTAUX_ACCUMULATOR_NONLOCAL_G_H
