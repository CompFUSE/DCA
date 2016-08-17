// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Description

#ifndef PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_BUILD_IN_TEST_H
#define PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_BUILD_IN_TEST_H

#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_template/mc_walker_bit.hpp"

#include <vector>

#include "dca/linalg/matrix.hpp"

#include "comp_library/function_library/include_function_library.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_auxilery_field_coefficients.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_hubbard_stratonovitch_configuration.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_structs/ctaux_vertex_singleton.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_typedefinitions.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G0_matrix_routines/ctaux_G0_matrix_routines_CPU.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_G_matrix_routines.h"
#include "phys_library/DCA+_step/cluster_solver/cluster_solver_mc_ctaux/ctaux_walker/ctaux_walker_tools/ctaux_N_matrix_routines.h"
#include "phys_library/domains/Quantum_domain/numerical_error_domain.h"

namespace DCA {
namespace QMCI {
// DCA::QMCI::

template <class parameters_type, class MOMS_type>
class MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type> {
  typedef vertex_singleton vertex_singleton_type;
  typedef CT_AUX_HS_configuration<parameters_type> configuration_type;

  using rng_type = typename parameters_type::random_number_generator;

  typedef
      typename MC_type_definitions<CT_AUX_SOLVER, parameters_type, MOMS_type>::profiler_type profiler_type;
  typedef typename MC_type_definitions<CT_AUX_SOLVER, parameters_type, MOMS_type>::concurrency_type
      concurrency_type;

public:
  MC_walker_BIT(parameters_type& parameters_ref, MOMS_type& MOMS_ref, int id);

  ~MC_walker_BIT();

  void initialize();

  FUNC_LIB::function<double, dmn_0<numerical_error_domain>>& get_error_distribution();

  template <dca::linalg::DeviceType device_t>
  void check_G0_matrices(configuration_type& configuration,
                         dca::linalg::Matrix<double, device_t>& G0_up,
                         dca::linalg::Matrix<double, device_t>& G0_dn);

  template <dca::linalg::DeviceType device_t>
  void check_N_matrices(configuration_type& configuration,
                        dca::linalg::Matrix<double, device_t>& G0_up,
                        dca::linalg::Matrix<double, device_t>& G0_dn,
                        dca::linalg::Matrix<double, device_t>& N_up,
                        dca::linalg::Matrix<double, device_t>& N_dn);

  template <dca::linalg::DeviceType device_t>
  void check_G_matrices(configuration_type& configuration,
                        dca::linalg::Matrix<double, device_t>& G0_up,
                        dca::linalg::Matrix<double, device_t>& G0_dn,
                        dca::linalg::Matrix<double, device_t>& N_up,
                        dca::linalg::Matrix<double, device_t>& N_dn,
                        dca::linalg::Matrix<double, device_t>& G_up,
                        dca::linalg::Matrix<double, device_t>& G_dn);

private:
  parameters_type& parameters;
  MOMS_type& MOMS;
  concurrency_type& concurrency;

  int thread_id;

  CV<parameters_type> CV_obj;

  G0_INTERPOLATION<dca::linalg::CPU, parameters_type> G0_CPU_tools_obj;
  N_TOOLS<dca::linalg::CPU, parameters_type> N_CPU_tools_obj;
  G_TOOLS<dca::linalg::CPU, parameters_type> G_CPU_tools_obj;

  dca::linalg::Matrix<double, dca::linalg::CPU> G0_up_CPU;
  dca::linalg::Matrix<double, dca::linalg::CPU> G0_dn_CPU;

  dca::linalg::Matrix<double, dca::linalg::CPU> N_up_CPU;
  dca::linalg::Matrix<double, dca::linalg::CPU> N_dn_CPU;

  dca::linalg::Matrix<double, dca::linalg::CPU> G_up_CPU;
  dca::linalg::Matrix<double, dca::linalg::CPU> G_dn_CPU;

  FUNC_LIB::function<double, dmn_0<numerical_error_domain>> error;
};

template <class parameters_type, class MOMS_type>
MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::MC_walker_BIT(parameters_type& parameters_ref,
                                                                        MOMS_type& MOMS_ref, int id)
    : parameters(parameters_ref),
      MOMS(MOMS_ref),
      concurrency(parameters.get_concurrency()),

      thread_id(id),

      CV_obj(parameters),

      G0_CPU_tools_obj(thread_id, parameters),
      N_CPU_tools_obj(thread_id, parameters, CV_obj),
      G_CPU_tools_obj(thread_id, parameters, CV_obj),

      G0_up_CPU("G0_up_CPU (MC_walker_BIT)"),
      G0_dn_CPU("G0_up_CPU (MC_walker_BIT)"),

      N_up_CPU("N_up_CPU (MC_walker_BIT)"),
      N_dn_CPU("N_up_CPU (MC_walker_BIT)"),

      G_up_CPU("G_up_CPU (MC_walker_BIT)"),
      G_dn_CPU("G_up_CPU (MC_walker_BIT)"),

      error("error") {}

template <class parameters_type, class MOMS_type>
MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::~MC_walker_BIT() {
  //     for(int l=0; l<numerical_error_domain::get_size(); l++)
  //       cout << numerical_error_domain::get_elements()[l] << "\t" << error(l) << endl;
}

template <class parameters_type, class MOMS_type>
FUNC_LIB::function<double, dmn_0<numerical_error_domain>>& MC_walker_BIT<
    CT_AUX_SOLVER, parameters_type, MOMS_type>::get_error_distribution() {
  return error;
}

template <class parameters_type, class MOMS_type>
void MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::initialize() {
  error = 0.;

  CV_obj.initialize(MOMS);
  G0_CPU_tools_obj.initialize(MOMS);

  //     for(int i=-16; i<=0; i++){
  //       for(double j=1; j<10; j+=1.){
  //        x.push_back(j*std::pow(10., i));
  //        y.push_back(0);
  //       }
  //     }
}

template <class parameters_type, class MOMS_type>
template <dca::linalg::DeviceType device_t>
void MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::check_G0_matrices(
    configuration_type& configuration, dca::linalg::Matrix<double, device_t>& G0_up,
    dca::linalg::Matrix<double, device_t>& G0_dn) {
  //     cout << __FUNCTION__ << endl;

  G0_CPU_tools_obj.build_G0_matrix(configuration, G0_up_CPU, e_UP);
  G0_CPU_tools_obj.build_G0_matrix(configuration, G0_dn_CPU, e_DN);

  G0_up_CPU.difference(G0_up);
  G0_dn_CPU.difference(G0_dn);
}

template <class parameters_type, class MOMS_type>
template <dca::linalg::DeviceType device_t>
void MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::check_N_matrices(
    configuration_type& configuration, dca::linalg::Matrix<double, device_t>& G0_up,
    dca::linalg::Matrix<double, device_t>& G0_dn, dca::linalg::Matrix<double, device_t>& N_up,
    dca::linalg::Matrix<double, device_t>& N_dn) {
  //     cout << __FUNCTION__ << endl;

  G0_up_CPU.difference(G0_up);
  G0_dn_CPU.difference(G0_dn);

  N_CPU_tools_obj.build_N_matrix(configuration, N_up_CPU, G0_up_CPU, e_UP);
  N_CPU_tools_obj.build_N_matrix(configuration, N_dn_CPU, G0_dn_CPU, e_DN);

  double err_up = N_up_CPU.difference(N_up);
  double err_dn = N_dn_CPU.difference(N_dn);

  std::vector<double>& x = numerical_error_domain::get_elements();
  for (size_t l = 0; l < x.size() - 1; l++)
    if ((err_up > x[l] and err_up < x[l + 1]) or (err_dn > x[l] and err_dn < x[l + 1]))
      error(l) += 1;
}

template <class parameters_type, class MOMS_type>
template <dca::linalg::DeviceType device_t>
void MC_walker_BIT<CT_AUX_SOLVER, parameters_type, MOMS_type>::check_G_matrices(
    configuration_type& configuration, dca::linalg::Matrix<double, device_t>& G0_up,
    dca::linalg::Matrix<double, device_t>& G0_dn, dca::linalg::Matrix<double, device_t>& N_up,
    dca::linalg::Matrix<double, device_t>& N_dn, dca::linalg::Matrix<double, device_t>& G_up,
    dca::linalg::Matrix<double, device_t>& G_dn) {
  //     cout << __FUNCTION__ << endl;

  G0_up_CPU.difference(G0_up);
  G0_dn_CPU.difference(G0_dn);

  N_up_CPU.difference(N_up);
  N_dn_CPU.difference(N_dn);

  G_CPU_tools_obj.build_G_matrix(configuration, N_up_CPU, G0_up_CPU, G_up_CPU, e_UP);
  G_CPU_tools_obj.build_G_matrix(configuration, N_dn_CPU, G0_dn_CPU, G_dn_CPU, e_DN);

  G_up_CPU.difference(G_up);
  G_dn_CPU.difference(G_dn);
}

}  // QMCI
}  // DCA

#endif  // PHYS_LIBRARY_DCA_STEP_CLUSTER_SOLVER_CLUSTER_SOLVER_MC_CTAUX_CTAUX_WALKER_CTAUX_WALKER_BUILD_IN_TEST_H
