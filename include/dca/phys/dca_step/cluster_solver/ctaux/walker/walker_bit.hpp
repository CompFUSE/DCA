// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// Built-in test (BIT) for CT-AUX walker.

#ifndef DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_WALKER_BIT_HPP
#define DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_WALKER_BIT_HPP

#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/linalg/matrix.hpp"
#include "dca/linalg/matrixop.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/ct_aux_hs_configuration.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/cv.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/vertex_singleton.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/ctaux_typedefs.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g0_interpolation/g0_interpolation.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/g_tools.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/walker/tools/n_tools.hpp"
#include "dca/phys/domains/quantum/numerical_error_domain.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctaux {
// dca::phys::solver::ctaux::

template <class Parameters, class Data, typename Real>
class WalkerBIT {
  typedef vertex_singleton vertex_singleton_type;
  typedef CT_AUX_HS_configuration<Parameters> configuration_type;

  using rng_type = typename Parameters::random_number_generator;

  typedef typename CtauxTypedefs<Parameters, Data>::profiler_type profiler_type;
  typedef typename CtauxTypedefs<Parameters, Data>::concurrency_type concurrency_type;

public:
  WalkerBIT(Parameters& parameters_ref, Data& MOMS_ref, int id);
  ~WalkerBIT();

  void initialize();

  func::function<Real, func::dmn_0<domains::numerical_error_domain>>& get_error_distribution();

  template <dca::linalg::DeviceType device_t>
  void check_G0_matrices(configuration_type& configuration,
                         dca::linalg::Matrix<Real, device_t>& G0_up,
                         dca::linalg::Matrix<Real, device_t>& G0_dn);

  template <dca::linalg::DeviceType device_t>
  void check_N_matrices(configuration_type& configuration, dca::linalg::Matrix<Real, device_t>& G0_up,
                        dca::linalg::Matrix<Real, device_t>& G0_dn,
                        dca::linalg::Matrix<Real, device_t>& N_up,
                        dca::linalg::Matrix<Real, device_t>& N_dn);

  template <dca::linalg::DeviceType device_t>
  void check_G_matrices(configuration_type& configuration, dca::linalg::Matrix<Real, device_t>& G0_up,
                        dca::linalg::Matrix<Real, device_t>& G0_dn,
                        dca::linalg::Matrix<Real, device_t>& N_up,
                        dca::linalg::Matrix<Real, device_t>& N_dn,
                        dca::linalg::Matrix<Real, device_t>& G_up,
                        dca::linalg::Matrix<Real, device_t>& G_dn);

private:
  Parameters& parameters;
  Data& MOMS;
  concurrency_type& concurrency;

  int thread_id;

  CV<Parameters> CV_obj;

  G0_INTERPOLATION<dca::linalg::CPU, Parameters, Real> G0_CPU_tools_obj;
  N_TOOLS<dca::linalg::CPU, Parameters, Real> N_CPU_tools_obj;
  G_TOOLS<dca::linalg::CPU, Parameters, Real> G_CPU_tools_obj;

  dca::linalg::Matrix<Real, dca::linalg::CPU> G0_up_CPU;
  dca::linalg::Matrix<Real, dca::linalg::CPU> G0_dn_CPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> N_up_CPU;
  dca::linalg::Matrix<Real, dca::linalg::CPU> N_dn_CPU;

  dca::linalg::Matrix<Real, dca::linalg::CPU> G_up_CPU;
  dca::linalg::Matrix<Real, dca::linalg::CPU> G_dn_CPU;

  func::function<Real, func::dmn_0<domains::numerical_error_domain>> error;
};

template <class Parameters, class Data, typename Real>
WalkerBIT<Parameters, Data, Real>::WalkerBIT(Parameters& parameters_ref, Data& MOMS_ref, int id)
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

template <class Parameters, class Data, typename Real>
WalkerBIT<Parameters, Data, Real>::~WalkerBIT() {
  //     for(int l=0; l<numerical_error_domain::get_size(); l++)
  //       cout << numerical_error_domain::get_elements()[l] << "\t" << error(l) << endl;
}

template <class Parameters, class Data, typename Real>
func::function<Real, func::dmn_0<domains::numerical_error_domain>>& WalkerBIT<
    Parameters, Data, Real>::get_error_distribution() {
  return error;
}

template <class Parameters, class Data, typename Real>
void WalkerBIT<Parameters, Data, Real>::initialize() {
  error = 0.;

  CV_obj.initialize(MOMS);
  G0_CPU_tools_obj.initialize(MOMS);

  //     for(int i=-16; i<=0; i++){
  //       for(Real j=1; j<10; j+=1.){
  //        x.push_back(j*std::pow(10., i));
  //        y.push_back(0);
  //       }
  //     }
}

template <class Parameters, class Data, typename Real>
template <dca::linalg::DeviceType device_t>
void WalkerBIT<Parameters, Data, Real>::check_G0_matrices(configuration_type& configuration,
                                                          dca::linalg::Matrix<Real, device_t>& G0_up,
                                                          dca::linalg::Matrix<Real, device_t>& G0_dn) {
  //     cout << __FUNCTION__ << endl;

  G0_CPU_tools_obj.build_G0_matrix(configuration, G0_up_CPU, e_UP);
  G0_CPU_tools_obj.build_G0_matrix(configuration, G0_dn_CPU, e_DN);

  linalg::matrixop::difference(G0_up_CPU, G0_up);
  linalg::matrixop::difference(G0_dn_CPU, G0_dn);
}

template <class Parameters, class Data, typename Real>
template <dca::linalg::DeviceType device_t>
void WalkerBIT<Parameters, Data, Real>::check_N_matrices(configuration_type& configuration,
                                                         dca::linalg::Matrix<Real, device_t>& G0_up,
                                                         dca::linalg::Matrix<Real, device_t>& G0_dn,
                                                         dca::linalg::Matrix<Real, device_t>& N_up,
                                                         dca::linalg::Matrix<Real, device_t>& N_dn) {
  //     cout << __FUNCTION__ << endl;

  linalg::matrixop::difference(G0_up_CPU, G0_up);
  linalg::matrixop::difference(G0_dn_CPU, G0_dn);

  N_CPU_tools_obj.build_N_matrix(configuration, N_up_CPU, G0_up_CPU, e_UP);
  N_CPU_tools_obj.build_N_matrix(configuration, N_dn_CPU, G0_dn_CPU, e_DN);

  Real err_up = linalg::matrixop::difference(N_up_CPU, N_up);
  Real err_dn = linalg::matrixop::difference(N_dn_CPU, N_dn);

  std::vector<Real>& x = domains::numerical_error_domain::get_elements();
  for (size_t l = 0; l < x.size() - 1; l++)
    if ((err_up > x[l] and err_up < x[l + 1]) or (err_dn > x[l] and err_dn < x[l + 1]))
      error(l) += 1;
}

template <class Parameters, class Data, typename Real>
template <dca::linalg::DeviceType device_t>
void WalkerBIT<Parameters, Data, Real>::check_G_matrices(configuration_type& configuration,
                                                         dca::linalg::Matrix<Real, device_t>& G0_up,
                                                         dca::linalg::Matrix<Real, device_t>& G0_dn,
                                                         dca::linalg::Matrix<Real, device_t>& N_up,
                                                         dca::linalg::Matrix<Real, device_t>& N_dn,
                                                         dca::linalg::Matrix<Real, device_t>& G_up,
                                                         dca::linalg::Matrix<Real, device_t>& G_dn) {
  //     cout << __FUNCTION__ << endl;

  linalg::matrixop::difference(G0_up_CPU, G0_up);
  linalg::matrixop::difference(G0_dn_CPU, G0_dn);

  linalg::matrixop::difference(N_up_CPU, N_up);
  linalg::matrixop::difference(N_dn_CPU, N_dn);

  G_CPU_tools_obj.build_G_matrix(configuration, N_up_CPU, G0_up_CPU, G_up_CPU, e_UP);
  G_CPU_tools_obj.build_G_matrix(configuration, N_dn_CPU, G0_dn_CPU, G_dn_CPU, e_DN);

  linalg::matrixop::difference(G_up_CPU, G_up);
  linalg::matrixop::difference(G_dn_CPU, G_dn);
}

}  // namespace ctaux
}  // namespace solver
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_DCA_STEP_CLUSTER_SOLVER_CTAUX_WALKER_WALKER_BIT_HPP
