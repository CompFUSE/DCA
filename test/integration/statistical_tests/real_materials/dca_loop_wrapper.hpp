// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// wrapper for accessing protected members of DCA_loop

#include "dca/phys/dca_loop/dca_loop.hpp"

#ifndef TEST_INTEGRATION_STATISTICAL_TESTS_REAL_MATERIALS_DCA_LOOP_WRAPPER_HPP
#define TEST_INTEGRATION_STATISTICAL_TESTS_REAL_MATERIALS_DCA_LOOP_WRAPPER_HPP

namespace dca {
namespace testing {
// dca::testing
template <class ParametersType, class DataType, class MonteCarloIntegrator>
class DcaLoopWrapper : public dca::phys::DcaLoop<ParametersType, DataType, MonteCarloIntegrator> {
private:
  int step_;

public:
  using BaseLoop = phys::DcaLoop<ParametersType, DataType, MonteCarloIntegrator>;
  using Concurrency = typename ParametersType::concurrency_type;

  DcaLoopWrapper(ParametersType& parameters_ref, DataType& MOMS_ref, Concurrency& concurrency_ref)
      : BaseLoop(parameters_ref, MOMS_ref, concurrency_ref),
        step_(BaseLoop::parameters.get_dca_iterations()) {}

  MonteCarloIntegrator& getSolver() {
    return BaseLoop::monte_carlo_integrator_;
  }

  void performPreIntegrationSteps() {
    BaseLoop::initialize();
    BaseLoop::adjust_chemical_potential();
    BaseLoop::perform_cluster_mapping();
    BaseLoop::adjust_coarsegrained_self_energy();
    BaseLoop::perform_cluster_exclusion_step();
  }

  void performIntegrationStep() {
    BaseLoop::monte_carlo_integrator_.initialize(step_);
    BaseLoop::monte_carlo_integrator_.integrate();
  }

  void performPostIntegrationStep() {
    BaseLoop::monte_carlo_integrator_.finalize();
    BaseLoop::adjust_impurity_self_energy();
    BaseLoop::perform_lattice_mapping();
    BaseLoop::update_DCA_loop_data_functions(step_);
  }

  void write(){
    BaseLoop::write();
    dca::io::HDF5Writer writer;
    writer.open_file(BaseLoop::parameters.get_filename_dca(), false);
    writer.open_group("additional_functions");
    writer.execute(BaseLoop::MOMS.Sigma_cluster);
    writer.close_group();
    writer.close_file();
  }
};

}  // testing
}  // dca
#endif  // TEST_INTEGRATION_STATISTICAL_TESTS_REAL_MATERIALS_DCA_LOOP_WRAPPER_HPP
