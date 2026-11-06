#ifndef DCA_PHYS_PARAMETERS_MAIN_PARAMETERS_HPP
#define DCA_PHYS_PARAMETERS_MAIN_PARAMETERS_HPP

#include "dca/config/dca.hpp"

extern template class dca::phys::params::Parameters<
    Concurrency, Threading, Profiler, Model, RandomNumberGenerator, solver_name,
    dca::NumericalTraits<dca::config::McOptions::MC_REAL,
                         typename dca::util::ScalarSelect<dca::config::McOptions::MC_REAL,
                                                          Model::lattice_type::complex_g0>::type>>;

#endif
