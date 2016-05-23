//-*-C++-*-

#ifndef COLLECTIVE_SUM_INTERFACE_H
#define COLLECTIVE_SUM_INTERFACE_H
#include "map"
#include <vector>
#include "comp_library/linalg/src/vector.h"
#include "comp_library/linalg/src/matrix.h"
namespace COMP_LIB {
/*!
 *  \author Peter Staar
 */
template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
class collective_sum_interface {
public:
  collective_sum_interface(processor_grouping<LIBRARY>& grouping_ref);
  ~collective_sum_interface();

  template <typename scalar_type>
  void sum(scalar_type& value);

  template <typename scalar_type>
  void sum(std::vector<scalar_type>& m);

  template <typename scalartype>
  void sum(std::map<std::string, std::vector<scalartype>>& m);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<scalar_type, domain>& f);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<scalar_type, domain>& f,
           FUNC_LIB::function<scalar_type, domain>& f_target);

  template <typename scalar_type, class domain>
  void sum(FUNC_LIB::function<std::vector<scalar_type>, domain>& f);

  template <typename scalar_type>
  void sum(LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& f);

  template <typename scalar_type>
  void sum(LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& f);

  template <typename some_type>
  void sum_and_average(some_type& obj, int size);

private:
  processor_grouping<LIBRARY>& grouping;
};

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
collective_sum_interface<LIBRARY>::collective_sum_interface(processor_grouping<LIBRARY>& grouping_ref)
    : grouping(grouping_ref) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
collective_sum_interface<LIBRARY>::~collective_sum_interface() {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(scalar_type& /*value*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(std::vector<scalar_type>& /*m*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(std::map<std::string, std::vector<scalar_type>>& /*m*/) {
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type, class domain>
void collective_sum_interface<LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& /*f*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type, class domain>
void collective_sum_interface<LIBRARY>::sum(FUNC_LIB::function<scalar_type, domain>& /*f*/,
                                            FUNC_LIB::function<scalar_type, domain>& /*f_target*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type, class domain>
void collective_sum_interface<LIBRARY>::sum(FUNC_LIB::function<std::vector<scalar_type>, domain>& /*f*/) {
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(LIN_ALG::vector<scalar_type, LIN_ALG::CPU>& /*f*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename scalar_type>
void collective_sum_interface<LIBRARY>::sum(LIN_ALG::matrix<scalar_type, LIN_ALG::CPU>& /*f*/) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename some_type>
void collective_sum_interface<LIBRARY>::sum_and_average(some_type& /*obj*/, int /*size*/) {}
}

#endif
