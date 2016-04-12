//-*-C++-*-

#ifndef PACKING_INTERFACE_HEADER
#define PACKING_INTERFACE_HEADER

#include <vector>
#include "comp_library/function_library/function.h"
namespace COMP_LIB {
/*!
 *  \author Peter Staar
 */
template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
class packing_interface {
public:
  packing_interface(processor_grouping<LIBRARY>& grouping_ref);
  ~packing_interface();

  /************************************
   ***  size
   ************************************/

  template <typename T>
  size_t get_buffer_size(T& item);

  size_t get_buffer_size(std::string str);

  template <typename T>
  size_t get_buffer_size(std::vector<T>& v);

  template <typename T, class dmn_type>
  size_t get_buffer_size(FUNC_LIB::function<T, dmn_type>& f);

  // Note: The packing and unpacking implementations are not provided by this class.

private:
  processor_grouping<LIBRARY>& grouping;
};

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
packing_interface<LIBRARY>::packing_interface(processor_grouping<LIBRARY>& grouping_ref)
    : grouping(grouping_ref) {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
packing_interface<LIBRARY>::~packing_interface() {}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename T>
size_t packing_interface<LIBRARY>::get_buffer_size(T& /*item*/) {
  return sizeof(T);
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
size_t packing_interface<LIBRARY>::get_buffer_size(std::string str) {
  return str.size() * sizeof(char);
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename T>
size_t packing_interface<LIBRARY>::get_buffer_size(std::vector<T>& v) {
  return v.size() * size(v[0]);
}

template <PARALLELIZATION_LIBRARY_NAMES LIBRARY>
template <typename T, class dmn_type>
size_t packing_interface<LIBRARY>::get_buffer_size(FUNC_LIB::function<T, dmn_type>& f) {
  return f.size() * size(f(0));
}
}

#endif
