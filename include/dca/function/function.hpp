// Copyright (C) 2009-2016 ETH Zurich
// Copyright (C) 2007?-2016 Center for Nanophase Materials Sciences, ORNL
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.txt for citation guidelines if you use this code for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//
// This class connects the function values to the domains.

#ifndef DCA_FUNCTION_FUNCTION_HPP
#define DCA_FUNCTION_FUNCTION_HPP

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "dca/function/copy_from.hpp"
#include "dca/function/scalar_cast.hpp"
#include "dca/function/set_to_zero.hpp"
#include "dca/util/type_utils.hpp"

namespace dca {
namespace func {
// dca::func::

template <typename scalartype, class domain>
class function {
public:
  typedef scalartype this_scalar_type;
  typedef domain this_domain_type;

  function();
  function(std::string name);
  function(const function<scalartype, domain>& other_one);
  // Destructive copy
  function(function<scalartype, domain>&& other_one);

  ~function();

  void reset();

  domain& get_domain() {
    return dmn;
  }
  const std::string& get_name() {
    return name_;
  }
  void set_name(const std::string& name) {
    name_ = name;
  }
  int signature() {
    return Nb_sbdms;
  }
  int size() {
    return Nb_elements;
  }
  // Returns the size of the subdomain with index 'index'.
  // Doesn't return function values!
  int operator[](int index) {
    return size_sbdm[index];
  }

  scalartype* values() {
    return fnc_values;
  }
  scalartype* values() const {
    return fnc_values;
  }

  void linind_2_subind(int i, int*& subind);
  void linind_2_subind(int linind, std::vector<int>& subind);

  void subind_2_linind(int* subind, int& i);
  template <typename T>
  int subind_2_linind(T i) {
    static_assert(std::is_integral<T>::value, "Index i must be an integer.");
    assert(i >= 0 && i < Nb_elements);
    return i;
  }
  template <typename... Ts>
  int subind_2_linind(Ts... indices) {
    // We need to cast all indices to the same type for dmn_variadic.
    return dmn(static_cast<int>(indices)...);
  }

  scalartype& operator()(int* subind);
  template <typename T>
  scalartype& operator()(T i) {
    static_assert(std::is_integral<T>::value, "Index i must be an integer.");
    assert(i >= 0 && i < Nb_elements);
    return fnc_values[i];
  }
  template <typename... Ts>
  scalartype& operator()(Ts... indices) {
    // We need to cast all indices to the same type for dmn_variadic.
    return fnc_values[dmn(static_cast<int>(indices)...)];
  }

  function<scalartype, domain>& operator=(function<scalartype, domain>& f_other);
  void operator+=(function<scalartype, domain>& f_other);
  void operator-=(function<scalartype, domain>& f_other);
  void operator*=(function<scalartype, domain>& f_other);
  void operator/=(function<scalartype, domain>& f_other);

  template <typename new_scalartype>
  void operator=(new_scalartype c);
  template <typename new_scalartype>
  void operator+=(new_scalartype c);
  template <typename new_scalartype>
  void operator-=(new_scalartype c);
  template <typename new_scalartype>
  void operator*=(new_scalartype c);
  template <typename new_scalartype>
  void operator/=(new_scalartype c);

  template <typename new_scalartype>
  void slice(int sbdm_index, int* subind, new_scalartype* fnc_vals);
  template <typename new_scalartype>
  void slice(int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals);
  template <typename new_scalartype>
  void distribute(int sbdm_index, int* subind, new_scalartype* fnc_vals);
  template <typename new_scalartype>
  void distribute(int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals);

  void print_fingerprint(std::ostream& stream);
  void print_fingerprint();
  void print_2_file(const char* file_name);

  template <typename concurrency_t>
  int get_buffer_size(const concurrency_t& concurrency);
  template <class concurrency_t>
  void pack(const concurrency_t& concurrency, int* buffer, int buffer_size, int& position);
  template <class concurrency_t>
  void unpack(const concurrency_t& concurrency, int* buffer, int buffer_size, int& position);

private:
  std::string name_;
  std::string function_type;

  domain dmn;
  int Nb_elements;

  int Nb_sbdms;
  std::vector<int>& size_sbdm;
  std::vector<int> step_sbdm;

  scalartype* fnc_values;
};

template <typename scalartype, class domain>
function<scalartype, domain>::function()
    : name_("no name"),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_elements(dmn.get_size()),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(Nb_sbdms, 1) {
  for (int i = 0; i < Nb_sbdms; i++)
    for (int j = 0; j < i; j++)
      step_sbdm[i] *= dmn.get_subdomain_size(j);

  fnc_values = new scalartype[Nb_elements];

  for (int linind = 0; linind < Nb_elements; linind++)
    set_to_zero::execute(fnc_values[linind]);
}

template <typename scalartype, class domain>
function<scalartype, domain>::function(std::string fnc_name)
    : name_(fnc_name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_elements(dmn.get_size()),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(Nb_sbdms, 1) {
  for (int i = 0; i < Nb_sbdms; i++)
    for (int j = 0; j < i; j++)
      step_sbdm[i] *= dmn.get_subdomain_size(j);

  fnc_values = new scalartype[Nb_elements];

  for (int linind = 0; linind < Nb_elements; linind++)
    set_to_zero::execute(fnc_values[linind]);
}

template <typename scalartype, class domain>
function<scalartype, domain>::function(const function<scalartype, domain>& other_one)
    : name_("no_name"),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_elements(dmn.get_size()),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(Nb_sbdms, 1) {
  for (int i = 0; i < Nb_sbdms; i++)
    for (int j = 0; j < i; j++)
      step_sbdm[i] *= dmn.get_subdomain_size(j);

  fnc_values = new scalartype[Nb_elements];

  copy_from<scalartype>::execute(Nb_elements, fnc_values, other_one.values());
}

template <typename scalartype, class domain>
function<scalartype, domain>::function(function<scalartype, domain>&& other_one)
    : name_("no_name"),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_elements(dmn.get_size()),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(Nb_sbdms, 1) {
  for (int i = 0; i < Nb_sbdms; i++)
    for (int j = 0; j < i; j++)
      step_sbdm[i] *= dmn.get_subdomain_size(j);

  fnc_values = other_one.fnc_values;
  other_one.fnc_values = nullptr;
}

template <typename scalartype, class domain>
function<scalartype, domain>::~function() {
  delete[] fnc_values;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::reset() {
  dmn.reset();

  for (int i = 0; i < Nb_sbdms; i++)
    size_sbdm[i] = dmn.get_subdomain_size(i);

  for (int i = 0; i < Nb_sbdms; i++) {
    step_sbdm[i] = 1;
    for (int j = 0; j < i; j++)
      step_sbdm[i] *= dmn.get_subdomain_size(j);
  }

  Nb_sbdms = dmn.get_leaf_domain_sizes().size();
  Nb_elements = dmn.get_size();

  delete[] fnc_values;
  fnc_values = new scalartype[Nb_elements];

  for (int linind = 0; linind < Nb_elements; linind++)
    set_to_zero::execute(fnc_values[linind]);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::linind_2_subind(int linind, int*& subind) {
  int tmp = linind;
  for (int i = 0; i < int(size_sbdm.size()); i++) {
    subind[i] = tmp % size_sbdm[i];
    tmp = (tmp - subind[i]) / size_sbdm[i];
  }
}

template <typename scalartype, class domain>
void function<scalartype, domain>::linind_2_subind(int linind, std::vector<int>& subind) {
  assert(int(subind.size()) == signature());

  int tmp = linind;
  for (int i = 0; i < int(size_sbdm.size()); i++) {
    subind[i] = tmp % size_sbdm[i];
    tmp = (tmp - subind[i]) / size_sbdm[i];
  }
}

template <typename scalartype, class domain>
void function<scalartype, domain>::subind_2_linind(int* subind, int& linind) {
  linind = 0;
  for (int i = 0; i < int(step_sbdm.size()); i++)
    linind += subind[i] * step_sbdm[i];
}

template <typename scalartype, class domain>
scalartype& function<scalartype, domain>::operator()(int* subind) {
  int linind;
  subind_2_linind(subind, linind);

  assert(linind >= 0 && linind < Nb_elements);
  return fnc_values[linind];
}

template <typename scalartype, class domain>
function<scalartype, domain>& function<scalartype, domain>::operator=(
    function<scalartype, domain>& f_other) {
  domain& dmn_other = f_other.get_domain();

  if (dmn.get_size() !=
      dmn_other.get_size())  // Domains were not initialized when function was created.
  {
    dmn.get_size() = dmn_other.get_size();
    dmn.get_branch_domain_sizes() = dmn_other.get_branch_domain_sizes();
    dmn.get_leaf_domain_sizes() = dmn_other.get_leaf_domain_sizes();

    for (int i = 0; i < Nb_sbdms; i++)
      size_sbdm[i] = dmn.get_subdomain_size(i);

    for (int i = 0; i < Nb_sbdms; i++)
      for (int j = 0; j < i; j++)
        step_sbdm[i] *= dmn.get_subdomain_size(j);

    Nb_sbdms = dmn.get_leaf_domain_sizes().size();
    Nb_elements = dmn.get_size();

    delete[] fnc_values;
    fnc_values = new scalartype[Nb_elements];
  }

  memcpy(fnc_values, f_other.values(), Nb_elements * sizeof(scalartype));

  return *this;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator+=(function<scalartype, domain>& f_other) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] += f_other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator-=(function<scalartype, domain>& f_other) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] -= f_other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator*=(function<scalartype, domain>& f_other) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] *= f_other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator/=(function<scalartype, domain>& f_other) {
  for (int linind = 0; linind < Nb_elements; linind++) {
    assert(ASSERT_NON_ZERO(f_other(linind)));
    fnc_values[linind] /= f_other(linind);
  }
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::operator=(new_scalartype c) {
  scalartype c_new(c);

  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] = c_new;
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::operator+=(new_scalartype c) {
  scalartype c_new(c);

  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] += c_new;
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::operator-=(new_scalartype c) {
  scalartype c_new(c);

  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] -= c_new;
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::operator*=(new_scalartype c) {
  scalartype c_new(c);

  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] *= c_new;
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::operator/=(new_scalartype c) {
  scalartype c_new(c);

  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] /= c_new;
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::slice(int sbdm_index, int* subind, new_scalartype* fnc_vals) {
  assert(sbdm_index >= 0);
  assert(sbdm_index < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index] = 0;
  subind_2_linind(subind, linind);

  for (int i = 0; i < size_sbdm[sbdm_index]; i++)
    fnc_vals[i] = ScalarCast<new_scalartype>::execute(fnc_values[linind + i * step_sbdm[sbdm_index]]);
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::slice(int sbdm_index_1, int sbdm_index_2, int* subind,
                                         new_scalartype* fnc_vals) {
  assert(sbdm_index_1 >= 0);
  assert(sbdm_index_2 >= 0);
  assert(sbdm_index_1 < Nb_sbdms);
  assert(sbdm_index_2 < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index_1] = 0;
  subind[sbdm_index_2] = 0;
  subind_2_linind(subind, linind);

  int size_sbdm_1 = size_sbdm[sbdm_index_1];
  int size_sbdm_2 = size_sbdm[sbdm_index_2];

  int step_sbdm_1 = step_sbdm[sbdm_index_1];
  int step_sbdm_2 = step_sbdm[sbdm_index_2];

  new_scalartype* fnc_ptr_left = NULL;
  new_scalartype* fnc_ptr_right = NULL;

  for (int j = 0; j < size_sbdm_2; j++) {
    fnc_ptr_left = &fnc_vals[0 + j * size_sbdm_1];
    fnc_ptr_right = &fnc_values[linind + j * step_sbdm_2];

    for (int i = 0; i < size_sbdm_1; i++)
      fnc_ptr_left[i] = fnc_ptr_right[i * step_sbdm_1];
    //       fnc_vals[i+j*size_sbdm[sbdm_index_1]] = fnc_values[linind + i*step_sbdm[sbdm_index_1] +
    //       j*step_sbdm[sbdm_index_2]];
  }
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::distribute(int sbdm_index, int* subind, new_scalartype* fnc_vals) {
  assert(sbdm_index >= 0);
  assert(sbdm_index < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index] = 0;
  subind_2_linind(subind, linind);

  for (int i = 0; i < size_sbdm[sbdm_index]; i++)
    fnc_values[linind + i * step_sbdm[sbdm_index]] = ScalarCast<scalartype>::execute(fnc_vals[i]);
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::distribute(int sbdm_index_1, int sbdm_index_2, int* subind,
                                              new_scalartype* fnc_vals) {
  assert(sbdm_index_1 >= 0);
  assert(sbdm_index_2 >= 0);
  assert(sbdm_index_1 < Nb_sbdms);
  assert(sbdm_index_2 < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index_1] = 0;
  subind[sbdm_index_2] = 0;
  subind_2_linind(subind, linind);

  for (int i = 0; i < size_sbdm[sbdm_index_1]; i++)
    for (int j = 0; j < size_sbdm[sbdm_index_2]; j++)
      fnc_values[linind + i * step_sbdm[sbdm_index_1] + j * step_sbdm[sbdm_index_2]] =
          fnc_vals[i + j * size_sbdm[sbdm_index_1]];
}

template <typename scalartype, class domain>
void function<scalartype, domain>::print_fingerprint(std::ostream& stream) {
  stream << std::endl << std::endl << "function : " << name_ << std::endl;

  stream << "*********************************" << std::endl;

  stream << "# subdomains        : " << Nb_sbdms << std::endl;

  dca::util::print_type<domain>::print(stream);

  stream << "size of subdomains  : " << std::endl;
  for (int i = 0; i < Nb_sbdms; i++)
    stream << size_sbdm[i] << "\t";
  stream << std::endl;

  stream << "memory step         : " << std::endl;
  for (int i = 0; i < Nb_sbdms; i++)
    stream << step_sbdm[i] << "\t";
  stream << std::endl;

  stream << "# elements          : " << Nb_elements << std::endl;
  stream << "# size              : " << Nb_elements * sizeof(scalartype) * (1.e-6)
         << " (mega-bytes)" << std::endl;
  stream << "*********************************" << std::endl;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::print_fingerprint() {
  print_fingerprint(std::cout);
}

template <typename scalartype, class domain>
template <typename concurrency_t>
int function<scalartype, domain>::get_buffer_size(const concurrency_t& concurrency) {
  int result = 0;
  result += concurrency.get_buffer_size(*this);
  return result;
}

template <typename scalartype, class domain>
template <class concurrency_t>
void function<scalartype, domain>::pack(const concurrency_t& concurrency, int* buffer,
                                        int buffer_size, int& position) {
  concurrency.pack(buffer, buffer_size, position, *this);
}

template <typename scalartype, class domain>
template <class concurrency_t>
void function<scalartype, domain>::unpack(const concurrency_t& concurrency, int* buffer,
                                          int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, *this);
}

}  // func
}  // dca

#endif  // DCA_FUNCTION_FUNCTION_HPP
