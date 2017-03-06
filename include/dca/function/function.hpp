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

#include <algorithm>
#include <cassert>
#include <cmath>    // for std::abs
#include <complex>  // for std::abs(std::complex)
#include <iostream>
#include <string>
#include <vector>

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
  function(const std::string& name);
  function(const function<scalartype, domain>& other_one);
  // Destructive copy
  function(function<scalartype, domain>&& other_one);

  ~function();

  void reset();

  const domain& get_domain() const {
    return dmn;
  }
  const std::string& get_name() const {
    return name_;
  }
  void set_name(const std::string& name) {
    name_ = name;
  }
  int signature() const {
    return Nb_sbdms;
  }
  int size() const {
    return Nb_elements;
  }
  // Returns the size of the subdomain with index 'index'.
  // Doesn't return function values!
  int operator[](int index) const {
    return size_sbdm[index];
  }

  scalartype* values() {
    return fnc_values;
  }
  const scalartype* values() const {
    return fnc_values;
  }

  void linind_2_subind(int i, int*& subind) const;
  void linind_2_subind(int linind, std::vector<int>& subind) const;

  void subind_2_linind(const int* subind, int& i) const;

  template <typename T>
  int subind_2_linind(T i) const {
    static_assert(std::is_integral<T>::value, "Index i must be an integer.");
    assert(i >= 0 && i < Nb_elements);
    return i;
  }
  template <typename... Ts>
  // Enable only if all arguments are integer.
  // INTERNAL: This is done to prevent subind_to_linind(int*, int) to resolve to
  //           subind_to_linind(int...) rather than subind_to_linind(const int*, int).
  std::enable_if_t<dca::util::if_all<std::is_integral<Ts>::value...>::value, int> subind_2_linind(
      Ts... indices) const {
    // We need to cast all indices to the same type for dmn_variadic.
    return dmn(static_cast<int>(indices)...);
  }

  scalartype& operator()(const int* subind);
  const scalartype& operator()(const int* subind) const;

  template <typename T>
  scalartype& operator()(const T i) {
    static_assert(std::is_integral<T>::value, "Index i must be an integer.");
    assert(i >= 0 && i < Nb_elements);
    return fnc_values[i];
  }

  template <typename T>
  const scalartype& operator()(const T i) const {
    static_assert(std::is_integral<T>::value, "Index i must be an integer.");
    assert(i >= 0 && i < Nb_elements);
    return fnc_values[i];
  }

  template <typename... Ts>
  scalartype& operator()(Ts... indices) {
    // We need to cast all indices to the same type for dmn_variadic.
    return fnc_values[dmn(static_cast<int>(indices)...)];
  }

  template <typename... Ts>
  const scalartype& operator()(Ts... indices) const {
    return fnc_values[dmn(static_cast<int>(indices)...)];
  }

  // Assignment operators copy (move) the function values, but not the name.
  function<scalartype, domain>& operator=(const function<scalartype, domain>& f_other);
  function<scalartype, domain>& operator=(function<scalartype, domain>&& f_other);

  void operator+=(const function<scalartype, domain>& f_other);
  void operator-=(const function<scalartype, domain>& f_other);
  void operator*=(const function<scalartype, domain>& f_other);
  void operator/=(const function<scalartype, domain>& f_other);

  template <typename new_scalartype>
  void operator=(new_scalartype c);
  template <typename new_scalartype>
  void operator+=(new_scalartype c);
  template <typename new_scalartype>
  void operator-=(new_scalartype c);
  template <typename new_scalartype>
  void operator*=(new_scalartype c);
  template <typename new_scalartype>
  void operator/=(new_scalartype );

  bool operator==(const function<scalartype, domain>& f_other) const;

  template <typename new_scalartype>
  void slice(int sbdm_index, int* subind, new_scalartype* fnc_vals) const;
  template <typename new_scalartype>
  void slice(int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals) const;
  template <typename new_scalartype>
  void distribute(int sbdm_index, int* subind, const new_scalartype* fnc_vals);
  template <typename new_scalartype>
  void distribute(int sbdm_index_1, int sbdm_index_2, int* subind, const new_scalartype* fnc_vals);

  // Prints the function's metadata.
  void print_fingerprint(std::ostream& stream = std::cout) const;
  // Prints the function's elements.
  void print_elements(std::ostream& stream = std::cout) const;

  template <typename concurrency_t>
  int get_buffer_size(const concurrency_t& concurrency) const;
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
  const std::vector<int>& size_sbdm;
  const std::vector<int>& step_sbdm;

  scalartype* fnc_values;
};

template <typename scalartype, class domain>
function<scalartype, domain>::function(const std::string& fnc_name)
    : name_(fnc_name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_elements(dmn.get_size()),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(dmn.get_leaf_domain_steps()) {
  fnc_values = new scalartype[Nb_elements];

  for (int linind = 0; linind < Nb_elements; linind++)
    setToZero(fnc_values[linind]);
}

template <typename scalartype, class domain>
function<scalartype, domain>::function() : function("no name") {}

template <typename scalartype, class domain>
function<scalartype, domain>::function(const function<scalartype, domain>& other_one)
    : name_(other_one.name_),
      function_type(__PRETTY_FUNCTION__),
      dmn(other_one.dmn),
      Nb_elements(other_one.Nb_elements),
      Nb_sbdms(other_one.Nb_sbdms),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(dmn.get_leaf_domain_steps()) {
  fnc_values = new scalartype[Nb_elements];
  std::copy_n(other_one.fnc_values, Nb_elements, fnc_values);
}

// Move contructor.
template <typename scalartype, class domain>
function<scalartype, domain>::function(function<scalartype, domain>&& other_one)
    : name_(std::move(other_one.name_)),
      function_type(__PRETTY_FUNCTION__),
      dmn(std::move(other_one.dmn)),
      Nb_elements(other_one.Nb_elements),
      Nb_sbdms(other_one.Nb_sbdms),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(dmn.get_leaf_domain_steps()) {
  fnc_values = other_one.fnc_values;
  other_one.fnc_values = nullptr;
  other_one.Nb_elements = 0;
}

template <typename scalartype, class domain>
function<scalartype, domain>::~function() {
  delete[] fnc_values;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::reset() {
  dmn.reset();

  Nb_sbdms = dmn.get_leaf_domain_sizes().size();
  Nb_elements = dmn.get_size();

  delete[] fnc_values;
  fnc_values = new scalartype[Nb_elements];

  for (int linind = 0; linind < Nb_elements; linind++)
    setToZero(fnc_values[linind]);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::linind_2_subind(int linind, int*& subind) const {
  int tmp = linind;
  for (int i = 0; i < int(size_sbdm.size()); i++) {
    subind[i] = tmp % size_sbdm[i];
    tmp = (tmp - subind[i]) / size_sbdm[i];
  }
}

template <typename scalartype, class domain>
void function<scalartype, domain>::linind_2_subind(int linind, std::vector<int>& subind) const {
  assert(int(subind.size()) == signature());

  int tmp = linind;
  for (int i = 0; i < int(size_sbdm.size()); i++) {
    subind[i] = tmp % size_sbdm[i];
    tmp = (tmp - subind[i]) / size_sbdm[i];
  }
}

template <typename scalartype, class domain>
void function<scalartype, domain>::subind_2_linind(const int* subind, int& linind) const {
  linind = 0;
  for (int i = 0; i < int(step_sbdm.size()); i++)
    linind += subind[i] * step_sbdm[i];
}

template <typename scalartype, class domain>
scalartype& function<scalartype, domain>::operator()(const int* subind) {
  int linind;
  subind_2_linind(subind, linind);

  assert(linind >= 0 && linind < Nb_elements);
  return fnc_values[linind];
}

template <typename scalartype, class domain>
const scalartype& function<scalartype, domain>::operator()(const int* subind) const {
  int linind;
  subind_2_linind(subind, linind);

  assert(linind >= 0 && linind < Nb_elements);
  return fnc_values[linind];
}

template <typename scalartype, class domain>
function<scalartype, domain>& function<scalartype, domain>::operator=(
    const function<scalartype, domain>& f_other) {
  const domain& dmn_other = f_other.get_domain();

  if (dmn.get_size() !=
      dmn_other.get_size())  // Domains were not initialized when function was created.
  {
    dmn = dmn_other;

    Nb_sbdms = f_other.Nb_sbdms;
    assert(Nb_sbdms == dmn.get_leaf_domain_sizes().size());
    Nb_elements = f_other.Nb_elements;
    assert(Nb_elements == dmn.get_size());

    delete[] fnc_values;
    fnc_values = new scalartype[Nb_elements];
  }

  std::copy_n(f_other.values(), Nb_elements, fnc_values);
  return *this;
}

// Move assignment.
template <typename scalartype, class domain>
function<scalartype, domain>& function<scalartype, domain>::operator=(
    function<scalartype, domain>&& f_other) {
  if (this != &f_other) {
    domain& dmn_other = f_other.dmn;
    if (dmn.get_size() !=
        dmn_other.get_size())  // Domains were not initialized when function was created.
    {
      std::swap(dmn, dmn_other);

      Nb_sbdms = f_other.Nb_sbdms;
      assert(Nb_sbdms == dmn.get_leaf_domain_sizes().size());
      Nb_elements = f_other.Nb_elements;
      assert(Nb_elements == dmn.get_size());
    }

    delete[] fnc_values;
    fnc_values = f_other.fnc_values;
    f_other.fnc_values = nullptr;
    f_other.Nb_elements = 0;
  }
  return *this;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator+=(const function<scalartype, domain>& f_other) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] += f_other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator-=(const function<scalartype, domain>& f_other) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] -= f_other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator*=(const function<scalartype, domain>& f_other) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] *= f_other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator/=(const function<scalartype, domain>& f_other) {
  for (int linind = 0; linind < Nb_elements; linind++) {
    assert(std::abs(f_other(linind)) > 1.e-16);
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
bool function<scalartype, domain>::operator==(const function<scalartype, domain>& f_other) const{
  assert(size() == f_other.size());
  for(int i = 0; i < Nb_elements; i++)
    if(f_other(i) != fnc_values[i])
      return false;

  return true;
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::slice(int sbdm_index, int* subind, new_scalartype* fnc_vals) const {
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
                                         new_scalartype* fnc_vals) const {
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
void function<scalartype, domain>::distribute(int sbdm_index, int* subind,
                                              const new_scalartype* fnc_vals) {
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
                                              const new_scalartype* fnc_vals) {
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
void function<scalartype, domain>::print_fingerprint(std::ostream& stream) const {
  stream << "****************************************\n";
  stream << "function: " << name_ << "\n";
  stream << "****************************************\n";

  stream << "# subdomains: " << Nb_sbdms << "\n";
  dca::util::print_type<domain>::print(stream);
  stream << "\n";

  stream << "size of subdomains:";
  for (int i = 0; i < Nb_sbdms; ++i)
    stream << "  " << size_sbdm[i];
  stream << "\n";

  stream << "# elements: " << Nb_elements << "\n";
  stream << "memory: " << Nb_elements * sizeof(scalartype) / (1024. * 1024.) << " MiB\n";
  stream << "****************************************\n" << std::endl;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::print_elements(std::ostream& stream) const {
  stream << "****************************************\n";
  stream << "function: " << name_ << "\n";
  stream << "****************************************\n";

  std::vector<int> subind(Nb_sbdms);
  for (int lindex = 0; lindex < Nb_elements; ++lindex) {
    linind_2_subind(lindex, subind);
    for (int index : subind)
      stream << index << "\t";
    stream << " \t" << fnc_values[lindex] << "\n";
  }

  stream << "****************************************\n" << std::endl;
}

template <typename scalartype, class domain>
template <typename concurrency_t>
int function<scalartype, domain>::get_buffer_size(const concurrency_t& concurrency) const {
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
