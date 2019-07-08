// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//
// This class connects the function values to the domains.
//
// TODO: Remame template parameter scalartype --> ElementType.

#ifndef DCA_FUNCTION_FUNCTION_HPP
#define DCA_FUNCTION_FUNCTION_HPP

#include <algorithm>  // std::copy_n
#include <cassert>
#include <cmath>    // std::abs
#include <complex>  // std::abs(std::complex)
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>  // std::is_integral
#include <utility>      // std::move, std::swap
#include <vector>

#include "dca/function/scalar_cast.hpp"
#include "dca/function/set_to_zero.hpp"
#include "dca/util/pack_operations.hpp"
#include "dca/util/type_utils.hpp"

namespace dca {
namespace func {
// dca::func::

template <typename scalartype, class domain>
class function {
  static const std::string default_name_;

public:
  typedef scalartype this_scalar_type;
  typedef domain this_domain_type;

  // Default constructor
  // Constructs the function with the name name.
  // Postcondition: All elements are set to zero.
  function(const std::string& name = default_name_);

  // Copy constructor
  // Constructs the function with the name name and a copy of the elements of other.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  function(const function<scalartype, domain>& other, const std::string& name = default_name_);

  // Move constructor
  // Constructs the function with the name name and the elements of other using move semantics.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  // Postcondition: The other function is in a non-specified state.
  function(function<scalartype, domain>&& other, const std::string& name = default_name_);

  // Copy assignment operator
  // Replaces the function's elements with a copy of the elements of other.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  // Postcondition: The function's name is unchanged.
  function<scalartype, domain>& operator=(const function<scalartype, domain>& other);

  // Move assignment operator
  // Replaces the function's elements with those of other using move semantics.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  // Postconditions: The function's name is unchanged.
  //                 The other function is in a non-specified state.
  function<scalartype, domain>& operator=(function<scalartype, domain>&& other);

  ~function();

  // Resets the function by resetting the domain object and reallocating the memory for the function
  // elements.
  // Postcondition: All elements are set to zero.
  void reset();

  const domain& get_domain() const {
    return dmn;
  }
  const std::string& get_name() const {
    return name_;
  }
  // TODO: Remove this method and use constructor parameter instead.
  void set_name(const std::string& name) {
    name_ = name;
  }
  int signature() const {
    return Nb_sbdms;
  }
  std::size_t size() const {
    return Nb_elements;
  }
  // Returns the size of the leaf domain with the given index.
  // Does not return function values!
  int operator[](const int index) const {
    return size_sbdm[index];
  }

  // Returns a pointer to the function's elements.
  scalartype* values() {
    return fnc_values;
  }
  const scalartype* values() const {
    return fnc_values;
  }
  scalartype* data() {
    return fnc_values;
  }
  const scalartype* data() const {
    return fnc_values;
  }

  //
  // Methods for index conversion
  //
  // Converts the linear index to the corresponding subindices of the leaf domains.
  // Pointer version
  // Precondition: The size of the array pointed to by subind must be equal to the number of leaf
  //               domains (Nb_sbdms).
  // TODO: Replace pointer version with std::array to be able to check subind's size.
  void linind_2_subind(int linind, int* subind) const;
  // std::vector version
  void linind_2_subind(int linind, std::vector<int>& subind) const;

  // Computes the linear index for the given subindices of the leaf domains.
  // Precondition: subind stores the the subindices of all LEAF domains.
  // TODO: Use std::array or std::vector to be able to check the size of subind.
  void subind_2_linind(const int* subind, int& linind) const;

  // Computes and returns the linear index for the given subindices of the branch or leaf domains,
  // depending on the size of subindices.
  // Enable only if all arguments are integral to prevent subind_to_linind(int*, int) to resolve to
  // subind_to_linind(int...) rather than subind_to_linind(const int* const, int).
  template <typename... Ts>
  std::enable_if_t<util::if_all<std::is_integral<Ts>::value...>::value, int> subind_2_linind(
      const Ts... subindices) const {
    // We need to cast all subindices to the same type for dmn_variadic.
    return dmn(static_cast<int>(subindices)...);
  }

  // TODO: Remove this method.
  template <typename T>
  int subind_2_linind(const T ind) const {
    static_assert(std::is_integral<T>::value, "Index ind must be an integer.");
    assert(ind >= 0 && ind < Nb_elements);
    return ind;
  }

  //
  // operator()
  //
  // TODO: Remove these two methods and use the variadic domains versions instead.
  scalartype& operator()(const int* subind);
  const scalartype& operator()(const int* subind) const;

  template <typename T>
  scalartype& operator()(const T linind) {
    static_assert(std::is_integral<T>::value, "Index linind must be an integer.");
    assert(linind >= 0 && linind < Nb_elements);
    return fnc_values[linind];
  }
  template <typename T>
  const scalartype& operator()(const T linind) const {
    static_assert(std::is_integral<T>::value, "Index linind must be an integer.");
    assert(linind >= 0 && linind < Nb_elements);
    return fnc_values[linind];
  }

  template <typename... Ts>
  scalartype& operator()(const Ts... subindices) {
    // We need to cast all indices to the same type for dmn_variadic.
    return fnc_values[dmn(static_cast<int>(subindices)...)];
  }
  template <typename... Ts>
  const scalartype& operator()(const Ts... subindices) const {
    return fnc_values[dmn(static_cast<int>(subindices)...)];
  }

  void operator+=(const function<scalartype, domain>& other);
  void operator-=(const function<scalartype, domain>& other);
  void operator*=(const function<scalartype, domain>& other);
  void operator/=(const function<scalartype, domain>& other);


  void operator=(scalartype c);
  void operator+=(scalartype c);
  void operator-=(scalartype c);
  void operator*=(scalartype c);
  void operator/=(scalartype c);

  // Equal-comparison opertor
  // Returns true if the function's elements (fnc_values) are equal to other's elements, false
  // otherwise.
  // TODO: Make the equal-comparison operator a non-member function.
  bool operator==(const function<scalartype, domain>& other) const;

  template <typename new_scalartype>
  void slice(int sbdm_index, int* subind, new_scalartype* fnc_vals) const;
  template <typename new_scalartype>
  void slice(int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals) const;
  template <typename new_scalartype>
  void distribute(int sbdm_index, int* subind, const new_scalartype* fnc_vals);
  template <typename new_scalartype>
  void distribute(int sbdm_index_1, int sbdm_index_2, int* subind, const new_scalartype* fnc_vals);

  //
  // Methods for printing
  //
  // Prints the function's metadata.
  void print_fingerprint(std::ostream& stream = std::cout) const;
  // Prints the function's elements.
  void print_elements(std::ostream& stream = std::cout) const;

  //
  // Methods for message passing concurrency
  //
  template <typename concurrency_t>
  int get_buffer_size(const concurrency_t& concurrency) const;
  template <class concurrency_t>
  void pack(const concurrency_t& concurrency, char* buffer, int buffer_size, int& position) const;
  // TODO: Make parameter buffer const correct (const char* const buffer).
  template <class concurrency_t>
  void unpack(const concurrency_t& concurrency, char* buffer, int buffer_size, int& position);

private:
  std::string name_;
  std::string function_type;

  domain dmn;  // TODO: Remove domain object?

  std::size_t Nb_elements;

  // The subdomains (sbdmn) represent the leaf domains, not the branch domains.
  int Nb_sbdms;
  const std::vector<std::size_t>& size_sbdm;  // TODO: Remove?
  const std::vector<std::size_t>& step_sbdm;  // TODO: Remove?

  scalartype* fnc_values;
};

template <typename scalartype, class domain>
const std::string function<scalartype, domain>::default_name_ = "no-name";

template <typename scalartype, class domain>
function<scalartype, domain>::function(const std::string& name)
    : name_(name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_elements(dmn.get_size()),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(dmn.get_leaf_domain_steps()),
      fnc_values(nullptr) {
  fnc_values = new scalartype[Nb_elements];
  for (int linind = 0; linind < Nb_elements; ++linind)
    setToZero(fnc_values[linind]);
}

template <typename scalartype, class domain>
function<scalartype, domain>::function(const function<scalartype, domain>& other,
                                       const std::string& name)
    : name_(name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_elements(dmn.get_size()),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(dmn.get_leaf_domain_steps()),
      fnc_values(nullptr) {
  if (dmn.get_size() != other.dmn.get_size())
    // The other function has not been resetted after the domain was initialized.
    throw std::logic_error("Copy construction from a not yet resetted function.");

  fnc_values = new scalartype[Nb_elements];
  std::copy_n(other.fnc_values, Nb_elements, fnc_values);
}

template <typename scalartype, class domain>
function<scalartype, domain>::function(function<scalartype, domain>&& other, const std::string& name)
    : name_(name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_elements(dmn.get_size()),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(dmn.get_leaf_domain_steps()),
      fnc_values(nullptr) {
  if (dmn.get_size() != other.dmn.get_size())
    // The other function has not been resetted after the domain was initialized.
    throw std::logic_error("Move construction from a not yet resetted function.");

  fnc_values = other.fnc_values;
  other.Nb_elements = 0;
  other.fnc_values = nullptr;
}

template <typename scalartype, class domain>
function<scalartype, domain>& function<scalartype, domain>::operator=(
    const function<scalartype, domain>& other) {
  if (this != &other) {
    if (dmn.get_size() != other.dmn.get_size()) {
      // Domain had not been initialized when the functions were created.
      // Reset this function and check again.
      reset();

      if (dmn.get_size() != other.dmn.get_size())
        // The other function has not been resetted after the domain was initialized.
        throw std::logic_error("Copy assignment from a not yet resetted function.");
    }

    std::copy_n(other.values(), Nb_elements, fnc_values);
  }

  return *this;
}

template <typename scalartype, class domain>
function<scalartype, domain>& function<scalartype, domain>::operator=(
    function<scalartype, domain>&& other) {
  if (this != &other) {
    if (dmn.get_size() != other.dmn.get_size()) {
      // Domain had not been initialized when the functions were created.
      // Reset this function and check again.
      reset();

      if (dmn.get_size() != other.dmn.get_size())
        // The other function has not been resetted after the domain was initialized.
        throw std::logic_error("Move assignment from a not yet resetted function.");
    }

    delete[] fnc_values;
    fnc_values = other.fnc_values;

    other.Nb_elements = 0;
    other.fnc_values = nullptr;
  }

  return *this;
}

template <typename scalartype, class domain>
function<scalartype, domain>::~function() {
  delete[] fnc_values;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::reset() {
  dmn.reset();

  Nb_elements = dmn.get_size();
  Nb_sbdms = dmn.get_leaf_domain_sizes().size();

  delete[] fnc_values;
  fnc_values = new scalartype[Nb_elements];

  for (int linind = 0; linind < Nb_elements; ++linind)
    setToZero(fnc_values[linind]);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::linind_2_subind(int linind, int* subind) const {
  for (int i = 0; i < int(size_sbdm.size()); ++i) {
    subind[i] = linind % size_sbdm[i];
    linind = (linind - subind[i]) / size_sbdm[i];
  }
}

// TODO: Resize vector if necessary.
template <typename scalartype, class domain>
void function<scalartype, domain>::linind_2_subind(int linind, std::vector<int>& subind) const {
  assert(int(subind.size()) == Nb_sbdms);

  for (int i = 0; i < int(size_sbdm.size()); ++i) {
    subind[i] = linind % size_sbdm[i];
    linind = (linind - subind[i]) / size_sbdm[i];
  }
}

template <typename scalartype, class domain>
void function<scalartype, domain>::subind_2_linind(const int* const subind, int& linind) const {
  linind = 0;
  for (int i = 0; i < int(step_sbdm.size()); ++i)
    linind += subind[i] * step_sbdm[i];
}

template <typename scalartype, class domain>
scalartype& function<scalartype, domain>::operator()(const int* const subind) {
  int linind;
  subind_2_linind(subind, linind);

  assert(linind >= 0 && linind < Nb_elements);
  return fnc_values[linind];
}

template <typename scalartype, class domain>
const scalartype& function<scalartype, domain>::operator()(const int* const subind) const {
  int linind;
  subind_2_linind(subind, linind);

  assert(linind >= 0 && linind < Nb_elements);
  return fnc_values[linind];
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator+=(const function<scalartype, domain>& other) {
  for (int linind = 0; linind < Nb_elements; ++linind)
    fnc_values[linind] += other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator-=(const function<scalartype, domain>& other) {
  for (int linind = 0; linind < Nb_elements; ++linind)
    fnc_values[linind] -= other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator*=(const function<scalartype, domain>& other) {
  for (int linind = 0; linind < Nb_elements; ++linind)
    fnc_values[linind] *= other(linind);
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator/=(const function<scalartype, domain>& other) {
  for (int linind = 0; linind < Nb_elements; ++linind) {
    assert(std::abs(other(linind)) > 1.e-16);
    fnc_values[linind] /= other(linind);
  }
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator=(const scalartype c) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] = c;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator+=(const scalartype c) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] += c;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator-=(const scalartype c) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] -= c;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator*=(const scalartype c) {
  for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] *= c;
}

template <typename scalartype, class domain>
void function<scalartype, domain>::operator/=(const scalartype c) {
    for (int linind = 0; linind < Nb_elements; linind++)
    fnc_values[linind] /= c;
}

template <typename scalartype, class domain>
bool function<scalartype, domain>::operator==(const function<scalartype, domain>& other) const {
  if (size() != other.size())
    // One of the function has not been resetted after the domain was initialized.
    throw std::logic_error("Comparing functions of different sizes.");

  for (int i = 0; i < Nb_elements; ++i)
    if (other(i) != fnc_values[i])
      return false;

  return true;
}

template <typename scalartype, class domain>
template <typename new_scalartype>
void function<scalartype, domain>::slice(const int sbdm_index, int* subind,
                                         new_scalartype* fnc_vals) const {
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
void function<scalartype, domain>::slice(const int sbdm_index_1, const int sbdm_index_2,
                                         int* subind, new_scalartype* fnc_vals) const {
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
void function<scalartype, domain>::distribute(const int sbdm_index, int* subind,
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
void function<scalartype, domain>::distribute(const int sbdm_index_1, const int sbdm_index_2,
                                              int* subind, const new_scalartype* fnc_vals) {
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
  util::print_type<domain>::print(stream);
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
void function<scalartype, domain>::pack(const concurrency_t& concurrency, char* buffer,
                                        const int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, *this);
}

template <typename scalartype, class domain>
template <class concurrency_t>
void function<scalartype, domain>::unpack(const concurrency_t& concurrency, char* buffer,
                                          const int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, *this);
}

}  // func
}  // dca

#endif  // DCA_FUNCTION_FUNCTION_HPP
