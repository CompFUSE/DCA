// Copyright (C) 2022 ETH Zurich
// Copyright (C) 2022 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Staar (taa@zurich.ibm.com)
//         Urs R. Haehner (haehneru@itp.phys.ethz.ch)
//         Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//         Peter Doak (doakpw@ornl.gov)
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
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <type_traits>  // std::is_integral
#include <utility>      // std::move, std::swap
#include <vector>
#include "dca/distribution/dist_types.hpp"
#include "dca/function/scalar_cast.hpp"
#include "dca/function/set_to_zero.hpp"
#include "dca/util/pack_operations.hpp"
#include "dca/util/integer_division.hpp"
#include "dca/util/type_utils.hpp"
#include "dca/util/to_string.hpp"

namespace dca {
namespace func {
// dca::func::

/** The tensor class used through DCA++.
 *  This is a pretty complex construct but probably helps keeping track of the large number of
 * different tensors over domains in the code. The memory layout is close packed so it is generally
 * necessary to copy slices of functions into matrices or vectors for efficient calculation.
 *
 *  First domain is fastest, indexes are in order of domain.
 *  example:
 *  func<double, dmn_variadic<dmn_0<4, double>, dmn_0<8, double>> a_func
 *  0  1  2  3
 *  4  5  6  7
 *  8  9  10 11
 *  12 13 14 15
 *  16 17 18 19
 *  20 21 22 23
 *  24 25 26 27
 *  28 29 30 21
 *
 *  the following is true
 *  a_func(0,3) == 12
 *  a_func(3,0) == 3
 *
 *  i.e. row major layout with column first indexing.
 */
template <typename scalartype, class domain, DistType DT = DistType::NONE>
class function;

template <typename scalartype, class domain, DistType DT>
class function {
  static const std::string default_name_;

public:
  static constexpr DistType dist = DT;
  typedef scalartype this_scalar_type;
  typedef domain this_domain_type;

  // Default constructor
  // Constructs the function with the name name.
  // Postcondition: All elements are set to zero.
  function(const std::string& name = default_name_);

  // Distributed function. Access with multi-index operator() is not safe.
  template <class Concurrency>
  function(const std::string& name, const Concurrency& concurrency);

  // Copy constructor
  // Constructs the function with the a copy of elements and name of other.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  function(const function<scalartype, domain, DT>& other);
  // Same as above, but with name change from name argument.
  function(const function<scalartype, domain, DT>& other, const std::string& name);

  // Move constructor
  // Constructs the function with elements and name of other using move semantics.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  // Postcondition: The other function is in a non-specified state.
  function(function<scalartype, domain, DT>&& other);
  function(function<scalartype, domain, DT>&& other, const std::string& name)
      : function(std::move(other)) {
    name_ = name;
  }

  // Initializer list constructor
  function(std::initializer_list<scalartype> init_list, const std::string& name = default_name_);

  // Copy assignment operator
  // Replaces the function's elements with a copy of the elements of other.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  // Postcondition: The function's name is unchanged.
  function<scalartype, domain, DT>& operator=(const function<scalartype, domain, DT>& other);
  template <typename Scalar2>
  function<scalartype, domain, DT>& operator=(const function<Scalar2, domain, DT>& other);

  // Move assignment operator
  // Replaces the function's elements with those of other using move semantics.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  // Postconditions: The function's name is unchanged.
  //                 The other function is in a non-specified state.
  function<scalartype, domain, DT>& operator=(function<scalartype, domain, DT>&& other);

  // Resets the function by resetting the domain object and reallocating the memory for the function
  // elements.
  // \todo These odd semantics need serious justification.
  // Postcondition: All elements are set to zero.
  template <class Concurrency>
  void reset(const Concurrency& conc);

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

  std::size_t get_start() const {
    return start_;
  }
  /** end in sense of last index not 1 past.
   */
  std::size_t get_end() const {
    assert(end_ > 0);
    return end_ - 1;
  }

  std::vector<size_t> get_start_subindex() const {
    return linind_2_subind(start_);
  }
  /** end in sense of last subindex
   */
  std::vector<size_t> get_end_subindex() const {
    assert(end_ > 0);
    return linind_2_subind(end_ - 1);
  }

  int signature() const {
    return Nb_sbdms;
  }
  std::size_t size() const {
    return fnc_values_.size();
  }

  // TODO: remove as it breaks class' invariant.
  void resize(std::size_t nb_elements_new) {
    fnc_values_.resize(nb_elements_new);
  }

  // void local_resize(std::size_t nb_elements) {}

  // Returns the size of the leaf domain with the given index.
  // Does not return function values!
  // Broken off except on rank 0.
  int operator[](const int index) const {
    return dmn.get_leaf_domain_sizes()[index];
  }

  const auto& getDomainSizes() const noexcept {
    return dmn.get_leaf_domain_sizes();
  }
  const std::vector<scalartype>& getValues() const noexcept {
    return fnc_values_;
  }

  // Begin and end methods for compatibility with range for loop.
  auto begin() {
    return fnc_values_.begin();
  }
  auto end() {
    return fnc_values_.end();
  }
  auto begin() const {
    return fnc_values_.begin();
  }
  auto end() const {
    return fnc_values_.end();
  }

  // Returns a pointer to the function's elements.
  scalartype* values() {
    return fnc_values_.data();
  }
  const scalartype* values() const {
    return fnc_values_.data();
  }
  scalartype* data() {
    return fnc_values_.data();
  }
  const scalartype* data() const {
    return fnc_values_.data();
  }

  //
  // Methods for index conversion
  //
  // Converts the linear index to the corresponding subindices of the leaf domains.
  // Pointer version
  // Precondition: The size of the array pointed to by subind must be equal to the number of leaf
  //               domains (Nb_sbdms).
  // \todo Replace pointer version with std::array to be able to check subind's size.
  // \todo validate or not usage of these for distributed (across MPI) functions, I strongly suspect they are
  //       not ok./
  void linind_2_subind(int linind, int* subind) const;
  // std::vector version
  void linind_2_subind(int linind, std::vector<int>& subind) const;

  template <std::size_t N>
  void linind_2_subind(int linind, std::array<int, N>& subind) const;

  // modern RVO version
  std::vector<size_t> linind_2_subind(int linind) const;

  // Computes the linear index for the given subindices of the leaf domains.
  // Precondition: subind stores the the subindices of all LEAF domains.
  // TODO: Use std::array or std::vector to be able to check the size of subind.
  void subind_2_linind(const int* subind, int& linind) const;

  // using standard vector and avoiding returning argument
  size_t subind_2_linind(const std::vector<int>& subind) const;

  // using standard vector and avoiding returning argument
  size_t branch_subind_2_linind(const std::vector<int>& subind) const;

  // Computes and returns the linear index for the given subindices of the branch or leaf domains,
  // depending on the size of subindices.
  // Enable only if all arguments are integral to prevent subind_to_linind(int*, int) to resolve to
  // subind_to_linind(int...) rather than subind_to_linind(const int* const, int).
  template <typename... Ts>
  std::enable_if_t<util::ifAll(std::is_integral_v<Ts>...), int> subind_2_linind(
      const Ts... subindices) const {
    // We need to cast all subindices to the same type for dmn_variadic.
    return dmn(static_cast<int>(subindices)...);
  }

  // TODO: Remove this method.
  template <typename T>
  int subind_2_linind(const T ind) const {
    static_assert(std::is_integral<T>::value, "Index ind must be an integer.");
    assert(ind >= 0 && ind < size());
    return ind;
  }

  //
  // operator()
  //
  // TODO: Remove these two methods and use the variadic domains versions instead.
  scalartype& operator()(const int* subind);
  const scalartype& operator()(const int* subind) const;

  template <std::size_t N>
  scalartype& operator()(const std::array<int, N>& subind) {
#ifndef NDEBUG
    auto linind = dmn.index_by_array(subind);
    assert(linind >= 0 && linind < size());
    return fnc_values_[linind];
#else
    return fnc_values_[dmn.index_by_array(subind)];
#endif
  }

  const scalartype& operator()(const std::vector<int>& subind) const;

  template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, void>::type>
  scalartype& operator()(const T linind) {
    static_assert(std::is_integral<T>::value, "Index linind must be an integer.");
    assert(linind >= 0 && linind < size());
    return fnc_values_[linind];
  }
  template <typename T, typename = typename std::enable_if<std::is_integral<T>::value, void>::type>
  const scalartype& operator()(const T linind) const {
    static_assert(std::is_integral<T>::value, "Index linind must be an integer.");
    assert(linind >= 0 && linind < size());
    return fnc_values_[linind];
  }

  template <typename T, typename... Ts,
            typename = typename std::enable_if<std::is_integral<T>::value, void>::type>
  scalartype& operator()(const T t, const Ts... subindices) {
    // We need to cast all indices to the same type for dmn_variadic.
    return fnc_values_[dmn(static_cast<int>(t), static_cast<int>(subindices)...)];
  }
  template <typename T, typename... Ts,
            typename = typename std::enable_if<std::is_integral<T>::value, void>::type>
  const scalartype& operator()(const T t, const Ts... subindices) const {
    return fnc_values_[dmn(static_cast<int>(t), static_cast<int>(subindices)...)];
  }

  void operator+=(const function<scalartype, domain, DT>& other);
  void operator-=(const function<scalartype, domain, DT>& other);
  void operator*=(const function<scalartype, domain, DT>& other);
  void operator/=(const function<scalartype, domain, DT>& other);

  void operator=(scalartype c);
  void operator+=(scalartype c);
  void operator-=(scalartype c);
  void operator*=(scalartype c);
  void operator/=(scalartype c);

  // Equal-comparison opertor
  // Returns true if the function's elements (fnc_values_) are equal to other's elements, false
  // otherwise.
  // TODO: Make the equal-comparison operator a non-member function.
  bool operator==(const function<scalartype, domain, DT>& other) const;

  /** slice across 1 leaf domain
   *  deprecated, really unsafe pointer use
   */
  template <typename new_scalartype>
  void slice(int sbdm_index, int* subind, new_scalartype* fnc_vals) const;
  /** slice across leaf domains or branch domains based on size of subind.
   *  \param[in] subind    subindex of the slice, index in the slice domain is ignored
   */
  template <typename new_scalartype>
  void slice(const int sbdm_index, std::vector<int> subind, new_scalartype* fnc_vals) const;

  template <typename new_scalartype>
  void slice(int sbdm_index_1, int sbdm_index_2, std::vector<int> subind,
             new_scalartype* fnc_vals) const;
  template <typename new_scalartype>
  void slice(int sbdm_index_1, int sbdm_index_2, int* subind, new_scalartype* fnc_vals) const;

  /** write a slice to a larger function.
   *  The thread safety of this is questionable especially if the sbdm_index > 0
   *  i.e. we are not writing a fastest slice.
   */
  template <typename new_scalartype>
  void distribute(int sbdm_index, std::vector<int> subind, const new_scalartype* fnc_vals);
  template <typename new_scalartype>
  void distribute(int sbdm_index, int* subind, const new_scalartype* fnc_vals);
  template <typename new_scalartype>
  void distribute(int sbdm_index_1, int sbdm_index_2, std::vector<int> subind,
                  const new_scalartype* fnc_vals);
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

  // Gather a function that was initialized as distributed.
  // Precondition: concurrency must be the same object used during construction.
  template <class Concurrency>
  inline function gather(const Concurrency& concurrency) const;

private:
  // For DistType::BLOCKED and DistType::LINEAR calculates the local_function_size and sets start_ and end_.
  template <class Concurrency>
  std::size_t calcDistribution(const Concurrency& concurrency);

  std::string name_;
  std::string function_type;

  domain dmn;  // TODO: Remove domain object?

  // The subdomains (sbdmn) represent the leaf domains, not the branch domains.
  int Nb_sbdms;

  std::vector<scalartype> fnc_values_;

  // These are the linear start and end indexes with respect to the complete function.
  std::size_t start_;
  std::size_t end_;
};

template <typename scalartype, class domain, DistType DT>
const std::string function<scalartype, domain, DT>::default_name_ = "no-name";

/** default constructor
 */
template <typename scalartype, class domain, DistType DT>
function<scalartype, domain, DT>::function(const std::string& name)
    : name_(name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()) {
  if constexpr (dist == DistType::BLOCKED || dist == DistType::LINEAR) {
    throw std::runtime_error(
        "function named constructor without concurrency reference may only be called for "
        "DistType::NONE");
  }
  start_ = 0;
  end_ = dmn.get_size();
  // will zero real or complex values
  fnc_values_.resize(dmn.get_size(), {});
}

/** copy constructor
 */
template <typename scalartype, class domain, DistType DT>
function<scalartype, domain, DT>::function(const function<scalartype, domain, DT>& other)
    : name_(other.name_),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      fnc_values_(other.fnc_values_) {
  start_ = other.start_;
  end_ = other.end_;
}

/** name change copy constructor
 */
template <typename scalartype, class domain, DistType DT>
function<scalartype, domain, DT>::function(const function<scalartype, domain, DT>& other,
                                           const std::string& name)
    : name_(name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      fnc_values_(other.fnc_values_) {
  start_ = other.start_;
  end_ = other.end_;
}

/** move constructor */
template <typename scalartype, class domain, DistType DT>
function<scalartype, domain, DT>::function(function<scalartype, domain, DT>&& other)
    : name_(std::move(other.name_)),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      fnc_values_(std::move(other.fnc_values_)) {
  if (dmn.get_size() != other.dmn.get_size())
    // The other function has not been reset after the domain was initialized.
    throw std::logic_error("Move construction from a not yet resetted function.");
  start_ = other.start_;
  end_ = other.end_;
}

/** initializer list constructor
 *  only needed for testing.
 */
template <typename scalartype, class domain, DistType DT>
function<scalartype, domain, DT>::function(std::initializer_list<scalartype> init_list,
                                           const std::string& name)
    : name_(name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()) {
  start_ = 0;
  end_ = dmn.get_size();
  fnc_values_.resize(dmn.get_size(), {});
  std::copy_n(init_list.begin(), init_list.size(), fnc_values_.begin());
}

template <typename scalartype, class domain, DistType DT>
template <class Concurrency>
std::size_t function<scalartype, domain, DT>::calcDistribution(const Concurrency& concurrency) {
  if constexpr (dist == DistType::NONE) {
    throw std::logic_error("calcDistritribution should not be called for DistType::NONE function.");
  }
  else if constexpr (dist == DistType::BLOCKED) {
    auto error_bad_block = [](size_t conc_size, domain& dmn) {
      std::ostringstream error_message;
      double block_size = (dmn.get_size() * sizeof(scalartype)) / 1024 / 1024 / 1024;  // gigabytes
      error_message << "Blocked concurrency is not possible. Concurrency size: " << conc_size
                    << " is not blockwise divisor of function with dimensions:\n"
                    << vectorToString(dmn.get_leaf_domain_sizes()) << '\n'
                    << "Total Blocked Function Size: " << block_size << "GB\n";
      throw std::runtime_error(error_message.str());
    };
    const std::size_t my_concurrency_id = concurrency.id();
    const std::size_t my_concurrency_size = concurrency.number_of_processors();
    std::size_t local_function_size = dca::util::ceilDiv(dmn.get_size(), my_concurrency_size);
    // This is a necessary but not sufficient proof of "regular blocking"
    if (local_function_size * my_concurrency_size != dmn.get_size()) {
      error_bad_block(my_concurrency_size, dmn);
    }
    start_ = local_function_size * my_concurrency_id;
    end_ = start_ + local_function_size;
    bool regular_local_function_size = false;
    size_t remaining_ranks = my_concurrency_size;
    for (int idim = (dmn.get_Nb_leaf_domains() - 1); idim >= 0; --idim) {
      if (remaining_ranks < dmn.get_subdomain_size(idim)) {
        if (dmn.get_subdomain_size(idim) % remaining_ranks != 0) {
          break;  // i.e no regular bricked blocking
        }
        else {
          regular_local_function_size = true;
          break;
        }
      }
      else {
        if (remaining_ranks % dmn.get_subdomain_size(idim) != 0) {
          break;  // i.e no regular bricked blocking
        }
        else {
          remaining_ranks /= dmn.get_subdomain_size(idim);
        }
      }
    }
    if (!regular_local_function_size) {
      error_bad_block(my_concurrency_size, dmn);
    }
    return local_function_size;
  }
  else if constexpr (dist == DistType::LINEAR) {
    const std::size_t my_concurrency_id = concurrency.id();
    const std::size_t my_concurrency_size = concurrency.number_of_processors();
    size_t local_function_size = dca::util::ceilDiv(dmn.get_size(), my_concurrency_size);
    size_t residue = dmn.get_size() % my_concurrency_size;
    start_ = local_function_size * my_concurrency_id;
    if (residue != 0 && my_concurrency_id > residue - 1) {
      start_ -= my_concurrency_id - residue;
      --local_function_size;
    }
    end_ = start_ + local_function_size;
    return local_function_size;
  }
}

/** distributed function constructor
 */
template <typename scalartype, class domain, DistType DT>
template <class Concurrency>
function<scalartype, domain, DT>::function(const std::string& name, const Concurrency& concurrency)
    : name_(name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()) {
  if constexpr (dist == DistType::NONE) {
    const std::size_t nb_elements = dmn.get_size();
    try {
      fnc_values_.resize(nb_elements);
    }
    catch (const std::exception& exc) {
      std::cout << "exception caught on resize of blocked function: " << exc.what() << '\n';
      throw(exc);
    }
    for (int linind = 0; linind < nb_elements; ++linind)
      setToZero(fnc_values_[linind]);
    start_ = 0;
    end_ = dmn.get_size();
  }
  else if constexpr (dist == DistType::LINEAR) {
    std::size_t local_function_size = calcDistribution(concurrency);
    try {
      fnc_values_.resize(local_function_size);
    }
    catch (const std::exception& exc) {
      std::cout << "exception caught on resize of linear distributed function: " << exc.what()
                << '\n';
      throw(exc);
    }
    for (int linind = 0; linind < local_function_size; ++linind)
      setToZero(fnc_values_[linind]);
  }
  else if constexpr (dist == DistType::BLOCKED) {
    std::size_t local_function_size = calcDistribution(concurrency);
    // Ok this can be a blocked function so we finally resize i.e. allocate.
    try {
      fnc_values_.resize(local_function_size);
    }
    catch (const std::exception& exc) {
      std::cout << "exception caught on resize of blocked function: " << exc.what() << '\n';
      throw(exc);
    }

    for (int linind = 0; linind < local_function_size; ++linind)
      setToZero(fnc_values_[linind]);
    double block_size = (local_function_size * sizeof(scalartype)) / 1024 / 1024 / 1024;  // gigabytes
    if (concurrency.id() == 0)
      std::cout << "Blocked function " << vectorToString(dmn.get_leaf_domain_sizes()) << '\n'
                << "allocated: " << block_size << "GB\n";
  }
}

template <typename scalartype, class domain, DistType DT>
function<scalartype, domain, DT>& function<scalartype, domain, DT>::operator=(
    const function<scalartype, domain, DT>& other) {
  if (this != &other) {
    if constexpr (dist == DistType::NONE) {
      if (dmn.get_size() != other.dmn.get_size() || size() != other.size()) {
        // Domain had not been initialized when the functions were created.
        // Reset this function and check again.
        reset();

        if (dmn.get_size() != other.dmn.get_size() || size() != other.size())
          // The other function has not been resetted after the domain was initialized.
          throw std::logic_error("Copy assignment from a not yet resetted function.");
      }
    }
    else if constexpr (dist == DistType::BLOCKED || dist == DistType::LINEAR) {
      Nb_sbdms = other.dmn.get_leaf_domain_sizes().size();
      start_ = other.start_;
      end_ = other.end_;
      fnc_values_.resize(other.size(), {});
    }
    fnc_values_ = other.fnc_values_;
  }
  return *this;
}

template <typename Scalar, class domain, DistType DT>
template <typename Scalar2>
inline function<Scalar, domain, DT>& function<Scalar, domain, DT>::operator=(
    const function<Scalar2, domain, DT>& other) {
  if (this != &other) {
    if constexpr (dist == DistType::NONE) {
      if (size() != other.size()) {
        throw(std::logic_error("Function size does not match."));
      }
    }
    else if constexpr (dist == DistType::LINEAR || dist == DistType::BLOCKED) {
      Nb_sbdms = other.dmn.get_leaf_domain_sizes().size();
      start_ = other.start_;
      end_ = other.end_;
      fnc_values_.resize(other.size(), {});
    }
    fnc_values_ = other.fnc_values_;
  }
  return *this;
}

template <typename scalartype, class domain, DistType DT>
inline function<scalartype, domain, DT>& function<scalartype, domain, DT>::operator=(
    function<scalartype, domain, DT>&& other) {
  if (this != &other) {
    if constexpr (dist == DistType::NONE) {
      if (dmn.get_size() != other.dmn.get_size()) {
        // Domain had not been initialized when the functions were created.
        // Reset this function and check again.
        reset();

        if (dmn.get_size() != other.dmn.get_size())
          // The other function has not been resetted after the domain was initialized.
          throw std::logic_error("Move assignment from a not yet resetted function.");
      }
    }
    else if constexpr (dist == DistType::LINEAR || dist == DistType::BLOCKED) {
      Nb_sbdms = other.dmn.get_leaf_domain_sizes().size();
      start_ = other.start_;
      end_ = other.end_;
    }
    fnc_values_ = std::move(other.fnc_values_);
  }

  return *this;
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::reset() {
  if constexpr (dist != DistType::NONE)
    throw std::logic_error("a distributed function must be reset with a concurrency reference");

  dmn.reset();
  const std::size_t nb_elements = dmn.get_size();
  try {
    fnc_values_.resize(nb_elements);
  }
  catch (const std::exception& exc) {
    std::cout << "exception caught on resize of non distributed function: " << exc.what() << '\n';
    throw(exc);
  }

  for (int linind = 0; linind < nb_elements; ++linind)
    setToZero(fnc_values_[linind]);
  start_ = 0;
  end_ = dmn.get_size();
}

template <typename scalartype, class domain, DistType DT>
template <class Concurrency>
void function<scalartype, domain, DT>::reset(const Concurrency& concurrency) {
  dmn.reset();
  if constexpr (dist == DistType::NONE) {
    const std::size_t nb_elements = dmn.get_size();
    try {
      fnc_values_.resize(nb_elements);
    }
    catch (const std::exception& exc) {
      std::cout << "exception caught on resize of non distributed function: " << exc.what() << '\n';
      throw(exc);
    }

    for (int linind = 0; linind < nb_elements; ++linind)
      setToZero(fnc_values_[linind]);
    start_ = 0;
    end_ = dmn.get_size();
  }
  else if constexpr (dist == DistType::LINEAR) {
    size_t local_function_size = calcDistribution();
    try {
      fnc_values_.resize(local_function_size);
    }
    catch (const std::exception& exc) {
      std::cout << "exception caught on resize of linearly distributed function: " << exc.what()
                << '\n';
      throw(exc);
    }
    for (int linind = 0; linind < local_function_size; ++linind)
      setToZero(fnc_values_[linind]);
  }
  else if constexpr (dist == DistType::BLOCKED) {
    size_t local_function_size = calcBlocking(concurrency);
    // Ok this can be a blocked function so we finally resize i.e. allocate.
    try {
      fnc_values_.resize(local_function_size);
    }
    catch (const std::exception& exc) {
      std::cout << "exception caught on resize of blocked function: " << exc.what() << '\n';
      throw(exc);
    }
    for (int linind = 0; linind < local_function_size; ++linind)
      setToZero(fnc_values_[linind]);
    double block_size = (local_function_size * sizeof(scalartype)) / 1024 / 1024 / 1024;  // gigabytes
    if (concurrency.id == 0) {
      std::cout << "Blocked function " << vectorToString(dmn.get_leaf_domain_sizes()) << '\n'
                << "allocated: " << block_size << "GB\n";
    }
  }
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::linind_2_subind(int linind, int* subind) const {
  auto& size_sbdm = dmn.get_leaf_domain_sizes();
  for (size_t i = 0; i < size_sbdm.size(); ++i) {
    subind[i] = linind % size_sbdm[i];
    linind = (linind - subind[i]) / size_sbdm[i];
  }
}

template <typename scalartype, class domain, DistType DT>
template <std::size_t N>
void function<scalartype, domain, DT>::linind_2_subind(int linind, std::array<int, N>& subind) const {
  assert(N == dmn.get_Nb_branch_domains() || dmn.get_Nb_leaf_domains());
  auto& size_sbdm = (N == dmn.get_Nb_branch_domains() ? dmn.get_branch_domain_sizes()
                                                      : dmn.get_leaf_domain_sizes());
  for (size_t i = 0; i < size_sbdm.size(); ++i) {
    subind[i] = linind % size_sbdm[i];
    linind = (linind - subind[i]) / size_sbdm[i];
  }
}

// TODO: Resize vector if necessary.
template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::linind_2_subind(int linind, std::vector<int>& subind) const {
  auto& size_sbdm = dmn.get_leaf_domain_sizes();
  assert(int(subind.size()) == Nb_sbdms);
  for (size_t i = 0; i < Nb_sbdms; ++i) {
    subind[i] = linind % size_sbdm[i];
    linind = (linind - subind[i]) / size_sbdm[i];
  }
}

template <typename scalartype, class domain, DistType DT>
std::vector<size_t> function<scalartype, domain, DT>::linind_2_subind(int linind) const {
  std::vector<size_t> subind(Nb_sbdms);
  auto& size_sbdm = dmn.get_leaf_domain_sizes();
  if constexpr (dist == DistType::NONE) {
    for (int i = 0; i < Nb_sbdms; ++i) {
      subind[i] = linind % size_sbdm[i];
      linind = (linind - subind[i]) / size_sbdm[i];
    }
  }
  else if constexpr (dist == DistType::LINEAR) {
    std::cout << "linind:" << linind << '\n';
    throw std::runtime_error("Subindices aren't valid accessors for DistType::LINEAR");
  }
  else if constexpr (dist == DistType::BLOCKED) {
    for (int i = 0; i < int(size_sbdm.size()); ++i) {
      subind[i] = linind % size_sbdm[i];
      linind = (linind - subind[i]) / size_sbdm[i];
    }
  }
  return subind;
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::subind_2_linind(const int* const subind, int& linind) const {
  auto& step_sbdm = dmn.get_leaf_domain_steps();
  linind = 0;
  for (int i = 0; i < int(step_sbdm.size()); ++i)
    linind += subind[i] * step_sbdm[i];
}

template <typename scalartype, class domain, DistType DT>
size_t function<scalartype, domain, DT>::subind_2_linind(const std::vector<int>& subind) const {
  auto& step_sbdm = dmn.get_leaf_domain_steps();
  assert(subind.size() == step_sbdm.size());
  int linind = 0;
  for (int i = 0; i < int(step_sbdm.size()); ++i)
    linind += subind[i] * step_sbdm[i];
  return linind;
}

template <typename scalartype, class domain, DistType DT>
size_t function<scalartype, domain, DT>::branch_subind_2_linind(const std::vector<int>& subind) const {
  auto& branch_sbdm = dmn.get_branch_domain_steps();
  assert(subind.size() == branch_sbdm.size());
  int linind = 0;
  for (int i = 0; i < int(branch_sbdm.size()); ++i)
    linind += subind[i] * branch_sbdm[i];
  return linind;
}

template <typename scalartype, class domain, DistType DT>
scalartype& function<scalartype, domain, DT>::operator()(const int* const subind) {
  auto& step_sbdm = dmn.get_leaf_domain_steps();
  int linind;
  subind_2_linind(subind, linind);

  assert(linind >= 0 && linind < size());
  return fnc_values_[linind];
}

template <typename scalartype, class domain, DistType DT>
const scalartype& function<scalartype, domain, DT>::operator()(const int* const subind) const {
  int linind;
  subind_2_linind(subind, linind);

  assert(linind >= 0 && linind < size());
  return fnc_values_[linind];
}

template <typename scalartype, class domain, DistType DT>
const scalartype& function<scalartype, domain, DT>::operator()(const std::vector<int>& subind) const {
  int linind = 0;  // silence warning
  if (subind.size() == Nb_sbdms) {
    linind = subind_2_linind(subind);
    assert(linind >= 0 && linind < size());
  }
  else if (subind.size() == dmn.get_Nb_branch_domains()) {
    linind = branch_subind_2_linind(subind);
  }
  else
    throw std::runtime_error("number of indicies matches neither branches or leaves");
  return fnc_values_[linind];
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::operator+=(const function<scalartype, domain, DT>& other) {
  for (int linind = 0; linind < size(); ++linind)
    fnc_values_[linind] += other(linind);
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::operator-=(const function<scalartype, domain, DT>& other) {
  for (int linind = 0; linind < size(); ++linind)
    fnc_values_[linind] -= other(linind);
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::operator*=(const function<scalartype, domain, DT>& other) {
  for (int linind = 0; linind < size(); ++linind)
    fnc_values_[linind] *= other(linind);
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::operator/=(const function<scalartype, domain, DT>& other) {
  for (int linind = 0; linind < size(); ++linind) {
    assert(std::abs(other(linind)) > 1.e-16);
    fnc_values_[linind] /= other(linind);
  }
}

template <typename scalartype, class domain, DistType DT>
inline void function<scalartype, domain, DT>::operator=(const scalartype c) {
  for (int linind = 0; linind < size(); linind++)
    fnc_values_[linind] = c;
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::operator+=(const scalartype c) {
  for (int linind = 0; linind < size(); linind++)
    fnc_values_[linind] += c;
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::operator-=(const scalartype c) {
  for (int linind = 0; linind < size(); linind++)
    fnc_values_[linind] -= c;
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::operator*=(const scalartype c) {
  for (int linind = 0; linind < size(); linind++)
    fnc_values_[linind] *= c;
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::operator/=(const scalartype c) {
  for (int linind = 0; linind < size(); linind++)
    fnc_values_[linind] /= c;
}

template <typename scalartype, class domain, DistType DT>
bool function<scalartype, domain, DT>::operator==(const function<scalartype, domain, DT>& other) const {
  if (size() != other.size())
    // One of the function has not been resetted after the domain was initialized.
    throw std::logic_error("Comparing functions of different sizes.");

  for (int i = 0; i < size(); ++i)
    if (other(i) != fnc_values_[i])
      return false;

  return true;
}

template <typename scalartype, class domain, DistType DT>
template <typename new_scalartype>
void function<scalartype, domain, DT>::slice(const int sbdm_index, std::vector<int> subind,
                                             new_scalartype* fnc_vals) const {
  assert(sbdm_index >= 0);
  assert(sbdm_index < Nb_sbdms);

  int linind = 0;

  if (subind.size() < Nb_sbdms && subind.size() == dmn.get_Nb_branch_domains()) {
    auto& size_sbdm = dmn.get_branch_domain_sizes();
    auto& step_sbdm = dmn.get_branch_domain_steps();

    subind[sbdm_index] = 0;
    linind = branch_subind_2_linind(subind);

    for (int i = 0; i < size_sbdm[sbdm_index]; i++)
      fnc_vals[i] =
          ScalarCast<new_scalartype>::execute(fnc_values_[linind + i * step_sbdm[sbdm_index]]);
  }
  else if (subind.size() == Nb_sbdms) {
    auto& size_sbdm = dmn.get_leaf_domain_sizes();
    auto& step_sbdm = dmn.get_leaf_domain_steps();

    subind[sbdm_index] = 0;
    linind = subind_2_linind(subind);

    for (int i = 0; i < size_sbdm[sbdm_index]; i++)
      fnc_vals[i] =
          ScalarCast<new_scalartype>::execute(fnc_values_[linind + i * step_sbdm[sbdm_index]]);
  }
  else
    throw std::runtime_error(
        "function::slice can only be called with subind of size of leaf or branch domains.");
}

template <typename scalartype, class domain, DistType DT>
template <typename new_scalartype>
void function<scalartype, domain, DT>::slice(const int sbdm_index, int* subind,
                                             new_scalartype* fnc_vals) const {
  assert(sbdm_index >= 0);
  assert(sbdm_index < Nb_sbdms);

  auto& size_sbdm = dmn.get_leaf_domain_sizes();
  auto& step_sbdm = dmn.get_leaf_domain_steps();

  int linind = 0;
  subind[sbdm_index] = 0;
  subind_2_linind(subind, linind);
  for (int i = 0; i < size_sbdm[sbdm_index]; i++)
    fnc_vals[i] =
        ScalarCast<new_scalartype>::execute(fnc_values_[linind + i * step_sbdm[sbdm_index]]);
}

template <typename scalartype, class domain, DistType DT>
template <typename new_scalartype>
void function<scalartype, domain, DT>::slice(const int sbdm_index_1, const int sbdm_index_2,
                                             int* subind, new_scalartype* fnc_vals) const {
  assert(sbdm_index_1 >= 0);
  assert(sbdm_index_2 >= 0);
  assert(sbdm_index_1 < Nb_sbdms);
  assert(sbdm_index_2 < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index_1] = 0;
  subind[sbdm_index_2] = 0;
  subind_2_linind(subind, linind);

  auto& size_sbdm = dmn.get_leaf_domain_sizes();
  auto& step_sbdm = dmn.get_leaf_domain_steps();

  int size_sbdm_1 = size_sbdm[sbdm_index_1];
  int size_sbdm_2 = size_sbdm[sbdm_index_2];

  int step_sbdm_1 = step_sbdm[sbdm_index_1];
  int step_sbdm_2 = step_sbdm[sbdm_index_2];

  new_scalartype* fnc_ptr_left = NULL;

  for (int j = 0; j < size_sbdm_2; j++) {
    fnc_ptr_left = &fnc_vals[0 + j * size_sbdm_1];
    const new_scalartype* fnc_ptr_right{fnc_values_.data() + linind + j * step_sbdm_2};

    for (int i = 0; i < size_sbdm_1; i++)
      fnc_ptr_left[i] = fnc_ptr_right[i * step_sbdm_1];
    //       fnc_vals[i+j*size_sbdm[sbdm_index_1]] = fnc_values_[linind + i*step_sbdm[sbdm_index_1]
    //       + j*step_sbdm[sbdm_index_2]];
  }
}

template <typename scalartype, class domain, DistType DT>
template <typename new_scalartype>
void function<scalartype, domain, DT>::slice(const int sbdm_index_1, const int sbdm_index_2,
                                             std::vector<int> subind, new_scalartype* fnc_vals) const {
  assert(sbdm_index_1 >= 0);
  assert(sbdm_index_2 >= 0);
  assert(sbdm_index_1 < Nb_sbdms);
  assert(sbdm_index_2 < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index_1] = 0;
  subind[sbdm_index_2] = 0;

  if (subind.size() < Nb_sbdms && subind.size() == dmn.get_Nb_branch_domains()) {
    linind = branch_subind_2_linind(subind);

    auto& size_sbdm = dmn.get_branch_domain_sizes();
    auto& step_sbdm = dmn.get_branch_domain_steps();

    int size_sbdm_1 = size_sbdm[sbdm_index_1];
    int size_sbdm_2 = size_sbdm[sbdm_index_2];

    int step_sbdm_1 = step_sbdm[sbdm_index_1];
    int step_sbdm_2 = step_sbdm[sbdm_index_2];

    new_scalartype* fnc_ptr_left = NULL;

    for (int j = 0; j < size_sbdm_2; j++) {
      fnc_ptr_left = &fnc_vals[0 + j * size_sbdm_1];
      const new_scalartype* fnc_ptr_right{fnc_values_.data() + linind + j * step_sbdm_2};

      for (int i = 0; i < size_sbdm_1; i++)
        fnc_ptr_left[i] = fnc_ptr_right[i * step_sbdm_1];
      //       fnc_vals[i+j*size_sbdm[sbdm_index_1]] = fnc_values_[linind + i*step_sbdm[sbdm_index_1]
      //       + j*step_sbdm[sbdm_index_2]];
    }
  }
  else {
    linind = subind_2_linind(subind);

    auto& size_sbdm = dmn.get_leaf_domain_sizes();
    auto& step_sbdm = dmn.get_leaf_domain_steps();
    int size_sbdm_1 = size_sbdm[sbdm_index_1];
    int size_sbdm_2 = size_sbdm[sbdm_index_2];

    int step_sbdm_1 = step_sbdm[sbdm_index_1];
    int step_sbdm_2 = step_sbdm[sbdm_index_2];

    new_scalartype* fnc_ptr_left = NULL;

    for (int j = 0; j < size_sbdm_2; j++) {
      fnc_ptr_left = &fnc_vals[0 + j * size_sbdm_1];
      const new_scalartype* fnc_ptr_right{fnc_values_.data() + linind + j * step_sbdm_2};

      for (int i = 0; i < size_sbdm_1; i++)
        fnc_ptr_left[i] = fnc_ptr_right[i * step_sbdm_1];
      //       fnc_vals[i+j*size_sbdm[sbdm_index_1]] = fnc_values_[linind + i*step_sbdm[sbdm_index_1]
      //       + j*step_sbdm[sbdm_index_2]];
    }
  }
}

template <typename scalartype, class domain, DistType DT>
template <typename new_scalartype>
void function<scalartype, domain, DT>::distribute(const int sbdm_index, std::vector<int> subind,
                                                  const new_scalartype* fnc_vals) {
  assert(sbdm_index >= 0);
  assert(sbdm_index < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index] = 0;

  if (subind.size() < Nb_sbdms && subind.size() == dmn.get_Nb_branch_domains()) {
    linind = branch_subind_2_linind(subind);

    auto& size_sbdm = dmn.get_branch_domain_sizes();
    auto& step_sbdm = dmn.get_branch_domain_steps();

    for (int i = 0; i < size_sbdm[sbdm_index]; i++)
      fnc_values_[linind + i * step_sbdm[sbdm_index]] = ScalarCast<scalartype>::execute(fnc_vals[i]);
  }
  else {
    linind = subind_2_linind(subind);

    auto& size_sbdm = dmn.get_leaf_domain_sizes();
    auto& step_sbdm = dmn.get_leaf_domain_steps();

    for (int i = 0; i < size_sbdm[sbdm_index]; i++)
      fnc_values_[linind + i * step_sbdm[sbdm_index]] = ScalarCast<scalartype>::execute(fnc_vals[i]);
  }
}

template <typename scalartype, class domain, DistType DT>
template <typename new_scalartype>
void function<scalartype, domain, DT>::distribute(const int sbdm_index, int* subind,
                                                  const new_scalartype* fnc_vals) {
  assert(sbdm_index >= 0);
  assert(sbdm_index < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index] = 0;
  subind_2_linind(subind, linind);

  auto& size_sbdm = dmn.get_leaf_domain_sizes();
  auto& step_sbdm = dmn.get_leaf_domain_steps();

  for (int i = 0; i < size_sbdm[sbdm_index]; i++)
    fnc_values_[linind + i * step_sbdm[sbdm_index]] = ScalarCast<scalartype>::execute(fnc_vals[i]);
}

template <typename scalartype, class domain, DistType DT>
template <typename new_scalartype>
void function<scalartype, domain, DT>::distribute(const int sbdm_index_1, const int sbdm_index_2,
                                                  std::vector<int> subind,
                                                  const new_scalartype* fnc_vals) {
  assert(sbdm_index_1 >= 0);
  assert(sbdm_index_2 >= 0);
  assert(sbdm_index_1 < Nb_sbdms);
  assert(sbdm_index_2 < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index_1] = 0;
  subind[sbdm_index_2] = 0;

  if (subind.size() < Nb_sbdms && subind.size() == dmn.get_Nb_branch_domains()) {
    linind = branch_subind_2_linind(subind);

    auto& size_sbdm = dmn.get_branch_domain_sizes();
    auto& step_sbdm = dmn.get_branch_domain_steps();

    for (int i = 0; i < size_sbdm[sbdm_index_1]; i++)
      for (int j = 0; j < size_sbdm[sbdm_index_2]; j++)
        fnc_values_[linind + i * step_sbdm[sbdm_index_1] + j * step_sbdm[sbdm_index_2]] =
            fnc_vals[i + j * size_sbdm[sbdm_index_1]];
  }
}

template <typename scalartype, class domain, DistType DT>
template <typename new_scalartype>
void function<scalartype, domain, DT>::distribute(const int sbdm_index_1, const int sbdm_index_2,
                                                  int* subind, const new_scalartype* fnc_vals) {
  assert(sbdm_index_1 >= 0);
  assert(sbdm_index_2 >= 0);
  assert(sbdm_index_1 < Nb_sbdms);
  assert(sbdm_index_2 < Nb_sbdms);

  int linind = 0;
  subind[sbdm_index_1] = 0;
  subind[sbdm_index_2] = 0;
  subind_2_linind(subind, linind);

  auto& size_sbdm = dmn.get_leaf_domain_sizes();
  auto& step_sbdm = dmn.get_leaf_domain_steps();

  for (int i = 0; i < size_sbdm[sbdm_index_1]; i++)
    for (int j = 0; j < size_sbdm[sbdm_index_2]; j++)
      fnc_values_[linind + i * step_sbdm[sbdm_index_1] + j * step_sbdm[sbdm_index_2]] =
          fnc_vals[i + j * size_sbdm[sbdm_index_1]];
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::print_fingerprint(std::ostream& stream) const {
  stream << "****************************************\n";
  stream << "function: " << name_ << "\n";
  stream << "****************************************\n";

  stream << "# subdomains: " << Nb_sbdms << "\n";
  util::print_type<domain>::print(stream);
  stream << "\n";

  auto& size_sbdm = dmn.get_leaf_domain_sizes();

  stream << "size of subdomains:";
  for (int i = 0; i < Nb_sbdms; ++i)
    stream << "  " << size_sbdm[i];
  stream << "\n";

  stream << "# elements: " << size() << "\n";
  stream << "memory: " << size() * sizeof(scalartype) / (1024. * 1024.) << " MiB\n";
  stream << "****************************************\n" << std::endl;
}

template <typename scalartype, class domain, DistType DT>
void function<scalartype, domain, DT>::print_elements(std::ostream& stream) const {
  stream << "****************************************\n";
  stream << "function: " << name_ << "\n";
  stream << "****************************************\n";

  std::vector<int> subind(Nb_sbdms);
  for (int lindex = 0; lindex < size(); ++lindex) {
    linind_2_subind(lindex, subind);
    for (int index : subind)
      stream << index << "\t";
    stream << " \t" << fnc_values_[lindex] << "\n";
  }

  stream << "****************************************\n" << std::endl;
}

template <typename scalartype, class domain, DistType DT>
template <typename concurrency_t>
int function<scalartype, domain, DT>::get_buffer_size(const concurrency_t& concurrency) const {
  int result = 0;
  result += concurrency.get_buffer_size(*this);
  return result;
}

template <typename scalartype, class domain, DistType DT>
template <class concurrency_t>
void function<scalartype, domain, DT>::pack(const concurrency_t& concurrency, char* buffer,
                                            const int buffer_size, int& position) const {
  concurrency.pack(buffer, buffer_size, position, *this);
}

template <typename scalartype, class domain, DistType DT>
template <class concurrency_t>
void function<scalartype, domain, DT>::unpack(const concurrency_t& concurrency, char* buffer,
                                              const int buffer_size, int& position) {
  concurrency.unpack(buffer, buffer_size, position, *this);
}

template <typename scalartype, class domain, DistType DT>
template <class Concurrency>
function<scalartype, domain, DT> function<scalartype, domain, DT>::gather(
    const Concurrency& concurrency) const {
  function result(name_, concurrency);
  concurrency.gather(*this, result, concurrency);
  return result;
}
}  // namespace func
}  // namespace dca
#endif  // DCA_FUNCTION_FUNCTION_HPP
