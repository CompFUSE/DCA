// Copyright (C) 2020 ETH Zurich
// Copyright (C) 2020 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Peter Doak (doakpw@ornl.gov)
//
// This is a partial specialization for blocked (i.e. regular blocks of tensor)
// of function.hpp
//
// TODO: Remame template parameter scalartype --> ElementType.

#ifndef DCA_FUNCTION_FUNCTION_DIST_BLOCKED_HPP
#define DCA_FUNCTION_FUNCTION_DIST_BLOCKED_HPP

namespace dca {
namespace func {
// dca::func::

template <typename scalartype, class domain>
class function<scalartype, domain, DistType::BLOCKED> {
  static const std::string default_name_;
  // This is a clumsy work around of the fact a constexpr DT shadows the template parameter in the
  // generic class so it can't be defined there.  This breaks later if constexpr so that needs to be
  // on DISTTYPE, but it is convenient for the specialized template classes declarations to "almost" match the
  // general class declaration.
  static constexpr auto DISTTYPE = DistType::BLOCKED;
  static constexpr auto DT = DistType::BLOCKED;

public:
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
  function(const function<scalartype, domain>& other);
  // Same as above, but with name = 'name'.
  function(const function<scalartype, domain>& other, const std::string& name) : function(other) {
    name_ = name;
  }

  // Move constructor
  // Constructs the function with elements and name of other using move semantics.
  // Precondition: The other function has been resetted, if the domain had been initialized after
  //               the other function's construction.
  // Postcondition: The other function is in a non-specified state.
  function(function<scalartype, domain>&& other);
  // Same as above, but with name = 'name'.
  function(function<scalartype, domain>&& other, const std::string& name)
      : function(std::move(other)) {
    name_ = name;
  }

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

  std::size_t get_start() const {
    return start_;
  }
  /** end in sense of last index not  1 past.
   */
  std::size_t get_end() const {
    return end_ - 1;
  }

  std::vector<int> get_start_subindex() const {
    return linind_2_subind(start_);
  }
  /** end in sense of last subindex
   */
  std::vector<int> get_end_subindex() const {
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

  // Returns the size of the leaf domain with the given index.
  // Does not return function values!
  int operator[](const int index) const {
    throw std::runtime_error("function::operator[] broken for concurrency id > 0, therefore do not use for blocked functions");
    return size_sbdm[index-start_];
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
  // modern RVO version
  std::vector<int> linind_2_subind(int linind) const;

  // Computes the linear index for the given subindices of the leaf domains.
  // Precondition: subind stores the the subindices of all LEAF domains.
  // TODO: Use std::array or std::vector to be able to check the size of subind.
  void subind_2_linind(const int* subind, int& linind) const;

  // using standard vector and avoiding returning argument
  int subind_2_linind(const std::vector<int>& subind) const;

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
    assert(ind >= 0 && ind < size());
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
    assert(linind >= 0 && linind < size());
    return fnc_values_[linind];
  }
  template <typename T>
  const scalartype& operator()(const T linind) const {
    static_assert(std::is_integral<T>::value, "Index linind must be an integer.");
    assert(linind >= 0 && linind < size());
    return fnc_values_[linind];
  }

  template <typename... Ts>
  scalartype& operator()(const Ts... subindices) {
    // We need to cast all indices to the same type for dmn_variadic.
    return fnc_values_[dmn(static_cast<int>(subindices)...)];
  }
  template <typename... Ts>
  const scalartype& operator()(const Ts... subindices) const {
    return fnc_values_[dmn(static_cast<int>(subindices)...)];
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

  // Gather a function that was initialized as distributed.
  // Precondition: concurrency must be the same object used during construction.
  template <class Concurrency>
  function gather(const Concurrency& concurrency) const;

private:
  std::string name_;
  std::string function_type;

  domain dmn;  // TODO: Remove domain object?

  // The subdomains (sbdmn) represent the leaf domains, not the branch domains.
  int Nb_sbdms;
  const std::vector<std::size_t>& size_sbdm;  // TODO: Remove?
  const std::vector<std::size_t>& step_sbdm;  // TODO: Remove?

  std::vector<scalartype> fnc_values_;

  std::size_t start_;
  std::size_t end_;
};

template <typename scalartype, class domain>
template <class Concurrency>
function<scalartype, domain, DistType::BLOCKED>::function(const std::string& name,
                                                          const Concurrency& concurrency)
    : name_(name),
      function_type(__PRETTY_FUNCTION__),
      dmn(),
      Nb_sbdms(dmn.get_leaf_domain_sizes().size()),
      size_sbdm(dmn.get_leaf_domain_sizes()),
      step_sbdm(dmn.get_leaf_domain_steps()) {
  const std::size_t conc_size = concurrency.number_of_processors();
  const std::size_t nb_elements = dca::util::ceilDiv(dmn.get_size(), conc_size);

  fnc_values_.resize(nb_elements);

  for (int linind = 0; linind < nb_elements; ++linind)
    setToZero(fnc_values_[linind]);

  int my_concurrency_id = concurrency.id();
  int my_concurrency_size = concurrency.number_of_processors();

  std::size_t local_function_size =
      dca::util::ceilDiv(dmn.get_size(), std::size_t(my_concurrency_size));
  start_ = local_function_size * my_concurrency_id;
  end_ = std::min(dmn.get_size(), start_ + local_function_size);
  // This is a necessary but not sufficient proof of "regular blocking"
  if (end_ != start_ + local_function_size)
    throw std::runtime_error(
        "Blocked concurrency is not possible if concurrency size is not blockwise divisor of "
        "functions size");
}

template <typename scalartype, class domain>
std::vector<int> function<scalartype, domain, DistType::BLOCKED>::linind_2_subind(int linind) const {
  std::vector<int> subind(Nb_sbdms);
  for (int i = 0; i < int(size_sbdm.size()); ++i) {
    subind[i] = linind % size_sbdm[i];
    linind = (linind - subind[i]) / size_sbdm[i];
  }
  return subind;
}

}  // namespace func
}  // namespace dca
#endif
