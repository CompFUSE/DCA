/*
 * dmn_variadic.h
 *
 *  Created on: Mar 10, 2016
 *      Author: John Biddiscombe
 */

#ifndef DMN_VARIADIC_H_
#define DMN_VARIADIC_H_

#include <sys/time.h>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <assert.h>
#include <utility>
#include <tuple>
#include <type_traits>
#include <iterator>
//
#include "domain.h"
#include "type_list.h"
#include "type_list_definitions.h"

//----------------------------------------------------------------------------
// Variadic domain class with a list of sub-domains
//----------------------------------------------------------------------------
template<typename... domain_list>
class dmn_variadic: public domain {

public:

//    typedef typename Append<domain_typelist_0,
//        typename Append<domain_typelist_1, domain_typelist_2>::Result>::Result this_type;

    /// constructor, responsible for initializing all domain indexing arrays
    dmn_variadic();

    /// return the size allocated by this domain.
    /// includes all space for all subdomains
    static int& dmn_size() {
        static int size = -1;
        std::cout << "Returning domain size " << size << std::endl << std::endl;
        return size;
    }

    /// reset the domain back to the original state at creation time
    void reset();

    /// primary operator to access elements by index.
    ///
    template <typename ...Args>
    int operator()(Args &&... args);

protected:

    template<typename... Args>
    int index_lookup(
        std::integral_constant<bool, true>::type,
        int branch_i0, Args... args
    );

    template<typename... Args>
    int index_lookup(
        std::integral_constant<bool, false>::type,
        int leaf_i0, Args... args
    );

    template <typename Tuple, std::size_t... DomainIndex>
    std::vector<int>
    init_branch_domain_sizes(Tuple &t, std::index_sequence<DomainIndex...>);

    template <typename Tuple, std::size_t... DomainIndex>
    std::vector<int>
    init_leaf_domain_sizes(Tuple &domaintuple, std::index_sequence<DomainIndex...>);


protected:
    std::tuple<domain_list...> domains;
};

//----------------------------------------------------------------------------
// utility function : takes a variadic list or arguments and returns nothing.
// It is used to convert a variadic pack expansion into a function so that
// arbitrary functions can be called with a pack expansion and drop the results.
//----------------------------------------------------------------------------
template<typename...Ts>
void ignore_returnvalues(Ts&&...) {}

//----------------------------------------------------------------------------
// Get the branch domain sizes for each of the domain template arguments
// returns a vector of the sizes
//----------------------------------------------------------------------------
template <typename... domain_list>
template <typename Tuple, std::size_t... DomainIndex>
std::vector<int>
dmn_variadic<domain_list...>::init_branch_domain_sizes(Tuple &domaintuple, std::index_sequence<DomainIndex...>)
{
    return { (std::get<DomainIndex>(domaintuple).get_size())... };
}

//----------------------------------------------------------------------------
// Get the number of leaf domains for each of the domain template arguments
// returns a vector of the counts
//----------------------------------------------------------------------------
template <typename... domain_list>
template <typename Tuple, std::size_t... DomainIndex>
std::vector<int>
dmn_variadic<domain_list...>::init_leaf_domain_sizes(Tuple &domaintuple, std::index_sequence<DomainIndex...>)
{
    return { (std::get<DomainIndex>(domaintuple).get_Nb_leaf_domains())... };
}

//----------------------------------------------------------------------------
// Invoke a function for each element of a tuple
//----------------------------------------------------------------------------
namespace detail
{
    template<typename T, typename F, std::size_t... Indices>
    void for_each(T&& t, F &&f, std::index_sequence<Indices...>)
    {
        auto l = { (f(std::get<Indices>(t)), 0)... };
    }
}

template<typename... Ts, typename F>
void for_each_in_tuple(std::tuple<Ts...> & t, F &&f)
{
    detail::for_each(t, std::forward<F>(f), std::make_index_sequence<sizeof...(Ts)>{});
}

//----------------------------------------------------------------------------
// helper function : for each domain, get the number of leaf domains
//----------------------------------------------------------------------------
struct leaf_domain_size_helper
{
    leaf_domain_size_helper(std::vector<int> &dest) : destination(dest) {}
    //
    template<typename Domain>
    void operator () (Domain &&d)
    {
        std::vector<int> temp = d.get_leaf_domain_sizes();
        std::cout << "Here with a domain of size " << d.get_size() << std::endl;
        std::copy(std::begin(temp), std::end(temp), std::ostream_iterator<int>(std::cout, ","));
        std::copy(std::begin(temp), std::end(temp), std::back_inserter(destination));
        std::cout << std::endl;
    }
    //
    std::vector<int> &destination;
};

//----------------------------------------------------------------------------
// modified an index sequence, offset/length
// for indexing operations we multiply arguments by branch/leaf step sizes
//----------------------------------------------------------------------------
template <std::size_t O, std::size_t ... Is>
std::index_sequence<(O + Is)...> add_offset(std::index_sequence<Is...>)
{
    return {};
}

template <std::size_t O, std::size_t N>
auto make_index_sequence_with_offset()
{
    return add_offset<O>(std::make_index_sequence<N>{});
}

//----------------------------------------------------------------------------
// index multiplication
//----------------------------------------------------------------------------
template<typename T>
T sum(T v) {
  return v;
}

template<typename T, typename... Args>
T sum(T first, Args... args) {
  return first + adder(args...);
}


template <typename ...Args, std::size_t ... Is>
int multiply_offsets(const std::vector<int> &multipliers, Args &&... offsets, std::index_sequence<Is...>)
{
    return sum((offsets*std::get<Is>(multipliers))...);
}
//----------------------------------------------------------------------------
// Constructor implementation
//----------------------------------------------------------------------------
template<typename... domain_list>
dmn_variadic<domain_list...>::dmn_variadic() : domain()
{
    dmn_variadic<domain_list...>::reset();
}

//----------------------------------------------------------------------------
// Get the number of leaf domains for each of the domain template arguments
// returns a vector of the counts
//----------------------------------------------------------------------------
template<typename... domain_list>
void dmn_variadic<domain_list...>::reset()
{
    domain::reset();

    // create an index sequence that indexes the domains we are templated on
    std::index_sequence_for<domain_list...> indices;

    // initialize a vector from the size of each top level domain (branch domain)
    branch_domain_sizes = init_branch_domain_sizes(domains,indices);

    std::cout << "Creating " << __PRETTY_FUNCTION__ << " " << branch_domain_sizes.size() << std::endl << "domain sizes : ";
    std::copy(branch_domain_sizes.begin(), branch_domain_sizes.end(), std::ostream_iterator<int>(std::cout,","));
    std::cout << std::endl;

    branch_domain_steps.resize(branch_domain_sizes.size(), 1);
    for (size_t i = 0; i < branch_domain_sizes.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            branch_domain_steps[i] *= branch_domain_sizes[j];
        }
    }
    std::cout << "Steps ";
    std::copy(branch_domain_steps.begin(), branch_domain_steps.end(), std::ostream_iterator<int>(std::cout,","));
    std::cout << std::endl;

    // generate the leaf domain sizes from each sub domain
    auto leaf_stuff = init_leaf_domain_sizes(domains, indices);
    std::cout << "Leaf stuff ";
    std::copy(leaf_stuff.begin(), leaf_stuff.end(), std::ostream_iterator<int>(std::cout, ","));
    std::cout << std::endl;

    for_each_in_tuple(domains, leaf_domain_size_helper(leaf_domain_sizes));
    std::cout << "leaf domain size ";
    std::copy(leaf_domain_sizes.begin(), leaf_domain_sizes.end(), std::ostream_iterator<int>(std::cout,","));
    std::cout << std::endl;

    leaf_domain_steps.resize(leaf_domain_sizes.size(), 1);
    for (size_t i = 0; i < leaf_domain_sizes.size(); i++)
        for (size_t j = 0; j < i; j++)
            leaf_domain_steps[i] *= leaf_domain_sizes[j];

    std::cout << "leaf step size ";
    std::copy(leaf_domain_steps.begin(), leaf_domain_steps.end(), std::ostream_iterator<int>(std::cout,","));
    std::cout << std::endl;

    size = 1;

    for (int i = 0; i < sizeof...(domain_list); i++)
        size *= branch_domain_sizes[i];

    dmn_size() = size;
    std::cout << "Domain size is " << size << std::endl;
}

//----------------------------------------------------------------------------
// indexing operator : access elements of the domain
// if sizeof(args) < sizeof(domains) : error
// if sizeof(args) == sizeof(domains) : call index via branch domains
// if sizeof(args) == sizeof(leaf domains) : call index via leaf domains
//----------------------------------------------------------------------------
template<typename... domain_list>
template <typename ...Args>
int dmn_variadic<domain_list...>::operator()(Args&&... args) {
    static_assert( sizeof...(Args) < sizeof...(domain_list), "not enough args");
    return index_lookup(
        std::integral_constant<bool, (sizeof...(Args) == sizeof...(domain_list))>::type,
        std::forward<Args>(args)...
    );
}

//----------------------------------------------------------------------------
// indexing operator : access elements of the domain via branches
// index_lookup is overloaded on std::integral_constant<bool, true>::type so
// that if sizeof...(Args) == sizeof...(domain_list) then this is called
//----------------------------------------------------------------------------
template<typename... domain_list>
template<typename... Args>
int dmn_variadic<domain_list...>::index_lookup(
    std::integral_constant<bool, true>::type,
    int branch_i0, Args... branch_indices)
{
    std::cout << "Branch overload " << std::endl;

    // create an index sequence that indexes the domains we are templated on
    std::index_sequence_for<domain_list...> indices;

    // create an index sequence starting from 1, with length sizeof...(args)-1
    auto seq = make_index_sequence_with_offset<1, sizeof...(Args)>();

    int N = multiply_offsets(branch_domain_steps, std::forward<Args...>(branch_indices)..., indices);
/*
    return branch_i0 + sum( * branch_indices)...;

    branch_index_helper bih(branch_domain_sizes);
    for_each_in_tuple(domains, bih);
    return bih.index;
*/



    /*
        assert(branch_domain_sizes.size() == 3);
        assert(branch_i0 >= 0 && branch_i0 < branch_domain_sizes[0]);
        assert(branch_i1 >= 0 && branch_i1 < branch_domain_sizes[1]);
        assert(branch_i2 >= 0 && branch_i2 < branch_domain_sizes[2]);

        return branch_i0
            + branch_domain_steps[1] * branch_i1
            + branch_domain_steps[2] * branch_i2;
    */
}

//----------------------------------------------------------------------------
// indexing operator : access elements of the domain via leaf domains
// index_lookup is overloaded on std::integral_constant<bool, true>::type so
// that if sizeof...(Args) == sizeof...(domain_list) then this is called
//----------------------------------------------------------------------------
template<typename... domain_list>
template<typename... Args>
int dmn_variadic<domain_list...>::index_lookup(
    std::integral_constant<bool, false>::type,
    int leaf_i0, Args... leaf_indices)
{
    std::cout << "Subdomain overload " << std::endl;
/*
    assert(leaf_domain_sizes.size() == 4);
    assert(sbdmn_i0 >= 0 && sbdmn_i0 < leaf_domain_sizes[0]);
    assert(sbdmn_i1 >= 0 && sbdmn_i1 < leaf_domain_sizes[1]);
    assert(sbdmn_i2 >= 0 && sbdmn_i2 < leaf_domain_sizes[2]);
    assert(sbdmn_i3 >= 0 && sbdmn_i3 < leaf_domain_sizes[3]);

    return sbdmn_i0
        + leaf_domain_steps[1] * sbdmn_i1
        + leaf_domain_steps[2] * sbdmn_i2
        + leaf_domain_steps[3] * sbdmn_i3;
*/
}

#endif

