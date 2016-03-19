//-*-C++-*-

#ifndef DMN_VARIADIC_H_
#define DMN_VARIADIC_H_

#include <sys/time.h>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <assert.h>
#include <utility>
#include <type_traits>
#include "domain.h"
#include "type_list.h"
#include "type_list_definitions.h"

//----------------------------------------------------------------------------
// Indices trick, from http://loungecpp.wikidot.com/tips-and-tricks%3aindices
//----------------------------------------------------------------------------

template<std::size_t... Indices>
struct indices {
};

// this assumes that function type F returns int
template<typename Tuple, std::size_t... Indices, typename F>
std::array<int, std::tuple_size<Tuple>::value>
apply(Tuple &&t, indices<Indices...>, F &&f) {
    return std::array<int, std::tuple_size<Tuple>::value> {
        {
            f(std::get<Indices>(std::forward<Tuple>(t)))...
        }
    };
}

template<std::size_t N, std::size_t... Is>
struct build_indices: build_indices<N - 1, N - 1, Is...> {
};

template<std::size_t... Is>
struct build_indices<0, Is...> : indices<Is...> {
};

template<typename Tuple>
using IndicesFor = build_indices<std::tuple_size<std::decay<Tuple>>::value>;


//----------------------------------------------------------------------------
// C++14 integer sequence
//----------------------------------------------------------------------------
template <typename Tuple, size_t... Indices>
std::array<int, sizeof...(Indices)>
call_f_detail(Tuple& tuple, std::index_sequence<Indices...> ) {
    return { f(std::get<Indices>(tuple))... };
}

template <typename Tuple>
std::array<int, std::tuple_size<Tuple>::value>
call_f(Tuple& tuple) {
    return call_f_detail(tuple,
        // make the sequence type sequence<0, 1, 2, ..., N-1>
        std::make_index_sequence<std::tuple_size<Tuple>::value>{}
        );
}
//----------------------------------------------------------------------------
using namespace TL;

template<typename... domain_list>
class dmn_variadic: public domain
{

public:

//    typedef typename Append<domain_typelist_0,
//        typename Append<domain_typelist_1, domain_typelist_2>::Result>::Result this_type;

    dmn_variadic();

    static int& dmn_size() {
        static int size = -1;
        return size;
    }

    void reset();

    template <typename ...Args>
    int operator()(Args &&... args);

    template <typename ...Args>
    int index_lookup(
        std::integral_constant<bool, true>::type,
        Args &&... args
    );

    template <typename ...Args>
    int index_lookup(
        std::integral_constant<bool, false>::type,
        Args &&... args
    );

/*
    template <bool, typename ...Args>
    struct operator_impl{
        int operator()(Args &&... args);
    };
*/
protected:
    std::tuple<domain_list...> domains;
};


template<typename... domain_list>
dmn_variadic<domain_list...>::dmn_variadic() : domain()
{
    std::index_sequence_for<domain_list...> indices;
    branch_domain_sizes = {
        std::get<indices>(domains).get_size()
    };

    branch_domain_steps.resize(branch_domain_sizes.size(), 1);
    for (size_t i = 0; i < branch_domain_sizes.size(); i++)
        for (size_t j = 0; j < i; j++)
            branch_domain_steps[i] *= branch_domain_sizes[j];
/*
    leaf_domain_sizes.insert(leaf_domain_sizes.end(),
        std::get<indices>(domains).get_leaf_domain_sizes().begin(),
        std::get<indices>(domains).get_leaf_domain_sizes().end());
    leaf_domain_sizes.insert(leaf_domain_sizes.end(),
        dmn_list_1.get_leaf_domain_sizes().begin(),
        dmn_list_1.get_leaf_domain_sizes().end());
    leaf_domain_sizes.insert(leaf_domain_sizes.end(),
        dmn_list_2.get_leaf_domain_sizes().begin(),
        dmn_list_2.get_leaf_domain_sizes().end());

    leaf_domain_steps.resize(leaf_domain_sizes.size(), 1);
    for (size_t i = 0; i < leaf_domain_sizes.size(); i++)
        for (size_t j = 0; j < i; j++)
            leaf_domain_steps[i] *= leaf_domain_sizes[j];
*/
    size = 1;

    for (int i = 0; i < 3; i++)
        size *= branch_domain_sizes[i];

    dmn_size() = size;
}

template<typename... domain_list>
void dmn_variadic<domain_list...>::reset()
{
    domain::reset();
    std::index_sequence_for<domain_list...> indices;
    std::get<indices>(domains).reset();

/*
    branch_domain_sizes.push_back(dmn_list_0.get_size());
    branch_domain_sizes.push_back(dmn_list_1.get_size());
    branch_domain_sizes.push_back(dmn_list_2.get_size());

    branch_domain_steps.resize(branch_domain_sizes.size(), 1);
    for (size_t i = 0; i < branch_domain_sizes.size(); i++)
        for (size_t j = 0; j < i; j++)
            branch_domain_steps[i] *= branch_domain_sizes[j];

    leaf_domain_sizes.insert(leaf_domain_sizes.end(),
        dmn_list_0.get_leaf_domain_sizes().begin(),
        dmn_list_0.get_leaf_domain_sizes().end());
    leaf_domain_sizes.insert(leaf_domain_sizes.end(),
        dmn_list_1.get_leaf_domain_sizes().begin(),
        dmn_list_1.get_leaf_domain_sizes().end());
    leaf_domain_sizes.insert(leaf_domain_sizes.end(),
        dmn_list_2.get_leaf_domain_sizes().begin(),
        dmn_list_2.get_leaf_domain_sizes().end());

    leaf_domain_steps.resize(leaf_domain_sizes.size(), 1);
    for (size_t i = 0; i < leaf_domain_sizes.size(); i++)
        for (size_t j = 0; j < i; j++)
            leaf_domain_steps[i] *= leaf_domain_sizes[j];
*/
    size = 1;

    for (int i = 0; i < 3; i++)
        size *= branch_domain_sizes[i];

    dmn_size() = size;
}

template<typename... domain_list>
template <typename ...Args>
int dmn_variadic<domain_list...>::operator()(Args&&... args) {
    static_assert( sizeof...(Args) < sizeof...(domain_list), "not enough args");
    return index_lookup(
        std::integral_constant<bool, (sizeof...(Args) == sizeof...(domain_list))>::type,
        std::forward<Args>(args)...
    );
}

template<typename... domain_list>
template<typename ...Args>
int dmn_variadic<domain_list...>::index_lookup(
    std::integral_constant<bool, true>::type,
    Args&&... args)
{
    std::cout << "Branch overload " << std::endl;
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

template<typename... domain_list>
template<typename ...Args>
int dmn_variadic<domain_list...>::index_lookup(
    std::integral_constant<bool, false>::type,
    Args&&... args)
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

