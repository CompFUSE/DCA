//-*-C++-*-

#ifndef TYPELIST_H_
#define TYPELIST_H_

template<typename parameters>
class dmn_0;
template<typename ... domain_list>
class dmn_variadic;

//----------------------------------------------------------------------------
template<typename T1, typename T2>
struct assert_same
{
    assert_same() {
        static_assert(std::is_same<T1, T2>::value, "Types must be equal");
    }
    static_assert(std::is_same<T1, T2>::value, "Types must be equal");
};
//----------------------------------------------------------------------------
class NullType {
};

//----------------------------------------------------------------------------
// utility function : takes a variadic list or arguments and returns nothing.
// It is used to convert a variadic pack expansion into a function so that
// arbitrary functions can be called with a pack expansion and drop the results.
//----------------------------------------------------------------------------
template<typename ...Ts>
void ignore_returnvalues(Ts&&...) {
}

//----------------------------------------------------------------------------
// From Simple C++11 metaprogramming, Peter Dimov
// http://pdimov.com/cpp2/simple_cxx11_metaprogramming.html
// interesting discussion of type transformations and examples
// mp_xxxx metafunctions taken from this source
//
// See also "TypeLists and a TypeList Toolbox via Variadic Templates"
// http://www.codeproject.com/Articles/1077852/TypeLists-and-a-TypeList-Toolbox-via-Variadic-Temp
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
/// generic mp_list type list
template<class ... T>
struct mp_list {
};

/// create an alias for backward naming compatibility with DCA+ typelists
template<class ... Ts>
using Typelist = mp_list<Ts...>;

//----------------------------------------------------------------------------
/// mp_rename A<...> -> B<...>
template<class A, template<class ...> class B>
struct mp_rename_impl;

template<template<class ...> class A, class... T, template<class...> class B>
struct mp_rename_impl<A<T...>, B>
{
    using type = B<T...>;
};

template<class A, template<class ...> class B>
using mp_rename = typename mp_rename_impl<A, B>::type;

//----------------------------------------------------------------------------
/// mp_size
template<class L>
struct mp_size_impl;

template<template<class ...> class L, class... T>
struct mp_size_impl<L<T...>>
{
    using type = std::integral_constant<std::size_t, sizeof...(T)>;
};

template<class L>
using mp_size = typename mp_size_impl<L>::type;

//----------------------------------------------------------------------------
/// mp_plus
template<class ... T>
struct mp_plus_impl;

template<class ... T>
using mp_plus = typename mp_plus_impl<T...>::type;

template<>
struct mp_plus_impl<>
{
    using type = std::integral_constant<int, 0>;
};

template<class T1, class ... T>
struct mp_plus_impl<T1, T...>
{
    static constexpr auto _v = T1::value + mp_plus<T...>::value;

    using type = std::integral_constant<
    typename std::remove_const<decltype(_v)>::type, _v>;
};
//----------------------------------------------------------------------------
/// mp_count
template<class L, class V> struct mp_count_impl;

template<template<class ...> class L, class... T, class V>
struct mp_count_impl<L<T...>, V>
{
    using type = mp_plus<std::is_same<T, V>...>;
};

template<class L, class V>
using mp_count = typename mp_count_impl<L, V>::type;

//----------------------------------------------------------------------------
/// mp_append
template<typename T, typename ...TL>
struct mp_append;

template<typename ...Ts>
struct mp_append<mp_list<Ts...>>
{
    typedef mp_list<Ts...> type;
};

template<typename T, typename ...Ts>
struct mp_append<mp_list<Ts...>, T>
{
    typedef mp_list<Ts..., T> type;
};

template<typename ...Ts1, typename ...Ts2>
struct mp_append<mp_list<Ts1...>, mp_list<Ts2...>>
{
    typedef mp_list<Ts1..., Ts2...> type;
};

template<typename ...Ts1,
    typename ...Ts2,
    typename ...Ts>
struct mp_append<mp_list<Ts1...>, mp_list<Ts2...>, Ts...>
{
    typedef typename mp_append<mp_list<Ts1..., Ts2...>, Ts...>::type type;
};

//----------------------------------------------------------------------------
/// mp_prepend
template<typename T, typename ...TL>
struct mp_prepend;

template<typename T, typename ...Ts>
struct mp_prepend<mp_list<Ts...>, T>
{
    typedef mp_list<T, Ts...> type;
};

template<typename ...Ts1, typename ...Ts2>
struct mp_prepend<mp_list<Ts1...>, mp_list<Ts2...>>
{
    typedef mp_list<Ts2..., Ts1...> type;
};

//----------------------------------------------------------------------------
// From : True Story: Efficient Packing
// http://talesofcpp.fusionfenix.com/post-22/true-story-efficient-packing
// get Nth type from typelist/tuple efficiently without recursion
//----------------------------------------------------------------------------
template<std::size_t I, typename T>
struct _indexed {
    using type = T;
};

template<typename Is, typename ...Ts>
struct _indexer;

template<std::size_t ...Is, typename ...Ts>
struct _indexer<std::index_sequence<Is...>, Ts...> : _indexed<Is, Ts> ...
{
};

template<std::size_t I, typename ...Ts>
struct _at_index {
    template<typename T>
    static _indexed<I, T> _select(_indexed<I, T>);

    using _impl = _indexer<std::index_sequence_for<Ts...>, Ts...>;
    using type = typename decltype(_select(_impl {}))::type;
};

// was tuple_element, but we want typelist_element
template<std::size_t I, typename Tuple>
struct mp_element;

template<std::size_t I, typename ...Ts>
struct mp_element<I, Typelist<Ts...>> : _at_index<I, Ts...>
{
};

//-----------------------------------------------------------------------------
/// Search a typelist for a first occurrence of the type T

// Implementation: has index as a template parameter
template<size_t idx, typename T, class List>
struct mp_index_of_impl;

template<size_t idx, typename T> /// The type T is not in the list
struct mp_index_of_impl<idx, T, mp_list<>>
{
    using type = std::integral_constant<int, -1>;
};

template<size_t idx, typename T, typename ... Ts>    ///> The type is found
struct mp_index_of_impl<idx, T, mp_list<T, Ts...>>
{
    using type = std::integral_constant<int, idx>;
};

template<size_t idx, typename T, typename H, typename ... Ts>  ///> Recursion
struct mp_index_of_impl<idx, T, mp_list<H, Ts...>>
{
    using type = typename mp_index_of_impl<idx + 1, T, mp_list<Ts...>>::type;
};

// Wrapping to supply initial index 0
template<typename T, class List>
struct mp_index_of {
    static constexpr int value = -1;
};

// Specializing for idx >= 0
template<typename T, typename ... Ts>
struct mp_index_of<T, mp_list<Ts...>>
{
    using type = typename mp_index_of_impl<0, T, mp_list<Ts...>>::type;
    using value_type = typename type::value_type;
    static constexpr value_type value = type::value;
};

// Specializing for idx >= 0
template<typename T, typename ... Ts>
struct mp_index_of<mp_list<Ts...>, T>
{
//  static_assert(false, "Parameter ordering incorrect");
};

//-----------------------------------------------------------------------------
/// swap first element out of a typelist
template<typename TList, typename T1, typename T2>
struct mp_swap {
};

template<typename T1, typename T2>
struct mp_swap<mp_list<>, T1, T2> {
    typedef mp_list<> type;
};

template<typename T1, typename T2, typename ...Ts>
struct mp_swap<mp_list<T1, Ts...>, T1, T2> {
    typedef mp_list<T2, Ts...> type;
};

template<typename T0, typename ...Ts, typename T1, typename T2>
struct mp_swap<mp_list<T0, Ts...>, T1, T2> {
    typedef typename mp_prepend<typename mp_swap<mp_list<Ts...>, T1, T2>::type, T0>::type type;
};

//----------------------------------------------------------------------------
// namespace TL
// reimplement DCA+ typelist operations using C++11 alias types
// on the new typelist operations mp_list, mp_rename, mp_plus, mp_count &etc
//----------------------------------------------------------------------------
namespace TL
{

    /***********************************
     * Length
     ***********************************/
    template<typename L>
    using Length = mp_size<L>;

    /***********************************
     *	 NumberOf
     ***********************************/
    template<class L, class V>
    using NumberOf = mp_count<L, V>;

    /***********************************
     * TypeAt
     ***********************************/
    template<int I, typename Ts>
    using TypeAt = mp_element<I, Ts>;

    /***********************************
     * 	IndexOf
     ***********************************/
    template<typename T, typename Ts>
    using IndexOf = mp_index_of<T,Ts>;

    /**********************************
     * 	Append
     **********************************/
    template<typename T1, typename T2>
    using Append = mp_append<T1,T2>;

    /**********************************
     * Prepend
     **********************************/
    template<typename T1, typename T2>
    using Prepend = mp_prepend<T1,T2>;

    /***********************************
     * 	Swap
     ***********************************/
    template<typename T1, typename T2, typename T3>
    using Swap = mp_swap<T1,T2,T3>;

    /***********************************
     * 	Print
     ***********************************/
    template<class Head>
    struct print_type {
        static void print() {
            std::cout << "\t" << __PRETTY_FUNCTION__ << std::endl;
            //std::cout << "\t" << Head::get_name() << std::endl;
        }

        static void print(std::ostream& ss) {
            //ss << "\t" << __PRETTY_FUNCTION__ << std::endl;
            ss << "\t" << Head::get_name() << std::endl;
        }

        template<class stream_type>
        static void to_JSON(stream_type& ss) {
            ss << __PRETTY_FUNCTION__;
            //ss << "\t" << Head::get_name() << std::endl;
        }
    };

    //----------------------------------------------------------------------------
    // PrintTL
    //----------------------------------------------------------------------------
    // basic print displays the type of the template instantiation
    template<typename D>
    struct printTL {
        static void print() {
            print(std::cout);
        }
        //
        static void print(std::ostream &stream) {
            stream << "\t" << __PRETTY_FUNCTION__ << "\n";
        }
        static void to_JSON(std::ostream &stream) {
            stream << "\t" << __PRETTY_FUNCTION__ << "\n";
        }
    };

    // dmn_0 override is actually the same as basic, but provided
    // for future customization
    template<typename Domain>
    struct printTL<dmn_0<Domain>> {
        static void print()
        {
            print(std::cout);
        }
        static void print(std::ostream &stream)
            {
            stream << "\t" << __PRETTY_FUNCTION__ << "\n";
        }
    };

    // dmn_variadic prints out all subdomains recursively via pack expansion
    template<typename ... Domains>
    struct printTL<dmn_variadic<Domains...>> {
        static void print() {
            print(std::cout);
        }
        // we can't expand a pack out without passing it as a parameter
        // so expand the pack as a parameter list, and drop dummy return values
        // use func(),0 because func() returns void
        static void print(std::ostream &s) {
            ignore_returnvalues((printTL<Domains>::print(std::cout),0)...);
        }
    };

    // dmn_variadic prints out all subdomains recursively via pack expansion
    template<typename Domain, typename ... Domains>
    struct printTL<Typelist<Domain, Domains...>> {
        static void print() {
            print(std::cout);
        }
        // we can't expand a pack out without passing it as a parameter
        // so expand the pack as a parameter list, and drop dummy return values
        // use func(),0 because func() returns void
        static void print(std::ostream &s) {
            ignore_returnvalues((printTL<Domains>::print(s),0)...);
        }

        static void to_JSON(std::ostream &s) {
            s << "\"";
            print_type<Domain>::to_JSON(s);
            if (sizeof...(Domains)==0) {
            s << "\"\n";
        }
        else {
            s << "\",\n";
            printTL<Typelist<Domains...>>::to_JSON(s);
        }
    }
};
}

#endif
