//-*-C++-*-

#ifndef TYPELIST_H_
#define TYPELIST_H_

//----------------------------------------------------------------------------
template <typename T1, typename T2>
struct assert_same
{
    assert_same() {
        static_assert(std::is_same<T1, T2>::value, "Types must be equal");
    }
    static_assert(std::is_same<T1, T2>::value, "Types must be equal");
};
//----------------------------------------------------------------------------
class NullType {};

//----------------------------------------------------------------------------
// utility function : takes a variadic list or arguments and returns nothing.
// It is used to convert a variadic pack expansion into a function so that
// arbitrary functions can be called with a pack expansion and drop the results.
//----------------------------------------------------------------------------
template<typename...Ts>
void ignore_returnvalues(Ts&&...) {}

//----------------------------------------------------------------------------
// From Simple C++11 metaprogramming, Peter Dimov
// http://pdimov.com/cpp2/simple_cxx11_metaprogramming.html
// interesting discussion of type transformations and examples
// mp_xxxx metafunctions taken from this source
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// generic mp_list type list
template<class... T>
struct mp_list {};

/// create an alias for backward naming compatibility with DCA+ typelists
template<class... Ts>
using Typelist = mp_list<Ts...>;

//----------------------------------------------------------------------------
// mp_rename A<...> -> B<...>
template<class A, template<class...> class B>
struct mp_rename_impl;

template<template<class...> class A, class... T, template<class...> class B>
struct mp_rename_impl<A<T...>, B>
{
    using type = B<T...>;
};

template<class A, template<class...> class B>
using mp_rename = typename mp_rename_impl<A, B>::type;

//----------------------------------------------------------------------------
// mp_size
template<class L>
struct mp_size_impl;

template<template<class...> class L, class... T>
struct mp_size_impl<L<T...>>
{
    using type = std::integral_constant<std::size_t, sizeof...(T)>;
};

template<class L>
using mp_size = typename mp_size_impl<L>::type;

//----------------------------------------------------------------------------
// mp_plus
template<class... T> struct mp_plus_impl;

template<class... T> using mp_plus = typename mp_plus_impl<T...>::type;

template<> struct mp_plus_impl<>
{
    using type = std::integral_constant<int, 0>;
};

template<class T1, class... T> struct mp_plus_impl<T1, T...>
{
    static constexpr auto _v = T1::value + mp_plus<T...>::value;

    using type = std::integral_constant<
        typename std::remove_const<decltype(_v)>::type, _v>;
};
//----------------------------------------------------------------------------
// mp_count
template<class L, class V> struct mp_count_impl;

template<template<class...> class L, class... T, class V>
struct mp_count_impl<L<T...>, V>
{
    using type = mp_plus<std::is_same<T, V>...>;
};

template<class L, class V>
using mp_count = typename mp_count_impl<L, V>::type;

//----------------------------------------------------------------------------
// From : True Story: Efficient Packing
// http://talesofcpp.fusionfenix.com/post-22/true-story-efficient-packing
// get Nth type from typelist/tuple efficiently without recursion
//----------------------------------------------------------------------------
template <std::size_t I, typename T>
struct _indexed {
  using type = T;
};

template <typename Is, typename ...Ts>
struct _indexer;

template <std::size_t ...Is, typename ...Ts>
struct _indexer<std::index_sequence<Is...>, Ts...> : _indexed<Is, Ts>...
{};

template <std::size_t I, typename ...Ts>
struct _at_index {
  template <typename T>
  static _indexed<I, T> _select(_indexed<I, T>);

  using _impl = _indexer<std::index_sequence_for<Ts...>, Ts...>;
  using type = typename decltype(_select(_impl{}))::type;
};

// was tuple_element, but we want typelist_element
template <std::size_t I, typename Tuple>
struct typelist_element;

template <std::size_t I, typename ...Ts>
struct typelist_element<I, Typelist<Ts...>> : _at_index<I, Ts...>
{};

//----------------------------------------------------------------------------
// namespace TL
// reimplement DCA+ typelist operations using C++11 alias types
// on the new typelist operations mp_list, mp_rename, mp_plus, mp_count &etc
//----------------------------------------------------------------------------
namespace TL
{

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
    template <int I, typename ...Ts>
    using TypeAt = typelist_element<I, Ts...>;

    /***********************************
     * 	IndexOf
     ***********************************/

    template <class TList, class T> struct IndexOf;

    template <class T>
    struct IndexOf<NullType, T>
    {
        enum { value = -1 };
    };

    template <class T, class Tail>
    struct IndexOf<Typelist<T, Tail>, T>
    {
        enum { value = 0 };
    };

    template <class Head, class Tail, class T>
    struct IndexOf<Typelist<Head, Tail>, T>
    {
    private:
        enum { temp = IndexOf<Tail, T>::value };
    public:
        enum { value = temp == -1 ? -1 : 1 + temp };
    };

    /***********************************
     * 	IndexOf_At
     ***********************************/

    template <class TList, class T, int N>
    struct IndexOf_At;

    template <class T, int N>
    struct IndexOf_At<NullType, T, N>
    {
        enum { value = -1 };
    };

    template <class Tail, class T>
    struct IndexOf_At<Typelist<T, Tail>, T, 0>
    {
    private:
        enum { temp = 0 };
    public:
        enum { value = temp == -1 ? -1 : 1 + temp };
    };

    template <class Tail, class T, int N>
    struct IndexOf_At<Typelist<T, Tail>, T, N>
    {
    private:
        enum { temp = IndexOf_At<Tail, T, N-1>::value };
    public:
        enum { value = temp == -1 ? -1 : 1 + temp };
    };

    template <class Head, class Tail, class T, int N>
    struct IndexOf_At<Typelist<Head, Tail>, T, N>
    {
    private:
        enum { temp = IndexOf_At<Tail, T, N>::value };
    public:
        enum { value = temp == -1 ? -1 : 1 + temp };
    };



    /**********************************
     * 	Append
     **********************************/
/*
    template <class TList, class T> struct Append;

    template <>
    struct Append<NullType, NullType>
    {
        typedef NullType Result;
    };

    template <class T>
    struct Append<NullType, T>
    {
        typedef TYPELIST_1(T) Result;
    };

    template <class Head, class Tail>
    struct Append<NullType, Typelist<Head, Tail> >
    {
        typedef Typelist<Head, Tail> Result;
    };

    template <class Head, class Tail, class T>
    struct Append<Typelist<Head, Tail>, T>
    {
        typedef Typelist<Head, typename Append<Tail, T>::Result> Result;
    };
*/
    /***********************************
     * 	Swap
     ***********************************/
/*
    template <class TList,class T1,class T2>
    struct Swap;

    template <class T1,class T2>
    struct Swap<NullType, T1, T2> {
        typedef NullType Result;
    };

    template <class T1, class T2,class Tail>
    struct Swap<Typelist<T1,Tail>,T1, T2> {
        typedef Typelist<T2,Tail> Result;
    };

    template <class Head, class Tail, class T1, class T2>
    struct Swap<Typelist<Head,Tail>,T1,T2> {
        typedef Typelist<Head, typename Swap<Tail,T1,T2>::Result> Result;
    };
*/

    /***********************************
     * 	Print
     ***********************************/
/*
    template <typename Head>
    struct print_typename {
        static void print() {
            //std::cout << "\t" << __PRETTY_FUNCTION__ << std::endl;
            std::cout << "\t" << Head::get_name() << std::endl;
        }
    };

    template <class Head>
    struct print_type {
        static void print() {
            std::cout << "\t" << __PRETTY_FUNCTION__ << std::endl;
            std::cout << "\t" << Head::get_name() << std::endl;
        }

        static void print(std::ostream &s) {
            s << "\t" << __PRETTY_FUNCTION__ << std::endl;
            //std::cout << "\t" << Head::get_name() << std::endl;
        }

        static void print(std::stringstream& ss) {
            //ss << "\t" << __PRETTY_FUNCTION__ << std::endl;
            ss << "\t" << Head::get_name() << std::endl;
        }

        template<class stream_type>
        static void to_JSON(stream_type& ss) {
            ss << __PRETTY_FUNCTION__;
            //ss << "\t" << Head::get_name() << std::endl;
        }
    };

    template <class Ts>
    struct printTL {
        static void print(std::ostream &s) {
//            print_type<typename TL::TypeAt<0,Ts>::type>::print(s);
//            printTL<typename TL::Tail>::print(s);
        }

        static void print(std::stringstream& ss) {
//            print_type<typename TL::Head>::print(ss);
//            printTL<typename TL::Tail>::print(ss);
        }

        template<class stream_type>
        static void to_JSON(stream_type& ss) {
            ss << "\"";
//            print_type<typename TL::Head>::to_JSON(ss);

//            if(Length<typename TL::Tail>::value != 0)
//                ss << "\",\n";

//            printTL<typename TL::Tail>::to_JSON(ss);
        }
    };

    template <>
    struct printTL<NullType> {
        static void print() {
            std::cout << "[end]" << std::endl;
        }

        static void print(std::stringstream& ss) {
            ss << std::endl;
        }

        template<class stream_type>
        static void to_JSON(stream_type& ss) {
            ss << "\"";
        }
    };
*/

}


#endif
