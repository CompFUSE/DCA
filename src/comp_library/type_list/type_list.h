//-*-C++-*-

#include "type_list_definitions.h"

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
class NullType {
};

struct EmptyType {};

template <class T, class U>
struct Typelist
{
  typedef T Head;
  typedef U Tail;
};

namespace TL
{
  /***********************************
   *	Length
   ***********************************/

  template <class TList> struct Length;

  template <> struct Length<NullType>
  {
    enum { value = 0 };
  };

  template <class T, class U>
  struct Length< Typelist<T, U> >
  {
    enum { value = 1 + Length<U>::value };
  };

  /***********************************
   *	NumberOf
   ***********************************/

  template <class TList, class T>
  struct NumberOf;

  template <class T>
  struct NumberOf<NullType,T>
  {
    enum { value = 0 };
  };

  template <class Tail, class T>
  struct NumberOf< Typelist<T,Tail> ,T >
  {
    enum { value = 1+NumberOf<Tail,T>::value };
  };

  template <class Head, class Tail, class T>
  struct NumberOf< Typelist<Head,Tail> ,T >
  {
    enum { value = NumberOf<Tail,T>::value };
  };

  /***********************************
   *	NumberOftypes
   ***********************************/

  template<class TList, class T, int SIZE>
  struct NUMBER_OF_TYPES;

  template<class T, int SIZE>
  struct NUMBER_OF_TYPES<NullType, T, SIZE>
  {
    enum { value = 0 };
  };

  template<class TList, class T>
  struct NUMBER_OF_TYPES<TList, T, 0>
  {
    enum { value = 0 };
  };

  template<class Tail, class T, int SIZE>
  struct NUMBER_OF_TYPES<Typelist<T, Tail>, T, SIZE>
  {
    enum { value = 1+NUMBER_OF_TYPES<Tail, T, SIZE-1>::value };
  };

  template <class Head, class Tail, class T, int SIZE>
  struct NUMBER_OF_TYPES<Typelist<Head,Tail>, T, SIZE>
  {
    enum { value = NUMBER_OF_TYPES<Tail, T, SIZE-1>::value };
  };

  /*
**********************************
*	TypeListAt
**********************************
*/

  template <class TList, unsigned int index>
  struct TypeListAt;

  template <class Head, class Tail>
  struct TypeListAt<Typelist<Head, Tail>, 0>
  {
    typedef Tail Result;
  };

  template <class Head, class Tail, unsigned int i>
  struct TypeListAt<Typelist<Head, Tail>, i>
  {
    typedef typename TypeListAt<Tail, i - 1>::Result Result;
  };


  /*
**********************************
*	TypeAt
**********************************
*/

  template <class TList, unsigned int index> struct TypeAt;

  template <class Head, class Tail>
  struct TypeAt<Typelist<Head, Tail>, 0>
  {
    typedef Head Result;
  };

  template <class Head, class Tail, unsigned int i>
  struct TypeAt<Typelist<Head, Tail>, i>
  {
    typedef typename TypeAt<Tail, i - 1>::Result Result;
  };


  /*
**********************************
*	IndexOf
**********************************
*/

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

  /*
**********************************
*	IndexOf_At
**********************************
*/


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
*	Append
**********************************/


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

  /*
**********************************
*	Swap
**********************************
*/

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

  /*
**********************************
*	SwapAll
**********************************
*/

  template <class TList,class T1,class T2>
  struct SwapAll;

  template <class T1,class T2>
  struct SwapAll<NullType, T1, T2> {
    typedef NullType Result;
  };

  template <class T1, class T2,class Tail>
  struct SwapAll<Typelist<T1,Tail>,T1, T2> {
    typedef Typelist<T2, typename SwapAll<Tail,T1,T2>::Result> Result;
  };

  template <class Head, class Tail, class T1, class T2>
  struct SwapAll<Typelist<Head,Tail>,T1,T2> {
    typedef Typelist<Head, typename SwapAll<Tail,T1,T2>::Result> Result;
  };

  /*
**********************************
*	Erase
**********************************
*/

  template <class TList,class T>
  struct Erase;

  template <class T>
  struct Erase<NullType, T> {
    typedef NullType Result;
  };

  template <class T,class Tail>
  struct Erase<Typelist<T,Tail>,T> {
    typedef Tail Result;
  };

  template <class Head,class Tail, class T>
  struct Erase<Typelist<Head,Tail>,T> {
    typedef Typelist<Head, typename Erase<Tail,T>::Result> Result;
  };

  /*
**********************************
*	Erase All
**********************************
*/

  template <class TList, class T>
  struct EraseAll;

  template <class T>
  struct EraseAll<NullType,T> {
    typedef NullType R;
  };

  template <class T, class Tail>
  struct EraseAll<Typelist<T,Tail>,T> {
    typedef typename EraseAll<Tail,T>::R R;
  };

  template <class Head, class Tail, class T>
  struct EraseAll<Typelist<Head,Tail>,T> {
    typedef Typelist<Head,typename EraseAll<Tail,T>::R> R;
  };

  /*
**********************************
*	Erase Duplicates
**********************************
*/

  template <class TList>
  struct NoDuplicates;

  template <>
  struct NoDuplicates<NullType> {
    typedef NullType R;
  };

  template <class Head, class Tail>
  struct NoDuplicates<Typelist<Head,Tail> > {
    typedef typename NoDuplicates<Tail>::R L1;
    typedef typename Erase<L1,Head>::R L2;
  public:
    typedef Typelist<Head,L2> R;
  };


  /*
**********************************
*	Is_Typelist
**********************************
*/

  template <typename T>
  struct IsTypelist {
    enum { r = false };
  };

  template <typename H, class T>
  struct IsTypelist< Typelist<H,T> > {
    enum { r = true };
  };

  template <>
  struct IsTypelist<NullType> {
    enum { r = true };
  };

  /*
**********************************
*	Print
**********************************
*/

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

  template <class TL>
  struct printTL {
    static void print() {
      print_type<typename TL::Head>::print();
      printTL<typename TL::Tail>::print();
    }

    static void print(std::stringstream& ss) {
      print_type<typename TL::Head>::print(ss);
      printTL<typename TL::Tail>::print(ss);
    }

    template<class stream_type>
    static void to_JSON(stream_type& ss) {
      ss << "\"";
      print_type<typename TL::Head>::to_JSON(ss);

      if(Length<typename TL::Tail>::value != 0)
	ss << "\",\n";

      printTL<typename TL::Tail>::to_JSON(ss);
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


}


#endif
