//-*-C++-*-

#include "type_list_definitions.h"



#ifndef DOMAIN_TYPE_OPERATIONS_H_
#define DOMAIN_TYPE_OPERATIONS_H_

namespace TL
{
/*
**********************************
*	SwapNext
**********************************
*/
  /*
  template <class T0, class T1,class T2>
  struct SwapNext<dmn_0<T0>, T1, T2> {
    typedef dmn_0<T0> Result;
  };

  template <class T1,class T2>
  struct SwapNext<dmn_0<T1>, T1, T2> {
    typedef dmn_0<T2> Result;
  };
  
  template <class D0, class T1, class T2>
  struct SwapNext<dmn_1<D0>,T1,T2> {
    typedef dmn_1<typename Swap<D0,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class T1, class T2>
    struct SwapNext<dmn_2<D0, D1>,T1,T2> {
    typedef dmn_2<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class D2, class T1, class T2>
    struct SwapNext<dmn_3<D0, D1, D2>,T1,T2> {
    typedef dmn_3<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result,
		  typename Swap<D2,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class D2, class D3, class T1, class T2>
  struct SwapNext<dmn_4<D0, D1, D2, D3>,T1,T2> {
    typedef dmn_4<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result,
		  typename Swap<D2,T1,T2>::Result,
		  typename Swap<D3,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class D2, class D3, class D4, class T1, class T2>
  struct SwapNext<dmn_5<D0, D1, D2, D3, D4>,T1,T2> {
    typedef dmn_5<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result,
		  typename Swap<D2,T1,T2>::Result,
		  typename Swap<D3,T1,T2>::Result,
		  typename Swap<D4,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class D2, class D3, class D4, class D5, class T1, class T2>
  struct SwapNext<dmn_6<D0, D1, D2, D3, D4, D5>,T1,T2> {
    typedef dmn_6<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result,
		  typename Swap<D2,T1,T2>::Result,
		  typename Swap<D3,T1,T2>::Result,
		  typename Swap<D4,T1,T2>::Result,
		  typename Swap<D5,T1,T2>::Result> Result;
  };
  */

/*
**********************************
*	Swap
**********************************
*/
  
  template <class T0, class T1,class T2>
  struct Swap<dmn_0<T0>, T1, T2> {
    typedef dmn_0<T0> Result;
  };

  template <class T1,class T2>
  struct Swap<dmn_0<T1>, T1, T2> {
    typedef dmn_0<T2> Result;
  };

  template <class T0, class T1,class T2>
  struct Swap<dmn_0<T0>, dmn_0<T1>, dmn_0<T2> > {
    typedef dmn_0<T0> Result;
  };

  template <class T1,class T2>
  struct Swap<dmn_0<T1>, dmn_0<T1>, dmn_0<T2> > {
    typedef dmn_0<T2> Result;
  };

  template <class D0, class T1, class T2>
  struct Swap<dmn_1<D0>,T1,T2> {
    typedef dmn_1<typename Swap<D0,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class T1, class T2>
    struct Swap<dmn_2<D0, D1>,T1,T2> {
    typedef dmn_2<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class D2, class T1, class T2>
    struct Swap<dmn_3<D0, D1, D2>,T1,T2> {
    typedef dmn_3<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result,
		  typename Swap<D2,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class D2, class D3, class T1, class T2>
  struct Swap<dmn_4<D0, D1, D2, D3>,T1,T2> {
    typedef dmn_4<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result,
		  typename Swap<D2,T1,T2>::Result,
		  typename Swap<D3,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class D2, class D3, class D4, class T1, class T2>
  struct Swap<dmn_5<D0, D1, D2, D3, D4>,T1,T2> {
    typedef dmn_5<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result,
		  typename Swap<D2,T1,T2>::Result,
		  typename Swap<D3,T1,T2>::Result,
		  typename Swap<D4,T1,T2>::Result> Result;
  };

  template <class D0, class D1, class D2, class D3, class D4, class D5, class T1, class T2>
  struct Swap<dmn_6<D0, D1, D2, D3, D4, D5>,T1,T2> {
    typedef dmn_6<typename Swap<D0,T1,T2>::Result,
                  typename Swap<D1,T1,T2>::Result,
		  typename Swap<D2,T1,T2>::Result,
		  typename Swap<D3,T1,T2>::Result,
		  typename Swap<D4,T1,T2>::Result,
		  typename Swap<D5,T1,T2>::Result> Result;
  };


/*
**********************************
*	Print
**********************************
*/

  template <class T0>
  struct printTL<dmn_0<T0> > {
    static void print() {
      std::cout <<  "\t" << /*typeid(T0).name()*/ __PRETTY_FUNCTION__ << "\n";
    }
  };

  template <class D0>
  struct printTL<dmn_1<D0> > {
    static void print() {
      printTL<D0>::print();
    }
  };

  template <class D0, class D1>
    struct printTL<dmn_2<D0, D1> > {
    static void print() {
      printTL<D0>::print();
      printTL<D1>::print();
    }
  };

  template <class D0, class D1, class D2>
    struct printTL<dmn_3<D0, D1, D2> > {
    static void print() {
      printTL<D0>::print();
      printTL<D1>::print();
      printTL<D2>::print();
    }
  };

  template <class D0, class D1, class D2, class D3>
    struct printTL<dmn_4<D0, D1, D2, D3> > {
    static void print() {
      printTL<D0>::print();
      printTL<D1>::print();
      printTL<D2>::print();
      printTL<D3>::print();
    }
  };

  template <class D0, class D1, class D2, class D3, class D4>
  struct printTL<dmn_5<D0, D1, D2, D3, D4> > {
    static void print() {
      printTL<D0>::print();
      printTL<D1>::print();
      printTL<D2>::print();
      printTL<D3>::print();
      printTL<D4>::print();
    }
  };

  template <class D0, class D1, class D2, class D3, class D4, class D5>
  struct printTL<dmn_6<D0, D1, D2, D3, D4, D5> > {
    static void print() {
      printTL<D0>::print();
      printTL<D1>::print();
      printTL<D2>::print();
      printTL<D3>::print();
      printTL<D4>::print();
      printTL<D5>::print();
    }
  };


}

#endif
