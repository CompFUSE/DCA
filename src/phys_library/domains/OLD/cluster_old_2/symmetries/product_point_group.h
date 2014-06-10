//-*-C++-*-

/*
 *      Author: Peter Staar
 */


#ifndef POINT_GROUP_PRODUCT_H_
#define POINT_GROUP_PRODUCT_H_

template<typename point_group_1, typename point_group_2>
class point_group_product
{};

template<typename point_group_1, typename point_group_2>
class point_group_product_left_2_right
{};

/*
 *   product_group_action< Typelist, Typelist >
 */

template<typename head_1, typename head_2, typename tail_1, typename tail_2>
class point_group_product<Typelist<head_1, tail_1>, Typelist<head_2, tail_2> >
{
public:

  typedef typename point_group_product_left_2_right<Typelist<head_1, tail_1>, Typelist<head_2, tail_2> >::Result product_ij;
  typedef typename point_group_product_left_2_right<Typelist<head_2, tail_2>, Typelist<head_1, tail_1> >::Result product_ji;

  typedef typename Append<product_ij, product_ji>::Result Result;
};

template<typename head_1, typename head_2, typename tail_1, typename tail_2>
class point_group_product_left_2_right<Typelist<head_1, tail_1>, Typelist<head_2, tail_2> >
{
public:

  typedef typename point_group_product_left_2_right <head_1, Typelist<head_2, tail_2> >::Result product_0j;
  typedef typename point_group_product_left_2_right <tail_1, Typelist<head_2, tail_2> >::Result product_ij;

  typedef typename Append<product_0j, product_ij>::Result Result;
};

template<typename head_1, typename head_2, typename tail_2>
class point_group_product_left_2_right<Typelist<head_2, tail_2>, TYPELIST_1(head_1)>
{
public:
  typedef typename point_group_product_left_2_right<Typelist<head_2, tail_2>, head_1>::Result Result;
};

template<typename head_1, typename head_2, typename tail_2>
class point_group_product_left_2_right<TYPELIST_1(head_1), Typelist<head_2, tail_2> >
{
public:
  typedef typename point_group_product_left_2_right<head_1, Typelist<head_2, tail_2> >::Result Result;
};

template<typename head_1, typename head_2>
class point_group_product_left_2_right<TYPELIST_1(head_1), TYPELIST_1(head_2)>
{
public:
  typedef product_group_action<head_1, head_2> product_t; 
  typedef TYPELIST_1(product_t)                Result;
};




/*
 *   product_group_action< T, Typelist >
 */

template<typename head, typename head_2, typename tail_2>
class point_group_product_left_2_right<head, Typelist<head_2, tail_2> >
{
public:
  typedef typename point_group_product_left_2_right<head, tail_2>::Result new_tail;

  typedef typename Swap<Typelist<head_2, new_tail>, head_2, product_group_action<head, head_2> >::Result Result;
};

template<typename head, typename last_head>
class point_group_product_left_2_right<head, TYPELIST_1(last_head) >
{
public:
  typedef product_group_action<head, last_head> product_t;
  typedef TYPELIST_1(product_t)                 Result;
};

/*
 *   product_group_action< Typelist, T >
 */

template<typename head, typename head_1, typename tail_1>
class point_group_product_left_2_right<Typelist<head_1, tail_1>, head>
{
public:
  typedef typename point_group_product_left_2_right<tail_1, head>::Result new_tail;
  typedef typename Swap<Typelist<head_1, new_tail>, head_1, product_group_action<head_1, head> >::Result Result;
};

template<typename head, typename last_head>
class point_group_product_left_2_right<TYPELIST_1(last_head), head>
{
public:
  typedef product_group_action<last_head, head> product_t;
  typedef TYPELIST_1(product_t)                 Result;
};
#endif
