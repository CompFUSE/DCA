//-*-C++-*-

#ifndef SQUARE_2D_H
#define SQUARE_2D_H

/*!
 *  \author: Peter Staar
 */

/*!
 *  group-actions
 */

typedef Cn_2D<1,4> Cn_2D_1_4_type;
typedef Cn_2D<2,4> Cn_2D_2_4_type;
typedef Cn_2D<3,4> Cn_2D_3_4_type;
typedef Cn_2D<4,4> Cn_2D_4_4_type;

typedef Sn_2D<0,8> Sn_2D_0_8_type;
typedef Sn_2D<1,8> Sn_2D_1_8_type;
typedef Sn_2D<2,8> Sn_2D_2_8_type;
typedef Sn_2D<3,8> Sn_2D_3_8_type;
typedef Sn_2D<4,8> Sn_2D_4_8_type;
typedef Sn_2D<5,8> Sn_2D_5_8_type;
typedef Sn_2D<6,8> Sn_2D_6_8_type;
typedef Sn_2D<7,8> Sn_2D_7_8_type;
typedef Sn_2D<8,8> Sn_2D_8_8_type;

/*
 *  pointgroup :: set of group-actions
 */

struct C4
{
  typedef Typelist<Cn_2D_1_4_type,
                     Cn_2D_2_4_type,
                     Cn_2D_3_4_type,
                     Cn_2D_4_4_type> point_group_type_list;
};

struct S4
{
  typedef Typelist<Sn_2D_0_8_type,
                     Sn_2D_2_8_type,
                     Sn_2D_4_8_type,
                     Sn_2D_6_8_type> point_group_type_list;
};

struct S4_plus
{
  typedef Typelist<identity_group_operation<2>,
                     Cn_2D_2_4_type,
                     Sn_2D_0_8_type,
                     Sn_2D_2_8_type> point_group_type_list;
};

struct S8
{
  typedef Typelist<Sn_2D_0_8_type,
                     Sn_2D_1_8_type,
                     Sn_2D_2_8_type,
                     Sn_2D_3_8_type,
                     Sn_2D_4_8_type,
                     Sn_2D_5_8_type,
                     Sn_2D_6_8_type,
                     Sn_2D_7_8_type,
                     identity_group_operation<2>> point_group_type_list;
};

struct D4_trial
{
  typedef point_group_product< C4::point_group_type_list,
                               S8::point_group_type_list>::Result point_group_type_list;
};

struct D4
{
  typedef Typelist<Cn_2D_1_4_type,
                      Cn_2D_2_4_type,
                      Cn_2D_3_4_type,
                      Sn_2D_0_8_type,
                      Sn_2D_1_8_type,
                      Sn_2D_2_8_type,
                      Sn_2D_3_8_type,
                      Sn_2D_4_8_type,
                      Sn_2D_5_8_type,
                      Sn_2D_6_8_type,
                      Sn_2D_7_8_type,
                      identity_group_operation<2>> point_group_type_list;
};

#endif
