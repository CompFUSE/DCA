//-*-C++-*-

/*
 *      Author: peterstaar
 */


#ifndef ROT_2_H_
#define ROT_2_H_

class Rot_2 : public point_group<2>
{
public:

  const static int size = 1;

  typedef Cn_2D<2,4> Cn_2D_2_4_type;

  typedef TYPELIST_1(Cn_2D_2_4_type) symmetry_list;
};


#endif
